"""
Fuzzy matcher module - Python replacement for nrgrep_coords C binary.

Uses the Python 'regex' module for approximate/fuzzy matching.
"""
import regex
from typing import List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class MatchResult:
    """Result of a fuzzy match."""
    start: int      # 0-indexed start position in file
    end: int        # 0-indexed end position in file (exclusive)
    matched_text: str


class FuzzyMatcher:
    """
    Fuzzy pattern matcher using Python's regex module.

    Replaces the nrgrep_coords C binary for approximate string matching.
    """

    def __init__(
        self,
        max_errors: int = 0,
        allow_insertions: bool = True,
        allow_deletions: bool = True,
        allow_substitutions: bool = True,
        case_insensitive: bool = True
    ):
        """
        Initialize the fuzzy matcher.

        Args:
            max_errors: Maximum number of errors (insertions + deletions + substitutions)
            allow_insertions: Allow insertion errors
            allow_deletions: Allow deletion errors
            allow_substitutions: Allow substitution errors
            case_insensitive: Perform case-insensitive matching
        """
        self.max_errors = max_errors
        self.allow_insertions = allow_insertions
        self.allow_deletions = allow_deletions
        self.allow_substitutions = allow_substitutions
        self.case_insensitive = case_insensitive

    def _build_fuzzy_pattern(self, pattern: str) -> str:
        """
        Build a fuzzy regex pattern with error specifications.

        The regex module supports fuzzy matching with the syntax:
        (?b)  - bestmatch flag
        {e<=N} - max N errors total
        {i<=N} - max N insertions
        {d<=N} - max N deletions
        {s<=N} - max N substitutions
        """
        if self.max_errors == 0:
            return pattern

        # Build error specification
        error_spec = []

        if self.allow_insertions and self.allow_deletions and self.allow_substitutions:
            # All error types allowed - use combined error limit
            error_spec.append(f'e<={self.max_errors}')
        else:
            # Specific error types
            if self.allow_insertions:
                error_spec.append(f'i<={self.max_errors}')
            if self.allow_deletions:
                error_spec.append(f'd<={self.max_errors}')
            if self.allow_substitutions:
                error_spec.append(f's<={self.max_errors}')

        if error_spec:
            error_str = ','.join(error_spec)
            # Use bestmatch to prefer matches with fewer errors
            return f'(?b)({pattern}){{{error_str}}}'

        return pattern

    def search(self, pattern: str, text: str, byte_offset: int = 0) -> List[MatchResult]:
        """
        Search for pattern matches in text.

        Args:
            pattern: The regex pattern to search for
            text: The text to search in
            byte_offset: Byte offset to add to match positions

        Returns:
            List of MatchResult objects
        """
        flags = regex.IGNORECASE if self.case_insensitive else 0

        fuzzy_pattern = self._build_fuzzy_pattern(pattern)

        results = []
        try:
            for match in regex.finditer(fuzzy_pattern, text, flags=flags, overlapped=True):
                results.append(MatchResult(
                    start=match.start() + byte_offset,
                    end=match.end() + byte_offset,
                    matched_text=match.group()
                ))
        except regex.error as e:
            # If fuzzy matching fails, try exact matching
            try:
                for match in regex.finditer(pattern, text, flags=flags, overlapped=True):
                    results.append(MatchResult(
                        start=match.start() + byte_offset,
                        end=match.end() + byte_offset,
                        matched_text=match.group()
                    ))
            except regex.error:
                pass

        return results


def search_file(
    pattern: str,
    file_path: str,
    max_errors: int = 0,
    error_types: str = 'ids',
    case_insensitive: bool = True,
    max_buffer_size: int = 1600000
) -> str:
    """
    Search a file for pattern matches, returning output in nrgrep format.

    This function mimics the output format of the nrgrep_coords binary:
    [start, end]: matched_text

    Args:
        pattern: The regex pattern to search for
        file_path: Path to the file to search
        max_errors: Maximum number of errors
        error_types: String containing error types to allow (i=insertions, d=deletions, s=substitutions)
        case_insensitive: Perform case-insensitive matching
        max_buffer_size: Maximum buffer size for reading file

    Returns:
        Output string in nrgrep format
    """
    allow_insertions = 'i' in error_types
    allow_deletions = 'd' in error_types
    allow_substitutions = 's' in error_types

    matcher = FuzzyMatcher(
        max_errors=max_errors,
        allow_insertions=allow_insertions,
        allow_deletions=allow_deletions,
        allow_substitutions=allow_substitutions,
        case_insensitive=case_insensitive
    )

    output_lines = []

    # Read file in chunks if it's large
    with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
        # For FASTA files, read the entire content
        content = f.read()

    # Perform search
    results = matcher.search(pattern, content)

    for result in results:
        # Output format matches nrgrep_coords: [start, end]: matched_text
        # nrgrep uses 1-indexed positions
        output_lines.append(f'[{result.start + 1}, {result.end}]: {result.matched_text}')

    return '\n'.join(output_lines)


def parse_mismatch_option(option: str) -> Tuple[int, str]:
    """
    Parse a mismatch option string like "2ids" into components.

    Args:
        option: Mismatch option string (e.g., "2ids", "1s", "0ids")

    Returns:
        Tuple of (max_errors, error_types)
    """
    if not option:
        return 0, 'ids'

    # Extract leading number
    num_str = ''
    error_types = ''

    for char in option:
        if char.isdigit():
            num_str += char
        else:
            error_types += char

    max_errors = int(num_str) if num_str else 0
    error_types = error_types if error_types else 'ids'

    return max_errors, error_types


def run_search(
    pattern: str,
    file_path: str,
    mismatch_option: str = '0ids',
    case_insensitive: bool = True
) -> str:
    """
    Run a search with the same interface as the nrgrep_coords binary.

    Args:
        pattern: The regex pattern to search for
        file_path: Path to the file to search
        mismatch_option: Mismatch option string (e.g., "2ids")
        case_insensitive: Perform case-insensitive matching

    Returns:
        Output string in nrgrep format
    """
    max_errors, error_types = parse_mismatch_option(mismatch_option)

    return search_file(
        pattern=pattern,
        file_path=file_path,
        max_errors=max_errors,
        error_types=error_types,
        case_insensitive=case_insensitive
    )
