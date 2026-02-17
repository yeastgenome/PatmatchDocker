"""
Pattern converter module - Python port of patmatch_to_nrgrep.pl

Converts Patmatch pattern expressions to regex patterns.

Class options:
- nucleotide: assume nucleotide pattern
- peptide: assume protein pattern
- complement: reverse complement nucleotide
"""
import re
from typing import Tuple

INFINITE = -1


class PatternConverter:
    """Converts Patmatch patterns to regex patterns."""

    # IUPAC codes for nucleotides
    NUCLEOTIDE_IUPAC = {
        'R': '[AG]',
        'Y': '[CT]',
        'S': '[GC]',
        'W': '[AT]',
        'M': '[AC]',
        'K': '[GT]',
        'V': '[ACG]',
        'H': '[ACT]',
        'D': '[AGT]',
        'B': '[CGT]',
    }

    # IUPAC codes for peptides
    PEPTIDE_IUPAC = {
        'J': '[IFVLWMAGCY]',  # Hydrophobic
        'O': '[TSHEDQNKR]',   # Hydrophilic
        'B': '[DN]',          # Aspartic acid or Asparagine
        'Z': '[EQ]',          # Glutamic acid or Glutamine
    }

    # Complement mapping for nucleotides
    COMPLEMENT_MAP = str.maketrans('ATCGRYSWMKVHDB', 'TAGCYRSWKMBDHV')

    def __init__(self, sequence_type: str = 'nucleotide'):
        """
        Initialize the converter.

        Args:
            sequence_type: One of 'nucleotide', 'peptide', or 'complement'
        """
        self.sequence_type = sequence_type.lower()
        if self.sequence_type not in ('nucleotide', 'peptide', 'complement'):
            raise ValueError(f"Invalid sequence type: {sequence_type}")

    def convert(self, pattern: str) -> str:
        """
        Convert a Patmatch pattern to a regex pattern.

        Args:
            pattern: The Patmatch pattern to convert

        Returns:
            The converted regex pattern
        """
        pattern = self._prepare_pattern(pattern)
        pattern = self._fix_wildcards(pattern)
        pattern = self._fix_repetitions(pattern)
        pattern = self._substitute_iupac_codes(pattern)
        pattern = self._finalize_pattern(pattern)
        return pattern

    def _prepare_pattern(self, pattern: str) -> str:
        """
        Prepare the pattern by removing spaces and converting to uppercase.
        Get reverse complement if necessary.
        """
        pattern = re.sub(r'\s', '', pattern)
        pattern = pattern.upper()
        if self.sequence_type == 'complement':
            pattern = self._get_reverse_complement(pattern)
        return pattern

    def _fix_wildcards(self, pattern: str) -> str:
        """Substitute Patmatch wildcards with regex wildcards."""
        if self.sequence_type == 'peptide':
            pattern = pattern.replace('X', '.')
        else:  # nucleotide or complement
            pattern = pattern.replace('N', '.').replace('X', '.')
        return pattern

    def _fix_repetitions(self, pattern: str) -> str:
        """
        Convert Patmatch repetitions {m}, {m,}, {,m}, {m,n} to regex format.
        """
        if '{' not in pattern:
            return pattern

        result = []
        chars = list(pattern)
        i = 0

        while i < len(chars):
            if chars[i] == '}':
                result = self._process_repetition(result)
                i += 1
            else:
                result.append(chars[i])
                i += 1

        return ''.join(result)

    def _process_repetition(self, chars: list) -> list:
        """Process a repetition specification found in the pattern."""
        # Extract repetition information (everything after '{')
        rep_info = self._extract_repetition_info(chars)
        # Extract the pattern to repeat (could be single char, [...], or (...))
        repeat_pattern = self._extract_repeat_pattern(chars)
        # Build and append the new repeat structure
        lower, upper = self._parse_repeat_info(rep_info)
        new_repeat = self._build_repeat(lower, upper, repeat_pattern)
        chars.append(new_repeat)
        return chars

    def _extract_repetition_info(self, chars: list) -> str:
        """Extract repetition info from chars, removing the {info part."""
        rep_info = []
        while chars:
            char = chars.pop()
            if char == '{':
                break
            rep_info.insert(0, char)
        return ''.join(rep_info)

    def _extract_repeat_pattern(self, chars: list) -> str:
        """Extract the pattern to be repeated from chars."""
        if not chars:
            return ''

        char = chars.pop()

        if char in (')', ']'):
            # It's a grouped pattern - need to find the matching opener
            bracket_stack = [char]
            left_bracket = '(' if char == ')' else '['
            right_bracket = char
            repeat = [char]

            while bracket_stack and chars:
                char = chars.pop()
                repeat.insert(0, char)
                if char == right_bracket:
                    bracket_stack.append(char)
                elif char == left_bracket:
                    bracket_stack.pop()

            return ''.join(repeat)
        else:
            return char

    def _parse_repeat_info(self, repeat_info: str) -> Tuple[int, int]:
        """Parse repeat info to get lower and upper bounds."""
        lower = 0
        upper = 0

        if re.match(r'^,\d+$', repeat_info):
            # {,m} - 0 to m
            upper = int(repeat_info.split(',')[1])
        elif re.match(r'^\d+,$', repeat_info):
            # {m,} - m to infinity
            lower = int(repeat_info.split(',')[0])
            upper = INFINITE
        elif re.match(r'^\d+$', repeat_info):
            # {m} - exactly m
            lower = int(repeat_info)
            upper = int(repeat_info)
        elif re.match(r'^\d+,\d+$', repeat_info):
            # {m,n} - m to n
            parts = repeat_info.split(',')
            lower = int(parts[0])
            upper = int(parts[1])

        return lower, upper

    def _build_repeat(self, lower: int, upper: int, pattern: str) -> str:
        """Build the regex repeat pattern."""
        result = []

        # Add the pattern 'lower' times (mandatory occurrences)
        for _ in range(lower):
            result.append(pattern)

        # Add optional occurrences
        if upper == INFINITE:
            result.append(f'{pattern}*')
        else:
            num_optional = upper - lower
            for _ in range(num_optional):
                result.append(f'{pattern}?')

        return ''.join(result)

    def _substitute_iupac_codes(self, pattern: str) -> str:
        """Substitute IUPAC wildcard characters with character classes."""
        if self.sequence_type == 'peptide':
            for code, replacement in self.PEPTIDE_IUPAC.items():
                pattern = pattern.replace(code, replacement)
        else:  # nucleotide or complement
            for code, replacement in self.NUCLEOTIDE_IUPAC.items():
                pattern = pattern.replace(code, replacement)

        pattern = self._remove_nested_brackets(pattern)
        return pattern

    def _remove_nested_brackets(self, pattern: str) -> str:
        """
        Remove nested brackets and duplicate characters within brackets.
        """
        result = []
        bracket_stack = []
        char_set = set()

        for char in pattern:
            if char == '[':
                if not bracket_stack:
                    result.append(char)
                    char_set = set()
                bracket_stack.append(char)
            elif char == ']':
                bracket_stack.pop()
                if not bracket_stack:
                    result.append(char)
            else:
                if not bracket_stack:
                    result.append(char)
                else:
                    if char not in char_set:
                        result.append(char)
                        char_set.add(char)

        return ''.join(result)

    def _finalize_pattern(self, pattern: str) -> str:
        """
        Place the entire pattern in parentheses.
        Convert anchors < and > to regex anchors ^ and $.
        """
        has_start_anchor = pattern.startswith('<')
        has_end_anchor = pattern.endswith('>')

        if has_start_anchor:
            pattern = pattern[1:]
        if has_end_anchor:
            pattern = pattern[:-1]

        if has_start_anchor and has_end_anchor:
            pattern = f'^({pattern})$'
        elif has_start_anchor:
            pattern = f'^({pattern})'
        elif has_end_anchor:
            pattern = f'({pattern})$'
        else:
            pattern = f'({pattern})'

        return pattern

    def _get_reverse_complement(self, pattern: str) -> str:
        """Get the reverse complement of a pattern."""
        pattern = self._complement_nucleotides(pattern)
        pattern = self._reverse_pattern(pattern)
        return pattern

    def _complement_nucleotides(self, pattern: str) -> str:
        """Get the complement of a nucleotide pattern."""
        pattern = pattern.translate(self.COMPLEMENT_MAP)

        # Swap anchors for complement
        if pattern.startswith('<'):
            pattern = '>' + pattern[1:]
        elif pattern.startswith('>'):
            pattern = '<' + pattern[1:]

        if pattern.endswith('>'):
            pattern = pattern[:-1] + '<'
        elif pattern.endswith('<'):
            pattern = pattern[:-1] + '>'

        return pattern

    def _reverse_pattern(self, pattern: str) -> str:
        """Reverse a pattern, keeping grouped elements together."""
        chars = list(pattern)
        result = []

        while chars:
            char = chars.pop()
            if char in (')', ']', '}'):
                group = self._extract_group(char, chars)
                result.append(group)
            else:
                result.append(char)

        return ''.join(result)

    def _extract_group(self, closer: str, chars: list) -> str:
        """Extract a group from the character list."""
        opener_map = {')': '(', ']': '[', '}': '{'}
        opener = opener_map[closer]

        group = [closer]
        internal_chars = []

        while chars:
            char = chars.pop()
            if char == opener:
                if opener != '{':
                    # For () and [], reverse internal chars
                    group.insert(0, ''.join(internal_chars))
                    group.insert(0, char)
                    break
                else:
                    # For {}, this is a repetition specifier
                    group.insert(0, char)
                    # Get the character/group being repeated
                    if chars:
                        repeater = chars.pop()
                        if repeater in (')', ']'):
                            repeater_group = self._extract_group(repeater, chars)
                            group.insert(0, repeater_group)
                        else:
                            group.insert(0, repeater)
                    break
            elif char in (')', ']', '}'):
                nested_group = self._extract_group(char, chars)
                internal_chars.append(nested_group)
            else:
                if closer == '}':
                    group.insert(0, char)
                else:
                    internal_chars.append(char)

        return ''.join(group)


def convert_pattern(pattern: str, sequence_type: str = 'nucleotide') -> str:
    """
    Convenience function to convert a Patmatch pattern.

    Args:
        pattern: The Patmatch pattern
        sequence_type: One of 'nucleotide', 'peptide', or 'complement'

    Returns:
        The converted regex pattern
    """
    converter = PatternConverter(sequence_type)
    return converter.convert(pattern)


def get_pattern_type_from_option(option: str) -> str:
    """
    Convert command-line option to sequence type.

    Args:
        option: One of '-n' (nucleotide), '-p' (peptide), '-c' (complement)

    Returns:
        The sequence type string
    """
    option_map = {
        '-n': 'nucleotide',
        '-p': 'peptide',
        '-c': 'complement',
    }
    return option_map.get(option, 'nucleotide')
