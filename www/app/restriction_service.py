"""
Restriction enzyme service - Python replacement for scan_for_matches C binary.

Provides restriction site search functionality.
"""
import os
import re
import regex
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class EnzymeInfo:
    """Information about a restriction enzyme."""
    name: str
    offset: int
    pattern: str
    overhang: int


@dataclass
class CutSite:
    """A cut site found in a sequence."""
    start: int
    end: int
    strand: str  # 'watson' or 'crick'


class RestrictionService:
    """Service for finding restriction enzyme cut sites."""

    # IUPAC ambiguity codes for DNA
    IUPAC_MAP = {
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
        'N': '[ACGT]',
    }

    def __init__(self, data_dir: str = '/data/restriction_mapper/'):
        """
        Initialize the restriction service.

        Args:
            data_dir: Directory containing enzyme data files
        """
        self.data_dir = data_dir
        self.fasta_file = os.path.join(data_dir, 'orf_genomic.seq')

    def _convert_pattern_to_regex(self, pattern: str) -> str:
        """Convert enzyme pattern with IUPAC codes to regex."""
        result = pattern.upper()
        for code, replacement in self.IUPAC_MAP.items():
            result = result.replace(code, replacement)
        return result

    def _get_complement(self, seq: str) -> str:
        """Get the complement of a DNA sequence."""
        complement_map = str.maketrans('ATCGRYSWMKVHDB', 'TAGCYRSWKMBDHV')
        return seq.upper().translate(complement_map)

    def _get_reverse_complement(self, seq: str) -> str:
        """Get the reverse complement of a DNA sequence."""
        return self._get_complement(seq)[::-1]

    def load_enzymes(self, enzyme_type: Optional[str] = None) -> List[EnzymeInfo]:
        """
        Load enzyme definitions from file.

        Args:
            enzyme_type: Type of enzymes to load (None for all)

        Returns:
            List of EnzymeInfo objects
        """
        enzyme_file = self._get_enzyme_file(enzyme_type)

        enzymes = []
        with open(enzyme_file, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    enzymes.append(EnzymeInfo(
                        name=parts[0],
                        offset=int(parts[1]),
                        pattern=parts[2],
                        overhang=int(parts[3])
                    ))

        return enzymes

    def _get_enzyme_file(self, enzyme_type: Optional[str]) -> str:
        """Get the appropriate enzyme file based on type."""
        if enzyme_type is None:
            return os.path.join(self.data_dir, 'rest_enzymes')

        if 'Six-base' in enzyme_type:
            return os.path.join(self.data_dir, 'rest_enzymes.6base')
        if 'blunt' in enzyme_type:
            return os.path.join(self.data_dir, 'rest_enzymes.blunt')
        if "3'" in enzyme_type or '3' in enzyme_type:
            return os.path.join(self.data_dir, 'rest_enzymes.3')
        if "5'" in enzyme_type or '5' in enzyme_type:
            return os.path.join(self.data_dir, 'rest_enzymes.5')

        return os.path.join(self.data_dir, 'rest_enzymes')

    def search_pattern(self, pattern: str, sequence: str, search_complement: bool = True) -> List[CutSite]:
        """
        Search for a pattern in a sequence.

        Args:
            pattern: The enzyme recognition pattern
            sequence: The DNA sequence to search
            search_complement: Whether to also search complement strand

        Returns:
            List of CutSite objects
        """
        regex_pattern = self._convert_pattern_to_regex(pattern)
        sites = []

        # Search Watson strand
        for match in regex.finditer(regex_pattern, sequence.upper(), flags=regex.IGNORECASE):
            sites.append(CutSite(
                start=match.start() + 1,  # 1-indexed
                end=match.end(),
                strand='watson'
            ))

        # Search Crick strand (reverse complement)
        if search_complement:
            rc_pattern = self._convert_pattern_to_regex(self._get_reverse_complement(pattern))
            for match in regex.finditer(rc_pattern, sequence.upper(), flags=regex.IGNORECASE):
                sites.append(CutSite(
                    start=match.end(),  # Reversed positions for complement
                    end=match.start() + 1,
                    strand='crick'
                ))

        return sites

    def get_sequence_from_fasta(self, name: str) -> Tuple[str, str]:
        """
        Get a sequence from the FASTA file.

        Args:
            name: Sequence name to look for

        Returns:
            Tuple of (defline, sequence)
        """
        name = name.replace('SGD:', '')

        defline = ""
        seq = ""

        with open(self.fasta_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    parts = line.split()
                    # Check various name formats
                    check_names = [
                        parts[0].replace('>', '').lower(),
                        parts[1].lower() if len(parts) > 1 else '',
                        parts[2].replace('SGDID:', '').replace(',', '').lower() if len(parts) > 2 else ''
                    ]

                    if name.lower() in check_names:
                        defline = line
                    continue
                elif defline:
                    seq = line
                if seq:
                    break

        defline = defline.replace('"', "'")
        return defline, seq

    def write_seq_file(self, defline: str, seq: str, seq_file: str) -> Tuple[str, str, int]:
        """
        Write a sequence to a temporary file.

        Args:
            defline: FASTA defline
            seq: Sequence string
            seq_file: Output file path

        Returns:
            Tuple of (seq_name, chr_coords, seq_length)
        """
        # Remove non-alphabet characters
        seq = re.sub(r'[^a-zA-Z]', '', seq)

        with open(seq_file, 'w') as f:
            f.write(defline + '\n')
            f.write(seq + '\n')

        seq_name = "Unnamed"
        chr_coords = ""

        # Parse defline for SGD format
        # >YAL067C SEO1 SGDID:S000000062, Chr I from 9016-7235, Genome Release ...
        if 'SGDID:' in defline and 'Genome Release' in defline:
            parts = defline.replace('>', '').split()
            systematic_name = parts[0]
            gene_name = parts[1] if len(parts) > 1 else ''

            comma_parts = defline.split(', ')
            if len(comma_parts) > 1:
                chr_coords = comma_parts[1]

            seq_name = systematic_name
            if gene_name:
                seq_name = f"{gene_name}/{systematic_name}"

        return seq_name, chr_coords, len(seq)

    def do_search(self, enzyme_type: Optional[str], sequence: str) -> Dict[str, List[CutSite]]:
        """
        Search for all enzyme cut sites in a sequence.

        Args:
            enzyme_type: Type of enzymes to search for
            sequence: The DNA sequence to search

        Returns:
            Dict mapping enzyme name to list of cut sites
        """
        enzymes = self.load_enzymes(enzyme_type)
        results = {}

        for enzyme in enzymes:
            sites = self.search_pattern(enzyme.pattern, sequence)
            if sites:
                results[enzyme.name] = sites

        return results


def get_enzyme_types(data_dir: str, enzyme_type: str) -> Dict[str, str]:
    """
    Load enzyme type classifications.

    Args:
        data_dir: Directory containing enzyme files
        enzyme_type: Type to load

    Returns:
        Dict mapping enzyme name to type string
    """
    enzyme_hash = {}

    enzyme_file = 'rest_enzymes.blunt'
    if '3' in enzyme_type:
        enzyme_file = 'rest_enzymes.3'
    elif '5' in enzyme_type:
        enzyme_file = 'rest_enzymes.5'

    file_path = os.path.join(data_dir, enzyme_file)
    if os.path.exists(file_path):
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    enzyme_hash[parts[0]] = enzyme_type

    return enzyme_hash


def run_restriction_search(
    sequence: str,
    enzyme_type: str,
    data_dir: str = '/data/restriction_mapper/'
) -> Tuple[Dict[str, Dict], List[str], str]:
    """
    Run a complete restriction site search.

    Args:
        sequence: The DNA sequence to search
        enzyme_type: Type of enzymes to search for
        data_dir: Directory containing enzyme data

    Returns:
        Tuple of (data_dict, not_cut_enzymes, error_message)
    """
    service = RestrictionService(data_dir)
    seq_len = len(sequence)

    # Load enzymes and search
    enzymes = service.load_enzymes(enzyme_type)
    data_hash = {}
    offset_map = {}
    overhang_map = {}
    recognition_seq_map = {}
    not_cut_enzymes = []

    for enzyme in enzymes:
        sites = service.search_pattern(enzyme.pattern, sequence)
        offset_map[enzyme.name] = str(enzyme.offset)
        overhang_map[enzyme.name] = str(enzyme.overhang)
        recognition_seq_map[enzyme.name] = enzyme.pattern

        if sites:
            coords = ':'.join([f"{s.start},{s.end}" for s in sites])
            data_hash[enzyme.name] = coords
        elif enzyme_type.lower() in ('all', '') or enzyme_type.lower().startswith('enzymes that do not'):
            if enzyme.name not in not_cut_enzymes:
                not_cut_enzymes.append(enzyme.name)

    not_cut_enzymes.sort()

    if enzyme_type.startswith('enzymes that do not'):
        return {}, not_cut_enzymes, ''

    # Apply cut limit filters
    if 'cut' in enzyme_type:
        cut_limit = 2 if 'twice' in enzyme_type else 1

        new_data_hash = {}
        for enzyme_name, coords in data_hash.items():
            coord_pairs = coords.split(':')
            w_cut = 0
            c_cut = 0

            for coord_pair in coord_pairs:
                parts = coord_pair.split(',')
                if len(parts) == 2:
                    beg = int(parts[0])
                    end = int(parts[1])
                    if beg < end:
                        w_cut += 1
                    else:
                        c_cut += 1

            if (c_cut == cut_limit and w_cut <= cut_limit) or (w_cut == cut_limit and c_cut <= cut_limit):
                new_data_hash[enzyme_name] = coords

        data_hash = new_data_hash

    # Load enzyme type classifications
    enzyme_type_map = {}
    for et in ["3' overhang", "5' overhang", "blunt end"]:
        enzyme_type_map.update(get_enzyme_types(data_dir, et))

    # Build final result data
    data = {}
    for enzyme_name in sorted(data_hash.keys()):
        if ('overhang' in enzyme_type or 'blunt' in enzyme_type) and \
           enzyme_type_map.get(enzyme_name, '') != enzyme_type:
            continue

        cut_w = []
        cut_c = []
        cut_positions = data_hash[enzyme_name].split(':')
        cut_all = []

        for position in cut_positions:
            parts = position.split(',')
            if len(parts) != 2:
                continue

            beg = int(parts[0])
            end = int(parts[1])
            offset = int(offset_map.get(enzyme_name, 0))
            overhang = int(overhang_map.get(enzyme_name, 0))

            if beg < end:  # Watson strand
                cut_site = beg + offset - 1
                if cut_site not in cut_w:
                    cut_w.append(cut_site)
            else:  # Crick strand
                beg, end = end, beg
                cut_site = beg + offset + overhang - 1
                if cut_site not in cut_c:
                    cut_c.append(cut_site)

            if cut_site not in cut_all:
                cut_all.append(cut_site)

        cut_all.append(seq_len)

        # Calculate fragments
        pre_cut_site = 0
        found = set()
        cut_fragments = []

        for cut_site in sorted(cut_all):
            cut_size = cut_site - pre_cut_site
            if cut_size != 0 and cut_size not in found:
                cut_fragments.append(cut_size)
                found.add(cut_size)
            pre_cut_site = cut_site

        cut_site_w = ', '.join(str(x) for x in sorted(cut_w))
        cut_site_c = ', '.join(str(x) for x in sorted(cut_c))
        fragments_real = ', '.join(str(x) for x in cut_fragments)
        fragments_sorted = ', '.join(str(x) for x in sorted(cut_fragments, reverse=True))

        data[enzyme_name] = {
            'cut_site_on_watson_strand': cut_site_w,
            'cut_site_on_crick_strand': cut_site_c,
            'fragment_size': fragments_sorted,
            'fragment_size_in_real_order': fragments_real,
            'offset': offset_map.get(enzyme_name, ''),
            'overhang': overhang_map.get(enzyme_name, ''),
            'recognition_seq': recognition_seq_map.get(enzyme_name, ''),
            'enzyme_type': enzyme_type_map.get(enzyme_name, '')
        }

    return data, not_cut_enzymes, ''
