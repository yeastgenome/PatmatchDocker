"""
FASTA utilities module - Python replacement for generate_sequence_index.pl

Provides functions for parsing FASTA files and generating sequence indices.
"""
from typing import Iterator, Tuple, Dict, List, Optional
import os


def generate_sequence_index(fasta_file: str) -> Tuple[List[int], Dict[int, str]]:
    """
    Generate a sequence index for a FASTA file.

    This is a Python port of generate_sequence_index.pl.
    Returns byte offsets for each sequence in the file.

    Args:
        fasta_file: Path to the FASTA file

    Returns:
        Tuple of (offset_list, offset_to_name_dict)
        - offset_list: List of byte offsets where sequences start
        - offset_to_name_dict: Dict mapping byte offset to sequence name
    """
    record_offset_list = []
    seq_name_for_offset = {}

    bytecount = 0

    with open(fasta_file, 'rb') as f:
        for line in f:
            # Decode line for pattern matching
            line_str = line.decode('utf-8', errors='replace')

            if line_str.startswith('>'):
                # Extract sequence name (first word after >)
                parts = line_str[1:].split()
                if parts:
                    seq_name = parts[0]

                    # Record header offset
                    record_offset_list.append(bytecount)
                    seq_name_for_offset[bytecount] = '>' + seq_name

                    # Record sequence start offset (right after header line)
                    seq_start = bytecount + len(line)
                    record_offset_list.append(seq_start)
                    seq_name_for_offset[seq_start] = seq_name

            bytecount += len(line)

    return record_offset_list, seq_name_for_offset


def get_name_offset(offset: int, record_offset_list: List[int]) -> int:
    """
    Perform binary search to find the sequence name offset for a given byte position.

    Args:
        offset: The byte offset to look up
        record_offset_list: Sorted list of record offsets

    Returns:
        The sequence start offset for the given position
    """
    if not record_offset_list:
        return 0

    low = 0
    high = len(record_offset_list) - 1

    while high > low:
        middle = (low + high) // 2

        if record_offset_list[middle] == offset:
            return offset
        elif high - low == 1:
            if offset >= record_offset_list[high]:
                return record_offset_list[high]
            else:
                return record_offset_list[low]
        elif record_offset_list[middle] < offset:
            low = middle
        else:  # record_offset_list[middle] > offset
            high = middle - 1

    return record_offset_list[low]


def parse_fasta(fasta_file: str) -> Iterator[Tuple[str, str, str]]:
    """
    Parse a FASTA file and yield sequences one at a time.

    Args:
        fasta_file: Path to the FASTA file

    Yields:
        Tuples of (sequence_name, defline, sequence)
    """
    current_name = None
    current_defline = None
    current_seq = []

    with open(fasta_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Yield previous sequence if exists
                if current_name is not None:
                    yield current_name, current_defline, ''.join(current_seq)

                # Start new sequence
                current_defline = line
                parts = line[1:].split()
                current_name = parts[0] if parts else ''
                current_seq = []
            else:
                current_seq.append(line)

    # Yield last sequence
    if current_name is not None:
        yield current_name, current_defline, ''.join(current_seq)


def get_sequence(dataset: str, seqname: str, data_dir: str = '/data/patmatch/') -> Dict[str, str]:
    """
    Get a specific sequence from a FASTA file.

    Args:
        dataset: Dataset name (with or without .seq extension)
        seqname: Name of the sequence to retrieve
        data_dir: Base directory for data files

    Returns:
        Dict with 'defline' and 'seq' keys
    """
    if '.seq' not in dataset:
        dataset = dataset + '.seq'
    if not dataset.startswith(data_dir) and 'patmatch' not in dataset:
        dataset = os.path.join(data_dir, dataset)

    defline = ""
    seq = ""
    found = False

    with open(dataset, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith('>' + seqname.lower()):
                found = True
                defline = line
                continue
            if not found:
                continue
            if found and line.startswith('>'):
                break
            seq += line

    defline = defline.replace('"', "'")

    return {'defline': defline, 'seq': seq}


def set_seq_length(data_file: str) -> Tuple[Dict[str, int], Dict[str, bool]]:
    """
    Calculate sequence lengths from a FASTA file.

    Args:
        data_file: Path to the FASTA file

    Returns:
        Tuple of (seq_name_to_length, seq_has_stop)
    """
    seq_name_to_length = {}
    seq_has_stop = {}

    current_name = None
    current_seq = []

    with open(data_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                # Process previous sequence
                if current_name is not None:
                    seq = ''.join(current_seq)
                    canon_name = current_name.rstrip(',')
                    has_stop = seq.endswith('*')
                    seq_has_stop[canon_name] = has_stop
                    if has_stop:
                        seq = seq[:-1]
                    seq_name_to_length[canon_name] = len(seq)

                # Start new sequence
                parts = line.replace('>', '').split(' ')
                current_name = parts[0].rstrip(',')
                current_seq = []
            else:
                current_seq.append(line.strip())

    # Process last sequence
    if current_name is not None:
        seq = ''.join(current_seq)
        canon_name = current_name.rstrip(',')
        has_stop = seq.endswith('*')
        seq_has_stop[canon_name] = has_stop
        if has_stop:
            seq = seq[:-1]
        seq_name_to_length[canon_name] = len(seq)

    return seq_name_to_length, seq_has_stop


def load_locus_data(data_dir: str = '/data/patmatch/') -> Dict[str, Tuple[str, str, str]]:
    """
    Load locus data from locus.txt file.

    Args:
        data_dir: Base directory for data files

    Returns:
        Dict mapping sequence name to (gene_name, sgdid, description)
    """
    locus_file = os.path.join(data_dir, 'locus.txt')
    name_to_data = {}

    if not os.path.exists(locus_file):
        return name_to_data

    with open(locus_file, 'r', encoding='utf-8') as f:
        for line in f:
            pieces = line.strip().split('\t')
            if len(pieces) >= 3:
                seq_name = pieces[0]
                gene_name = pieces[1]
                sgdid = pieces[2]
                desc = pieces[3] if len(pieces) > 3 else ''
                name_to_data[seq_name] = (gene_name, sgdid, desc)

    return name_to_data


def load_not_feature_data(data_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Load chromosome and ORF data from a NotFeature FASTA file.

    Args:
        data_file: Path to the NotFeature FASTA file

    Returns:
        Tuple of (seq_name_to_chr, seq_name_to_orfs)
    """
    seq_name_to_chr = {}
    seq_name_to_orfs = {}

    with open(data_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                # Parse header like:
                # >A:2170-2479, Chr I from 2170-2479, Genome Release 64-3-1, between YAL068C and YAL067W-A
                pieces = line.strip().replace('>', '').split(' ')
                seq_name = pieces[0].replace(',', '')

                if len(pieces) > 2:
                    chr_name = pieces[2]
                    seq_name_to_chr[seq_name] = chr_name

                if 'between ' in line:
                    orfs = line.strip().split('between ')[1]
                    orfs = orfs.replace('and', '-')
                    seq_name_to_orfs[seq_name] = orfs

    return seq_name_to_chr, seq_name_to_orfs
