"""
Patmatch service module - Core pattern matching logic.

This is the main service that orchestrates pattern matching using the
pattern converter, fuzzy matcher, and FASTA utilities.
"""
import os
import re
import json
import hashlib
import time
import threading
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

from pattern_converter import PatternConverter
from fuzzy_matcher import run_search as fuzzy_search
from fasta_utils import (
    generate_sequence_index,
    get_name_offset,
    get_sequence,
    set_seq_length,
    load_locus_data,
    load_not_feature_data
)

# Try to import boto3 for S3 uploads
try:
    import boto3
    HAS_BOTO3 = True
except ImportError:
    boto3 = None
    HAS_BOTO3 = False


# Constants
MAX_BUFFER_SIZE = 1600000
MIN_TOKEN = 3
MINHITS = 500
MAXHITS = 100000
DEFAULT_MAXHITS = 500
DAY_IN_SECONDS = 86400

# Default directories (can be overridden by environment variables)
DATA_DIR = os.environ.get('DATA_DIR', '/data/patmatch/')
TMP_DIR = os.environ.get('TMP_DIR', '/var/www/tmp/')
CONF_DIR = os.environ.get('CONF_DIR', '/var/www/conf/')


def get_config(conf: Optional[str] = None) -> Dict:
    """
    Load configuration from JSON file.

    Args:
        conf: Configuration file name (without .json extension)

    Returns:
        Configuration dictionary
    """
    if conf is None:
        conf = 'patmatch'
    if not conf.endswith('.json'):
        conf = conf + '.json'

    config_path = os.path.join(CONF_DIR, conf)
    with open(config_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def clean_up_temp_files():
    """Remove temp files older than one day."""
    now = time.time()
    for f in os.listdir(TMP_DIR):
        file_path = os.path.join(TMP_DIR, f)
        if os.path.isfile(file_path) and os.stat(file_path).st_mtime < now - DAY_IN_SECONDS:
            try:
                os.remove(file_path)
            except OSError:
                pass


def upload_file_to_s3_async(file, filename: str):
    """Upload file to S3 asynchronously."""
    try:
        upload_file_to_s3(file, filename)
    except Exception as e:
        print(f"Error uploading file: {e}")
    finally:
        try:
            file.close()
        except Exception:
            pass


def upload_file_to_s3(file, filename: str) -> str:
    """
    Upload file to S3 if boto3 and S3_BUCKET are available.

    Returns:
        S3 URL or empty string
    """
    s3_bucket = os.environ.get('S3_BUCKET')

    if not HAS_BOTO3 or not s3_bucket:
        clean_up_temp_files()
        return ""

    s3_key = 'patmatch/' + filename

    s3 = boto3.client('s3')
    file.seek(0)
    s3.upload_fileobj(file, s3_bucket, s3_key, ExtraArgs={'ACL': 'public-read'})
    clean_up_temp_files()

    return f"https://{s3_bucket}.s3.amazonaws.com/{s3_key}"


def get_download_url(tmp_file: str) -> str:
    """
    Get download URL for a result file.

    Args:
        tmp_file: Temporary file name

    Returns:
        Download URL (S3 or local)
    """
    download_file = os.path.join(TMP_DIR, tmp_file)
    this_file = Path(download_file)

    if not this_file.exists():
        return ""

    # Calculate MD5 hash for unique filename
    with this_file.open('rb') as fh:
        md5sum = hashlib.md5(fh.read()).hexdigest()

    if md5sum:
        new_filename = md5sum + '.txt'
        new_file_path = os.path.join(TMP_DIR, new_filename)
        os.rename(download_file, new_file_path)
        tmp_file = new_filename

    file_path = os.path.join(TMP_DIR, tmp_file)
    f = open(file_path, 'rb')

    s3_bucket = os.environ.get('S3_BUCKET')
    if HAS_BOTO3 and s3_bucket:
        thread = threading.Thread(target=upload_file_to_s3_async, args=(f, tmp_file))
        thread.start()
        return f"https://{s3_bucket}.s3.amazonaws.com/patmatch/{tmp_file}"
    else:
        f.close()
        return ""


def cleanup_pattern(pattern: str) -> str:
    """Decode common URL escapes in pattern."""
    return (pattern
            .replace('%28', '(').replace('%29', ')')
            .replace('%7B', '{').replace('%7D', '}')
            .replace('%5B', '[').replace('%5D', ']')
            .replace('%2C', ',')
            .replace('%5E', '^'))


def check_pattern(pattern: str, seqtype: str) -> str:
    """
    Validate a pattern for the given sequence type.

    Returns:
        Error message or empty string if valid
    """
    if seqtype in ('pep', 'protein'):
        if 'u' in pattern.lower():
            return 'Invalid peptide character found in pattern.'
    else:
        invalid_chars = ('E', 'F', 'I', 'J', 'L', 'O', 'P', 'Q', 'Z')
        if any(x in pattern.upper() for x in invalid_chars):
            return 'Invalid nucleotide character found in pattern.'

    # Count tokens (residues) in pattern
    tokens = 0
    counting_mode = True

    for x in pattern:
        if x in ('(', '[', '{'):
            if counting_mode:
                tokens += 1
            counting_mode = False
        elif x in (')', ']', '}'):
            counting_mode = True
        elif counting_mode:
            tokens += 1

    # Skip min token check if pattern has ranges
    if '{' in pattern:
        return ''

    if tokens < MIN_TOKEN:
        return f"Your pattern is shorter than the minimum number of {MIN_TOKEN} residues."

    return ''


def process_pattern(
    pattern: str,
    seqtype: Optional[str],
    strand: Optional[str],
    insertion: Optional[str],
    deletion: Optional[str],
    substitution: Optional[str],
    mismatch: Optional[int]
) -> Tuple[str, str, str]:
    """
    Process and convert a pattern for searching.

    Returns:
        Tuple of (converted_pattern, complement_pattern, mismatch_option)
    """
    # Determine sequence type option
    if seqtype is None:
        seqtype = 'pep'

    if seqtype in ('pep', 'protein'):
        seq_type = 'peptide'
    elif strand and 'complement' in strand.lower():
        seq_type = 'complement'
    else:
        seq_type = 'nucleotide'

    # Convert pattern
    converter = PatternConverter(seq_type)
    converted_pattern = converter.convert(pattern)

    # Get complement pattern for DNA if searching both strands
    comp_pattern = ""
    if seqtype.lower() in ('dna', 'nuc') and (strand is None or strand.startswith('Both')):
        comp_converter = PatternConverter('complement')
        comp_pattern = comp_converter.convert(pattern)

    # Build mismatch option
    mismatch_option = ""

    if insertion and insertion.startswith('insertion'):
        mismatch_option += 'i'
    if deletion and deletion.startswith('deletion'):
        mismatch_option += 'd'
    if substitution and substitution.startswith('substitution'):
        mismatch_option += 's'

    if not mismatch_option:
        mismatch_option = 'ids'

    if mismatch is None:
        mismatch = 0

    mismatch_option = str(mismatch) + mismatch_option

    return converted_pattern, comp_pattern, mismatch_option


def find_exclusion_offset(pattern: str) -> Optional[int]:
    """
    Find the position of the first exclusion character class [^...].

    Returns:
        Offset of first exclusion or None if not found
    """
    token_re = r'\[[^\]]+\]|.(?:[*+?]|\{\d*(?:,\d*)?\})?'
    tokens = re.findall(token_re, pattern)

    # Find first negated character class
    excl_idx = None
    for i, t in enumerate(tokens):
        if t.startswith('[^'):
            excl_idx = i
            break

    if excl_idx is None:
        return None

    offset = 0
    for tok in tokens[:excl_idx]:
        if tok.startswith('['):
            offset += 1
        else:
            if len(tok) == 1:
                min_repeats = 1
            else:
                quant = tok[1:]
                if quant == '*':
                    min_repeats = 0
                elif quant == '+':
                    min_repeats = 1
                elif quant == '?':
                    min_repeats = 0
                elif quant.startswith('{'):
                    quant = quant.strip('{}').split(',')
                    try:
                        min_repeats = int(quant[0]) if quant[0] else 0
                    except (ValueError, IndexError):
                        min_repeats = 0
                else:
                    min_repeats = 1
                offset += min_repeats

    return offset


def process_output(
    record_offset_list: List[int],
    seq_name_for_offset: Dict[int, str],
    output: str,
    datafile: str,
    maxhits: Optional[int],
    beg_match: int,
    end_match: int,
    download_file: str,
    original_pattern: str
) -> Tuple[List[Dict], int, int, str]:
    """
    Process search output and generate results.

    Returns:
        Tuple of (data, unique_hits, total_hits, error_message)
    """
    # Load sequence lengths if needed
    seq_name_to_length = {}
    if end_match == 1:
        seq_name_to_length, _ = set_seq_length(datafile)

    # Track exclusion positions
    exclusion_positions = []
    for m in re.finditer(r'\[\^([^\]]+)\]', original_pattern):
        excl_pos = find_exclusion_offset(m.string[:m.start()])
        exclusion_positions.append((excl_pos, set(m.group(1))))

    # Load locus data if needed
    name_to_data = {}
    if 'orf_' in datafile:
        name_to_data = load_locus_data(DATA_DIR)

    # Load NotFeature data if needed
    seq_name_to_chr = {}
    seq_name_to_orfs = {}
    if 'Not' in datafile:
        seq_name_to_chr, seq_name_to_orfs = load_not_feature_data(datafile)

    data = []
    total_hits = 0
    unique_hits = 0
    hit_count_for_seq = {}

    # Parse maxhits
    if maxhits is None:
        maxhits = DEFAULT_MAXHITS
    elif str(maxhits).isdigit():
        maxhits = int(maxhits)
    elif str(maxhits).lower() in ('no limit', 'no+limit'):
        maxhits = MAXHITS
    else:
        maxhits = DEFAULT_MAXHITS

    # Process output lines
    for line in output.split('\n'):
        if not line.startswith('['):
            continue

        # Parse line: [start, end]: matched_text
        line = line.replace('[', '').replace(']', '')
        line = line.replace(':', '').replace(',', '')
        pieces = line.split()

        if len(pieces) < 3:
            continue

        beg = int(pieces[0])
        end = int(pieces[1])
        matching_pattern = pieces[2]

        # Check exclusion positions
        exclude_match = False
        for pos, chars in exclusion_positions:
            if pos is not None and pos < len(matching_pattern):
                if matching_pattern[pos] in chars:
                    exclude_match = True
                    break
        if exclude_match:
            continue

        # Find sequence name from offset
        offset = get_name_offset(beg, record_offset_list)
        if not str(offset).isdigit():
            continue

        seq_beg = beg - offset + 1
        seq_end = end - offset
        seq_name = seq_name_for_offset.get(offset)

        if seq_name is None:
            continue
        if beg_match == 1 and seq_beg != 1:
            continue
        if end_match == 1:
            length = seq_name_to_length.get(seq_name)
            if length is None or seq_end != length:
                continue
        if seq_name.startswith('>'):
            continue

        if seq_name.endswith(','):
            seq_name = seq_name.rstrip(',')

        # Handle NotFeature dataset
        if 'Not' in datafile:
            pieces = seq_name.split(':')
            if len(pieces) < 2:
                continue
            num = int(pieces[1].split('-')[0])
            seq_beg = seq_beg + num - 1
            seq_end = seq_end + num - 1

            if seq_name not in seq_name_to_chr or seq_name not in seq_name_to_orfs:
                continue

            row = (f"{seq_name_to_orfs.get(seq_name)}\t{seq_beg}\t{seq_end}\t"
                   f"{matching_pattern}\t{seq_name_to_chr.get(seq_name)}\t{seq_name}")
        else:
            gene, sgdid, desc = name_to_data.get(seq_name, ('', '', ''))
            row = (f"{seq_name}\t{seq_beg}\t{seq_end}\t{matching_pattern}\t"
                   f"{gene}\t{sgdid}\t{desc}")

        if seq_name not in hit_count_for_seq:
            unique_hits += 1
        if total_hits >= maxhits:
            break

        if seq_name in hit_count_for_seq:
            hit_count_for_seq[seq_name] += 1
        else:
            hit_count_for_seq[seq_name] = 1

        total_hits += 1
        data.append(row)

    # Generate output file
    file_content = []
    error_message = ''

    if 'Not' in datafile:
        header = "Chromosome\tBetweenORFtoORF\tHitNumber\tMatchPattern\tMatchStartCoord\tMatchStopCoord\n"
    elif 'orf_' in datafile:
        header = "Feature Name\tGene Name\tHitNumber\tMatchPattern\tMatchStartCoord\tMatchStopCoord\tLocusInfo\n"
    else:
        header = "Sequence Name\tHitNumber\tMatchPattern\tMatchStartCoord\tMatchStopCoord\n"

    file_content.append(header)

    new_data = []
    data.sort()

    for row in data:
        try:
            if 'Not' in datafile:
                orfs, beg, end, match_pattern, chr_name, seq_name = row.split('\t')
                count = hit_count_for_seq[seq_name]
                orfs = orfs.strip()
                new_data.append({
                    'orfs': orfs,
                    'chr': chr_name,
                    'beg': beg,
                    'end': end,
                    'count': count,
                    'seqname': seq_name,
                    'matchingPattern': match_pattern
                })
                line = f"{chr_name}\t{orfs}\t{count}\t{match_pattern}\t{beg}\t{end}\n"
            else:
                seq_name, beg, end, match_pattern, gene, sgdid, desc = row.split('\t')
                count = hit_count_for_seq.get(seq_name, 0)

                if sgdid:
                    if gene == seq_name:
                        gene = ""
                    new_data.append({
                        'seqname': seq_name,
                        'beg': beg,
                        'end': end,
                        'count': count,
                        'matchingPattern': match_pattern,
                        'gene_name': gene,
                        'sgdid': sgdid,
                        'desc': desc
                    })
                    line = f"{seq_name}\t{gene}\t{count}\t{match_pattern}\t{beg}\t{end}\t{desc}\n"
                else:
                    new_data.append({
                        'seqname': seq_name,
                        'gene_name': gene,
                        'sgdid': sgdid,
                        'beg': beg,
                        'end': end,
                        'count': count,
                        'matchingPattern': match_pattern,
                        'desc': desc
                    })
                    line = f"{seq_name}\t{count}\t{match_pattern}\t{beg}\t{end}\n"

            file_content.append(line)

        except (IndexError, ValueError) as e:
            error_message += f"Error processing row: {row}, error: {e}\n"
            continue
        except Exception as e:
            error_message += f"Unexpected error for row: {row}, error: {e}\n"
            continue

    # Write output file
    try:
        with open(download_file, 'w', encoding='utf-8') as fw:
            fw.writelines(file_content)
    except Exception as e:
        error_message += f"Error writing to file {download_file}: {e}\n"

    return new_data, unique_hits, total_hits, error_message


def run_patmatch(
    pattern: str,
    dataset: Optional[str] = None,
    seqtype: str = 'pep',
    seqname: Optional[str] = None,
    strand: Optional[str] = None,
    insertion: Optional[str] = None,
    deletion: Optional[str] = None,
    substitution: Optional[str] = None,
    mismatch: int = 0,
    max_hits: int = 500,
    request_id: str = '0'
) -> Dict[str, Any]:
    """
    Run a pattern matching search.

    Args:
        pattern: The search pattern
        dataset: Dataset name
        seqtype: Sequence type ('pep', 'protein', 'dna', 'nuc')
        seqname: Specific sequence name to retrieve
        strand: Strand option for DNA
        insertion: Allow insertions
        deletion: Allow deletions
        substitution: Allow substitutions
        mismatch: Number of allowed mismatches
        max_hits: Maximum number of hits to return
        request_id: Unique request identifier

    Returns:
        Result dictionary
    """
    tmp_file = f"patmatch.{request_id}"
    download_file = os.path.join(TMP_DIR, tmp_file)

    # Build dataset path
    if dataset:
        dataset = dataset + '.seq'
    else:
        if seqtype and seqtype in ('dna', 'nuc'):
            dataset = 'orf_dna.seq'
        else:
            dataset = 'orf_pep.seq'

    datafile = os.path.join(DATA_DIR, dataset)

    # Handle sequence retrieval
    if seqname:
        return get_sequence(datafile, seqname, DATA_DIR)

    # Clean up pattern
    pattern = cleanup_pattern(pattern)

    # Check for anchors
    beg_match = 0
    end_match = 0
    if pattern.startswith('<'):
        beg_match = 1
        pattern = pattern.replace('<', '')
    elif pattern.endswith('>'):
        end_match = 1
        pattern = pattern.replace('>', '')

    # Validate pattern
    error = check_pattern(pattern, seqtype)
    if error:
        return {"error": error}

    # Process pattern
    converted_pattern, comp_pattern, mismatch_option = process_pattern(
        pattern, seqtype, strand, insertion, deletion, substitution, mismatch
    )

    # Run search
    output = fuzzy_search(converted_pattern, datafile, mismatch_option)

    # Search complement strand if needed
    if comp_pattern:
        output2 = fuzzy_search(comp_pattern, datafile, mismatch_option)
        output = output + '\n' + output2

    # Get sequence index
    record_offset_list, seq_name_for_offset = generate_sequence_index(datafile)

    # Process output
    data, unique_hits, total_hits, error_message = process_output(
        record_offset_list, seq_name_for_offset, output,
        datafile, max_hits, beg_match, end_match, download_file, converted_pattern
    )

    # Get download URL
    download_url = ''
    if unique_hits > 0:
        try:
            download_url = get_download_url(tmp_file)
        except Exception as e:
            error_message = (error_message or '') + f" Error generating download URL: {e}"

    return {
        "hits": data,
        "uniqueHits": unique_hits,
        "totalHits": total_hits,
        "downloadUrl": download_url,
        "error_message": error_message
    }
