# -*- coding: utf-8 -*-
import traceback
import json
import os
import hashlib
from pathlib import Path
import boto3
import time
import threading
import re

from flask import send_from_directory, Response

MAX_BUFFER_SIZE = 1600000
MIN_TOKEN = 3
MINHITS = 500
MAXHITS = 100000
DEFAULT_MAXHITS = 500

binDir = '/var/www/bin/'
dataDir = '/data/patmatch/'
tmpDir = '/var/www/tmp/'
config_dir = '/var/www/conf/'
seqIndexCreateScript = binDir + 'generate_sequence_index.pl'
patternConvertScript = binDir + 'patmatch_to_nrgrep.pl'
searchScript = binDir + 'nrgrep_coords'
day = 1  ## delete temp files that are one day old

def set_download_file(filename):

    return send_from_directory(tmpDir, filename, as_attachment=True, mimetype='application/text', attachment_filename=(str(filename)))

def clean_up_temp_files():
    
    now = time.time()
    for f in os.listdir(tmpDir):
        file = os.path.join(tmpDir, f)
        if os.stat(file).st_mtime < now - day * 86400:
            if os.path.isfile(file):
                os.remove(file)


def upload_file_to_s3_async(file, filename):
    """
    Uploads file to S3 asynchronously.
    """
    try:
        upload_file_to_s3(file, filename)
    except Exception as e:
        print("Error uploading file:", e)
    finally:
        file.close()

def upload_file_to_s3(file, filename):

    filename = 'patmatch/' + filename
    
    S3_BUCKET = os.environ['S3_BUCKET']
    s3 = boto3.client('s3')
    file.seek(0)
    s3.upload_fileobj(file, S3_BUCKET, filename, ExtraArgs={'ACL': 'public-read'})
    clean_up_temp_files()
    return "https://" + S3_BUCKET + ".s3.amazonaws.com/" + filename


def get_downloadUrl(tmpFile):
    downloadFile = tmpDir + tmpFile
    thisFile = Path(str(downloadFile))

    # If the file doesn't exist, return empty URL (or raise a controlled error)
    if not thisFile.exists():
        return ""

    with thisFile.open(mode="rb") as fh:
        md5sum = hashlib.md5(fh.read()).hexdigest()

    newFileName = downloadFile
    if md5sum:
        tmpFile = md5sum + ".txt"
        newFileName = tmpDir + tmpFile
        os.rename(downloadFile, newFileName)

    f = open(newFileName, "rb")
    thread = threading.Thread(target=upload_file_to_s3_async, args=(f, tmpFile))
    thread.start()

    S3_BUCKET = os.environ['S3_BUCKET']
    return f"https://{S3_BUCKET}.s3.amazonaws.com/patmatch/{tmpFile}"


def get_downloadUrl_old(tmpFile):

    downloadFile = tmpDir + tmpFile
    thisFile = Path(str(downloadFile))
    md5sum = None
    with thisFile.open(mode="rb") as fh:
        md5sum = hashlib.md5(fh.read()).hexdigest()
    newFileName = downloadFile
    if md5sum:
        tmpFile = md5sum + ".txt"
        newFileName = tmpDir + tmpFile
        os.rename(downloadFile, newFileName)

    file = open(newFileName, "rb")
       
    # s3_url = upload_file_to_s3(file, tmpFile)

    # Start a new thread for file upload
    thread = threading.Thread(target=upload_file_to_s3_async, args=(file, tmpFile))
    thread.start()

    # Return the expected S3 URL immediately
    S3_BUCKET = os.environ['S3_BUCKET']
    return "https://" + S3_BUCKET + ".s3.amazonaws.com/" + 'patmatch/' + tmpFile


def get_config(conf):

    if conf is None:
        conf = 'patmatch'
    if not conf.endswith('.json'):
        conf = conf + '.json'
    f = open(config_dir + conf, encoding="utf-8")
    data = ''
    for line in f:
        data = data + line.strip()
    f.close()
    return json.loads(data)


def get_record_offset(datafile):

    recordOffSetList = []
    seqNm4offSet = {}
    
    cmd = seqIndexCreateScript + " < " + datafile 

    out = os.popen(cmd).read()
    
    for line in out.split('\n'):
        pieces = line.strip().split('\t')
        if len(pieces) < 2:
            continue
        offSet = int(pieces[0])
        seqNm = pieces[1]
        recordOffSetList.append(offSet)
        seqNm4offSet[offSet] = seqNm
        
    return (recordOffSetList, seqNm4offSet)


def get_name_offset(offSet, recordOffSetList):

    ### perform binary search                                                                          
    low = 0
    high = len(recordOffSetList) - 1
    while high > low:
        ### pre-condition: offSet is in recordOffSetList[$low .. $high]
        middle = int((low+high)/2)
        if recordOffSetList[middle] == offSet:
            return offSet
        elif high - low == 1:
            if offSet >= recordOffSetList[high]:
                return recordOffSetList[high]
            else:
                return recordOffSetList[low]
        elif recordOffSetList[middle] < offSet:
            low = middle
        elif recordOffSetList[middle] > offSet:
            high = middle - 1

    return recordOffSetList[low]
                     
            
def check_pattern(pattern, seqtype):

    if seqtype in ['pep', 'protein']:
        if 'u' in pattern.lower():
            return 'Invalid peptide character found in pattern.'
    else:
        if any(x in pattern.upper() for x in ('E', 'F', 'I', 'J', 'L', 'O', 'P', 'Q', 'Z')):
            return 'Invalid nucleotide character found in pattern.'

    tokens = 0
    countingMode = 1
    for x in pattern:
        if x in ['(', '[', '{']:
            if countingMode:
                tokens = tokens + 1
            countingMode = 0
        elif x in [')', ']', '}']:
            countingMode = 1
        elif countingMode:
            tokens = tokens + 1

    if '{' in pattern or '{' in pattern:
        return ''

    if tokens < MIN_TOKEN:
        return "Your pattern is shorter than the minimum number of " + str(MIN_TOKEN) + " residues."
    return ''
    
def process_pattern(pattern, seqtype, strand, insertion, deletion, substitution, mismatch):

    mismatch_option = ""

    match_characters = ['{', '}', '[', ']', '(', ')']
    for x in match_characters:
        if x in pattern:
            check_pattern(pattern, seqtype)
            break
    
    option = ''

    if seqtype is None:
        seqtype = 'pep'
    if seqtype in ['pep', 'protein']:
        option = '-p'
    elif strand and 'complement' in strand.lower():
        option = '-c'
    else:
        option = '-n'
        
    cmd = patternConvertScript + " " + option + " '" + pattern + "'"
    pattern = os.popen(cmd).read()
    
    comp_pattern = ""    
    if seqtype.lower() in ['dna', 'nuc'] and (strand is None or strand.startswith('Both')):
        cmd2 = patternConvertScript + " -c " + "'" + pattern + "'"
        comp_pattern = os.popen(cmd2).read()
    
    if insertion and insertion.startswith('insertion'): 
        mismatch_option = mismatch_option + 'i'

    if deletion and deletion.startswith('deletion'):
        mismatch_option = mismatch_option + 'd'

    if substitution and substitution.startswith('substitution'):
        mismatch_option = mismatch_option + 's'

    if mismatch_option == '':
        mismatch_option = 'ids'

    if mismatch is None:
        mismatch = 0
    
    mismatch_option = str(mismatch) + mismatch_option

    return (pattern, comp_pattern, mismatch_option)


def get_sequence(dataset, seqname):

    if '.seq' not in dataset:
        dataset = dataset + ".seq"
    if 'patmatch' not in dataset:
        dataset = dataDir + dataset
    f = open(dataset, encoding="utf-8")
    
    found = 0
    seq = ""
    defline = ""

    for line in f:
        line = line.strip()
        if line.lower().startswith('>'+seqname.lower()):
            found = 1
            defline = line
            continue
        if found == 0:
            continue
        if found == 1 and line.startswith('>'):
            break
        seq = seq + line
    
    f.close()

    defline = defline.replace('"', "'")

    return { 'defline': defline,
             'seq': seq }


def get_param(request, name, default=None):

    """
    p = request.args
    f = request.form

    return f.get(name) if f.get(name) else p.get(name)
    """

    # Check if the parameter is in the query string
    value = request.args.get(name)
    
    # If not in the query string, check the form data
    if value is None:
        value = request.form.get(name)

    # Return the value if found, otherwise return the default value
    return value if value is not None else default

    
def cleanup_pattern(pattern):
    """
    pattern = pattern.replace('%28', '(').replace('%29', ')')
    pattern = pattern.replace('%7B', '{').replace('%7D', '}')
    pattern = pattern.replace('%5B', '[').replace('%5D', ']')
    patteern = pattern.replace('%2C', ',')
    """
    # decode common URL‐escapes
    pattern = (pattern
               .replace('%28', '(').replace('%29', ')')
               .replace('%7B', '{').replace('%7D', '}')
               .replace('%5B', '[').replace('%5D', ']')
               .replace('%2C', ',')
               .replace('%5E', '^'))
    return pattern


def set_seq_length(seqNm2length, datafile):
    seqHasStop = {}
    with open(datafile, encoding="utf-8") as f:
        seq = ''
        preSeqNm = ''
        for line in f:
            if line.startswith('>'):
                if preSeqNm != '':
                    seqHasStop[preSeqNm] = seq.endswith('*')
                    # keep your original behavior: biological length (strip '*')
                    if seqHasStop[preSeqNm]:
                        seq = seq[:-1]
                    seqNm2length[preSeqNm] = len(seq)
                preSeqNm = line.replace('>', '').split(' ')[0]
                seq = ''
            else:
                seq += line.strip()
        if preSeqNm and seq:
            seqHasStop[preSeqNm] = seq.endswith('*')
            if seqHasStop[preSeqNm]:
                seq = seq[:-1]
            seqNm2length[preSeqNm] = len(seq)
    return seqHasStop


def find_exclusion_offset(pattern):
    # Improved regex to tokenize character classes, quantified atoms, and literals
    token_re = r'\[[^\]]+\]|.(?:[*+?]|\{\d*(?:,\d*)?\})?'
    tokens = re.findall(token_re, pattern)
    
    # Find the first negated character class [^...]
    excl_idx = next((i for i, t in enumerate(tokens) if t.startswith('[^')), None)
    if excl_idx is None:
        return None
    
    offset = 0
    for tok in tokens[:excl_idx]:
        if tok.startswith('['):
            # Character class always matches 1 character
            offset += 1
        else:
            # Handle quantified atoms and literals
            if len(tok) == 1:
                # Single character with no quantifier
                offset += 1
            else:
                # Extract quantifier and determine minimal repeats
                quant = tok[1:]
                if quant == '*':
                    min_repeats = 0
                elif quant == '+':
                    min_repeats = 1
                elif quant == '?':
                    min_repeats = 0
                elif quant.startswith('{'):
                    # Parse {n}, {n,m}, {n,}, {,m} cases
                    quant = quant.strip('{}').split(',')
                    try:
                        # Minimum is first number or 0 if empty
                        min_repeats = int(quant[0]) if quant[0] else 0
                    except (ValueError, IndexError):
                        # Invalid format, assume minimum 0
                        min_repeats = 0
                else:
                    # Non-quantifier suffix (shouldn't occur with proper tokenization)
                    min_repeats = 1
                offset += min_repeats
                
    return offset

    
def process_output(recordOffSetList, seqNm4offSet, output, datafile, maxhits, begMatch, endMatch, downloadFile, original_pattern):

    seqNm2length = {}
    if endMatch == 1:
        set_seq_length(seqNm2length, datafile)

    # excluded_offset = find_exclusion_offset(pattern)

    # Track exclusion positions and their characters
    exclusion_positions = []
    for m in re.finditer(r'\[\^([^\]]+)\]', original_pattern):
        exclusion_pos_pos = find_exclusion_offset(m.string[:m.start()])  # Position of this [^...]
        exclusion_positions.append( (exclusion_pos_pos, set(m.group(1))) )

    # for m in re.finditer(r'\[\^([^\]]+)\]', original_pattern):
    #    excluded_chars.update(m.group(1))

    name2data = {}
    if 'orf_' in datafile:
        with open(dataDir + "locus.txt", encoding="utf-8") as f:
            for line in f:
                pieces = line.strip().split('\t')
                seqName = pieces[0]
                geneName = pieces[1]
                sgdid = pieces[2]
                desc = ''
                if len(pieces) > 3:
                    desc = pieces[3]
                name2data[seqName] = (geneName, sgdid, desc)

    seqNm2chr = {}
    seqNm2orfs = {}
    if 'Not' in datafile:
        with open(datafile, encoding="utf-8") as f:
            for line in f:
                if line.startswith('>'):
                    # >A:2170-2479, Chr I from 2170-2479, Genome Release 64-3-1, between YAL068C and YAL067W-A
                    # /^>([^ ]+)\, Chr ([^ ]+) from .+ between ([^ ]+ and [^ ]+)/)
                    pieces = line.strip().replace('>', '').split(' ')
                    seqName = pieces[0].replace(',', '')
                    chr = pieces[2]
                    orfs = line.strip().split('between ')[1]
                    seqNm2chr[seqName] = chr
                    orfs = orfs.replace('and', '-')
                    seqNm2orfs[seqName] = orfs;
        
    data = []

    totalHits = 0
    uniqueHits = 0
    hitCount4seqNm = {}

    if maxhits is None:
        maxhits = DEFAULT_MAXHITS
    elif str(maxhits).isdigit():
        maxhits = int(maxhits)
    elif str(maxhits).lower() in ['no limit', 'no+limit']:
        maxhits = MAXHITS
    else:
        # Log unexpected value of maxhits, if needed
        maxhits = DEFAULT_MAXHITS
    
    for line in output.split('\n'):
        
        if line.startswith('['):

            line = line.replace('[', '').replace(']', '')
            line = line.replace(':', '').replace(',', '')
            pieces = line.split(' ')
            if len(pieces) < 3:
                continue
            beg = int(pieces[0])
            end = int(pieces[1])
            matchingPattern = pieces[2]

            # Check all exclusion positions
            exclude_match = False
            for pos, chars in exclusion_positions:
                if pos is not None and pos < len(matchingPattern):
                    if matchingPattern[pos] in chars:
                        exclude_match = True
                        break
            if exclude_match:
                continue
                        
            offSet = get_name_offset(beg, recordOffSetList);
            if not str(offSet).isdigit():
                continue
            seqBeg = beg - offSet + 1
            seqEnd = end - offSet
            seqNm = seqNm4offSet.get(offSet, None)
            if seqNm is None:
                continue
            if begMatch == 1 and seqBeg != 1:
                continue
            if endMatch == 1 and seqEnd != seqNm2length[seqNm]:
                continue
            if seqNm.startswith('>'):
                ## match to the fasta header line 
                continue

            if seqNm.endswith(','):
                seqNm = seqNm.rstrip(seqNm[-1])
                
            if 'Not' in datafile:
                # num = int(seqNm.split(':')[1].split('-')[0])
                pieces = seqNm.split(':')
                if len(pieces) < 2:
                    continue
                num = int(pieces[1].split('-')[0])
                seqBeg = seqBeg + num -1
                seqEnd = seqEnd + num -1
                if seqNm not in seqNm2chr or seqNm not in seqNm2orfs:
                    continue

                row = str(seqNm2orfs.get(seqNm)) + "\t" + str(seqBeg) +  "\t" + str(seqEnd) + "\t" + matchingPattern + "\t" + str(seqNm2chr.get(seqNm)) + "\t" + seqNm
                
            else:
                (gene, sgdid, desc) = name2data.get(seqNm, ('', '', ''))
                row = seqNm + "\t" + str(seqBeg) + "\t" + str(seqEnd) + "\t" + matchingPattern + "\t" + gene + "\t" + sgdid + "\t" + desc

            if seqNm not in hitCount4seqNm:
                uniqueHits = uniqueHits + 1
            if totalHits >= maxhits:
                break

            if seqNm in hitCount4seqNm:
                hitCount4seqNm[seqNm] = hitCount4seqNm[seqNm] + 1
            else:
                hitCount4seqNm[seqNm] = 1
            totalHits = totalHits + 1

            data.append(row)

    file_content = []
    header_line = ""

    if 'Not' in datafile:
        header_line = "Chromosome\tBetweenORFtoORF\tHitNumber\tMatchPattern\tMatchStartCoord\tMatchStopCoord\n"    
    elif 'orf_' in datafile:
        header_line = "Feature Name\tGene Name\tHitNumber\tMatchPattern\tMatchStartCoord\tMatchStopCoord\tLocusInfo\n"
    else:
        header_line = "Sequence Name\tHitNumber\tMatchPattern\tMatchStartCoord\tMatchStopCoord\n"

    file_content.append(header_line)

    newData = []

    data.sort()

    error_message = ''
    
    for row in data:
        try:
            if 'Not' in datafile:
                [orfs, beg, end, matchPattern, chr, seqNm] = row.split('\t')
                count = hitCount4seqNm[seqNm]
                orfs = orfs.strip()
                newData.append({ 'orfs': orfs,
                             'chr': chr,
                             'beg': beg,
                             'end': end,
                             'count': count,
                             'seqname': seqNm,
                             'matchingPattern': matchPattern })
                line = chr + "\t" + orfs + "\t" + str(count) + "\t" + matchPattern + "\t" + beg + "\t" + end + "\n"
            else:

                [seqNm, beg, end, matchPattern, gene, sgdid, desc] = row.split('\t')
        
                count = hitCount4seqNm.get(seqNm, 0)
        
                if sgdid != "":
                    if gene == seqNm:
                        gene = ""
                    newData.append({ 'seqname': seqNm,
                                 'beg': beg,
                                 'end': end,
                                 'count': count,
                                 'matchingPattern': matchPattern,
                                 'gene_name': gene,
                                 'sgdid': sgdid,
                                 'desc': desc })
                    line = seqNm + "\t" + gene + "\t" + str(count) + "\t" + matchPattern + "\t" + beg + "\t" + end + "\t" + desc + "\n"
                else:
                    newData.append({ 'seqname': seqNm,
                                 'gene_name': gene,
                                 'sgdid': sgdid,
                                 'beg': beg,
                                 'end': end,
                                 'count': count,
                                 'matchingPattern': matchPattern,
                                 'desc': desc })
                    line = seqNm + "\t" + str(count) + "\t" + matchPattern + "\t" + beg + "\t" + end + "\n"
                file_content.append(line)
        except MemoryError as e:
            error_message += "Memory Error: " + str(e) + "\n"
            break
        except OSError as e:
            error_message += "OS Error: " + str(e) + "\n"
            continue
        except (IndexError, ValueError) as e:
            # Handle specific errors like IndexError or ValueError
            error_message += "Error processing row: " + str(row) + "error: " + str(e) + "\n"
            continue
        except Exception as e:
            # Catch any other unforeseen errors
            error_message += "Unexpected error for row: " + str(row) + "error: " + str(e) + "\n"
            error_message += "Traceback: " + str(traceback.format_exc()) + "\n"
            continue
    try:
        with open(downloadFile, "w", encoding='utf-8') as fw:
            fw.writelines(file_content)
    except MemoryError as e:
        error_message += "Memory Error during file writing: " + str(e) + "\n"
    except OSError as e:
        error_message += "OS Error during file writing: " + str(e) + "\n"
    except UnicodeEncodeError as e:
        error_message += "Unicode Encoding Error: " + str(e) + "\n"
    except Exception as e:
        error_message += "Error writing to file " + downloadFile + ":" + str(e)
        error_message += "Traceback: " + str(traceback.format_exc()) + "\n"

    return (newData, uniqueHits, totalHits, error_message)


def run_patmatch(request, id):

    tmpFile = "patmatch." + id
    downloadFile = tmpDir + tmpFile

    p = request.args
    f = request.form
    
    dataset = get_param(request, 'dataset')
    seqtype = get_param(request, 'seqtype')
    seqname = get_param(request, 'seqname')

    if seqtype is None:
        seqtype = 'pep'
    
    if dataset:
        dataset = dataset + ".seq"
    else:
        if seqtype and seqtype in ['dna', 'nuc']:
            dataset = "orf_dna.seq"
        else:
            dataset = "orf_pep.seq"
	
    datafile = dataDir + dataset

    if seqname:
        data = get_sequence(datafile, seqname)
        return data

    pattern = cleanup_pattern(get_param(request, 'pattern'))

    begMatch = 0
    endMatch = 0
    if pattern.startswith('<'):
        begMatch = 1
        pattern = pattern.replace('<', '')
    elif pattern.endswith('>'):
        endMatch = 1
        pattern = pattern.replace('>', '')
    
    error = check_pattern(pattern, seqtype)
    if error:
        return { "error": error }
    
    (pattern, comp_pattern, option) = process_pattern(pattern,
                                                      get_param(request, 'seqtype'),
                                                      get_param(request, 'strand'),
                                                      get_param(request, 'insertion'),
                                                      get_param(request, 'deletion'),
                                                      get_param(request, 'substitution'),
                                                      get_param(request, 'mismatch'))

    maxBufferSize = MAX_BUFFER_SIZE
    
    nrgrep = searchScript + " -i -b " + str(maxBufferSize) + " -k " + option + " '" + pattern + "' '" + datafile + "'"
    output = os.popen(nrgrep).read()

    nrgrep2 = ''
    output2 = ''
    if comp_pattern:
        nrgrep2 = searchScript + " -i -b " +str( maxBufferSize) + " -k " + option + " '" + comp_pattern +"' '" + datafile + "'"
        output2 = os.popen(nrgrep2).read()
        output = output + "\n" + output2

    # return { "nrgrep": nrgrep,
    #         "nrgrep2": nrgrep2,
    #         "output": output }
        
    (recordOffSetList, seqNm4offSet) = get_record_offset(datafile)

    # return { "nrgrep": nrgrep,
    #         "nrgrep2": nrgrep2,
    #         "recordOffSetlist": recordOffSetList,
    #         "seqNm4offSet": seqNm4offSet }
    
    (data, uniqueHits, totalHits, error_message) = process_output(recordOffSetList, seqNm4offSet, output,
                                                                  datafile, get_param(request, 'max_hits'),
                                                                  begMatch, endMatch, downloadFile, pattern)

    downloadUrl = ''
    if uniqueHits > 0:
        downloadUrl = get_downloadUrl(tmpFile)
        
    return { "hits": data,
             "uniqueHits": uniqueHits,
             "totalHits": totalHits,
             "downloadUrl": downloadUrl,
             "error_message": error_message }
    
    

















