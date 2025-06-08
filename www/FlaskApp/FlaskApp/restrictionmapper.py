import json
import os
import re
from patmatch import get_downloadUrl

# Directories and executables
enzymedir = '/data/restriction_mapper/'
binDir = '/var/www/bin/'
dataDir = enzymedir
tmpDir = "/var/www/tmp/"
scan4matches = os.path.join(binDir, "scan_for_matches")
fastafile = os.path.join(dataDir, "orf_genomic.seq")


def normalize_mode(enzymetype: str) -> str:
    """
    Derive a canonical mode string for filtering:
      - "all"
      - "six-base"
      - "3' overhang"
      - "5' overhang"
      - "blunt end"
    """
    et = enzymetype.lower().replace('+', ' ').replace('%27', "'").strip()
    if et in ('', 'all'):
        return 'all'
    if 'six-base' in et:
        return 'six-base'
    if et.startswith("3'") or et.startswith('3') or '3 overhang' in et:
        return "3' overhang"
    if et.startswith("5'") or et.startswith('5') or '5 overhang' in et:
        return "5' overhang"
    if 'blunt' in et:
        return 'blunt end'
    # fallback to all\ n    return 'all'


def get_downloadURLs(cutSiteFile, notCutFile):
    return (get_downloadUrl(cutSiteFile), get_downloadUrl(notCutFile))


def get_sequence(name):
    name = name.replace('SGD:', '')
    seq = ''
    defline = ''
    with open(fastafile, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                pieces = line.split(' ')
                if (pieces[0].lstrip('>').lower() == name.lower() or
                    any(p.lower() == name.lower() for p in pieces[1:3])):
                    defline = line
                continue
            if defline:
                seq = line
                break
    defline = defline.replace('"', "'") or '>Unnamed sequence'
    return defline, seq


def write_seqfile(defline, seq, seqfile):
    # sanitize sequence: keep only letters
    seq_clean = re.sub('[^a-zA-Z]', '', seq)
    with open(seqfile, 'w') as fw:
        fw.write(defline + "\n")
        fw.write(seq_clean + "\n")

    seqNm = 'Unnamed'
    chrCoords = ''
    if 'SGDID:' in defline and 'Genome Release' in defline:
        parts = defline.lstrip('>').split(' ')
        systematic = parts[0]
        gene = parts[1]
        chrCoords = defline.split(', ')[1]
        seqNm = f"{gene}/{systematic}" if gene else systematic

    return seqNm, chrCoords, len(seq_clean)


def set_enzyme_file(enzymetype):
    mode = normalize_mode(enzymetype)
    if mode == 'six-base':
        filename = 'rest_enzymes.6base'
    elif mode == "3' overhang":
        filename = 'rest_enzymes.3'
    elif mode == "5' overhang":
        filename = 'rest_enzymes.5'
    elif mode == 'blunt end':
        filename = 'rest_enzymes.blunt'
    else:
        filename = 'rest_enzymes.all'
    return os.path.join(dataDir, filename)


def do_search(enzymefile, patfile, outfile, seqfile):
    if not os.path.isfile(enzymefile):
        return f"Error: enzyme file not found: {enzymefile}"

    # clear output
    open(outfile, 'w').close()
    error_msg = ''

    with open(enzymefile, encoding='utf-8') as f:
        for line in f:
            enzyme, offset, pat, overhang = line.strip().split(' ')
            # write pattern
            with open(patfile, 'w') as fw:
                fw.write(pat + '\n')
            with open(outfile, 'a') as fw:
                fw.write(f">>{enzyme}: {offset} {overhang} {pat}\n")
            cmd = f"{scan4matches} -c {patfile} < {seqfile} >> {outfile}"
            if os.system(cmd) < 0:
                error_msg = f"RestrictionMapper: problem running {scan4matches} returned {err}"
                break

    if error_msg:
        return error_msg
    return '' if os.path.isfile(outfile) else f"No {outfile} generated in do_search!"


def set_enzyme_types(enzymeHash, enzymeType):
    data_files = {
        "3' overhang": 'rest_enzymes.3',
        "5' overhang": 'rest_enzymes.5',
        'blunt end': 'rest_enzymes.blunt',
    }
    fname = data_files.get(enzymeType)
    if not fname:
        return
    path = os.path.join(dataDir, fname)
    if not os.path.isfile(path):
        return
    with open(path, encoding='utf-8') as f:
        for line in f:
            key = line.strip().split(' ')[0]
            enzymeHash[key] = enzymeType


def process_data(seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut):
    # parse results
    dataHash = {}
    notCutEnzyme = []
    offset = {}
    overhang = {}
    recognition_seq = {}

    with open(outfile, encoding='utf-8') as f:
        preLine = ''
        enzyme = ''
        for line in f:
            if line.startswith('>>'):
                parts = line.strip().split(' ')
                enzyme = parts[0].lstrip('>>').rstrip(':')
                offset[enzyme] = parts[1]
                overhang[enzyme] = parts[2]
                recognition_seq[enzyme] = parts[3]
                if enzymetype.lower().startswith('enzymes that do not'):
                    if preLine.startswith('>>'):
                        prev = preLine.lstrip('>>').rstrip(':').split(' ')[0]
                        if prev not in notCutEnzyme:
                            notCutEnzyme.append(prev)
            elif line.startswith('>'):
                coords = line.split('[')[-1].rstrip(']\n').replace(',', ',')
                dataHash.setdefault(enzyme, coords)
                if enzyme in dataHash:
                    dataHash[enzyme] += ':' + coords
            preLine = line.strip()

    # finalize not-cut list
    if notCutEnzyme and not fake: pass  # removed for brevity

    mode = normalize_mode(enzymetype)

    # if non-cut mode
    if enzymetype.lower().startswith('enzymes that do not'):
        with open(downloadfile4notCut, 'w') as fw:
            for e in sorted(notCutEnzyme): fw.write(e + '\n')
        return {}, notCutEnzyme

    # apply cut count filter for "cut" modes
    if 'cut' in enzymetype.lower():
        limit = 2 if 'twice' in enzymetype.lower() else 1
        filtered = {}
        for key, coords_str in dataHash.items():
            wCuts, cCuts = 0, 0
            for pair in coords_str.split(':'):
                b, e = map(int, pair.split(','))
                if b < e: wCuts += 1
                else: cCuts += 1
            if (wCuts == limit and cCuts <= limit) or (cCuts == limit and wCuts <= limit):
                filtered[key] = coords_str
        dataHash = filtered

    # build enzyme_type map
    enzyme_type = {}
    for etype in ("3' overhang", "5' overhang", 'blunt end'):
        set_enzyme_types(enzyme_type, etype)

    # write cut-site output header
    with open(downloadfile4cutSite, 'w') as fw:
        fw.write(
            "Enzyme\toffset (bp)\toverhang (bp)\trecognition sequence\t"
            "enzyme type\tnumber of cuts\tordered fragment size\t"
            "sorted fragment size\tcut site on watson strand\tcut site on crick strand\n"
        )

    result = {}
    cutList = []
    with open(downloadfile4cutSite, 'a') as fw:
        for enzyme in sorted(dataHash):
            # filter by mode if needed
            if mode not in ('all', 'six-base') and enzyme_type.get(enzyme) != mode:
                continue

            coords_all = []
            cutW = []
            cutC = []
            for pair in dataHash[enzyme].split(':'):
                b, e = map(int, pair.split(','))
                if b < e:
                    site = b + int(offset[enzyme]) - 1
                    cutW.append(site)
                else:
                    site = min(b, e) + int(offset[enzyme]) + int(overhang[enzyme]) - 1
                    cutC.append(site)
                coords_all.append(site)

            coords_all.append(seqLen)
            coords_all.sort()

            # fragment sizes
            fragments = []
            prev = 0
            for c in coords_all:
                size = c - prev
                if size > 0:
                    fragments.append(size)
                prev = c

            ordered = ', '.join(map(str, fragments))
            sorted_frags = ', '.join(map(str, sorted(fragments, reverse=True)))
            fw.write(
                f"{enzyme}\t{offset[enzyme]}\t{overhang[enzyme]}\t"
                f"{recognition_seq[enzyme]}\t{enzyme_type[enzyme]}\t"
                f"{len(fragments)-1}\t{ordered}\t{sorted_frags}\t"
                f"{', '.join(map(str, cutW))}\t{', '.join(map(str, cutC))}\n"
            )

            result[enzyme] = {
                'offset': offset[enzyme],
                'overhang': overhang[enzyme],
                'recognition_seq': recognition_seq[enzyme],
                'enzyme_type': enzyme_type[enzyme],
                'cut_site_on_watson_strand': ', '.join(map(str, cutW)),
                'cut_site_on_crick_strand': ', '.join(map(str, cutC)),
                'fragment_size_in_real_order': ordered,
                'fragment_size': sorted_frags,
            }

    return result, notCutEnzyme


def run_restriction_site_search(request, id):
    patfile = os.path.join(tmpDir, f"patfile.{id}.txt")
    outfile = os.path.join(tmpDir, f"outfile.{id}.txt")
    seqfile = os.path.join(tmpDir, f"seqfile.{id}.txt")

    cutSiteName = f"restrictionmapper.{id}"
    notCutName = f"restrictionmapper_not_cut_enzyme.{id}"
    downloadfile4cutSite = os.path.join(tmpDir, cutSiteName)
    downloadfile4notCut = os.path.join(tmpDir, notCutName)

    p = request.args
    f = request.form
    seq = f.get('seq') or p.get('seq')
    name = f.get('name') or p.get('name')
    enzymetype = (f.get('type') or p.get('type', 'all'))

    mode = normalize_mode(enzymetype)

    if seq:
        defline = ">Unnamed sequence"
    else:
        defline, seq = get_sequence(name)

    seqNm, chrCoords, seqLen = write_seqfile(defline, seq, seqfile)
    enzymefile = set_enzyme_file(enzymetype)

    err = do_search(enzymefile, patfile, outfile, seqfile)
    if err:
        return {
            "ERROR": err,
            "seqName": seqNm,
            "chrCoords": chrCoords,
            "seqLength": seqLen,
            "notCutEnzyme": [],
            "downloadUrl": '',
            "downloadUrl4notCutEnzyme": ''
        }

    data, notCutList = process_data(seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut)
    url_cut, url_not = get_downloadURLs(cutSiteName, notCutName)

    return {
        "data": data,
        "seqName": seqNm,
        "chrCoords": chrCoords,
        "seqLength": seqLen,
        "notCutEnzyme": notCutList,
        "downloadUrl": url_cut,
        "downloadUrl4notCutEnzyme": url_not
    }
