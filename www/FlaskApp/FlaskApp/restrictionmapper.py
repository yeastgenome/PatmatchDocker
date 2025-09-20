import json
import os
import re
from patmatch import get_downloadUrl

binDir = '/var/www/bin/'
dataDir = '/data/restriction_mapper/'
tmpDir = "/var/www/tmp/"

scan4matches = binDir + "scan_for_matches"
fastafile = dataDir + "orf_genomic.seq"

def get_downloadURLs(cutSiteFile, notCutFile):

    return (get_downloadUrl(cutSiteFile), get_downloadUrl(notCutFile))
    
def get_sequence(name):

    name = name.replace('SGD:', '')

    f = open(fastafile, encoding="utf-8")
    
    seq = ""
    defline = ""
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            pieces = line.split(' ')
            if pieces[0].replace('>', '').lower() == name.lower() or pieces[1].lower() == name.lower() or pieces[2].replace('SGDID:', '').replace(',', '').lower() == name.lower():
                defline = line
            continue
        elif defline != '':
            seq = line
        if seq != '':
            break

    f.close()

    defline = defline.replace('"', "'")
    
    return (defline, seq)

def write_seqfile(defline, seq, seqfile):
    
    fw = open(seqfile, "w")

    ## remove all non-alphabet chars from seq string
    regex = re.compile('[^a-zA-Z]')
    seq = regex.sub('', seq)
        
    fw.write(defline + "\n")
    fw.write(seq + "\n")
    fw.close()

    seqNm = "Unnamed"
    chrCoords = ""
    # >YAL067C SEO1 SGDID:S000000062, Chr I from 9016-7235, Genome Release ...
    if "SGDID:" in defline and 'Genome Release' in defline:
        pieces = defline.replace('>', '').split(' ')
        systematic_name = pieces[0]
        gene_name = pieces[1]
        chrCoords = defline.split(', ')[1]
        seqNm = systematic_name
        if gene_name:
            seqNm = gene_name + "/" + systematic_name

    return (seqNm, chrCoords, len(seq))

def set_enzyme_file(enzymetype):

    if enzymetype is None:
        return dataDir + 'rest_enzymes'

    if "Six-base" in enzymetype:
        return dataDir + 'rest_enzymes.6base'

    if "blunt" in enzymetype:
        return dataDir + 'rest_enzymes.blunt'
    
    if "3" in enzymetype:
        return dataDir + 'rest_enzymes.3'

    if "5" in enzymetype:
        return dataDir + 'rest_enzymes.5'
    
    return dataDir + 'rest_enzymes'

def do_search(enzymefile, patfile, outfile, seqfile):

    f = open(enzymefile, encoding="utf-8")

    ## reset file
    fw = open(outfile, "w")
    fw.close()
    
    error_msg = ""
    for line in f:

        pieces = line.strip().split(' ')
        enzyme = pieces[0]
        offset = pieces[1]
        pat = pieces[2]
        overhang = pieces[3]
        fw = open(patfile, "w")
        fw.write(pat + "\n")
        fw.close()
        fw = open(outfile, "a")
        fw.write(">>" + enzyme + ": " + str(offset) + " " + overhang + " " + pat + "\n")
        fw.close()

        cmd = scan4matches + " -c " + patfile + " < " + seqfile + " >> " + outfile
        err = os.system(cmd)
        
        if err < 0: 
            error_msg = "RestrmictionMapper: problem running " + scan4matches + " returned " + str(err)
            break
        
    f.close()

    if error_msg:
        return error_msg
    
    if os.path.isfile(outfile):
        return ""
    else:
        return "No " + outfile + " generated in do_search!"

def set_enzyme_types(enzymeHash, enzymeType):

    enzymeFile = "rest_enzymes.blunt"
    if '3' in enzymeType:
        enzymeFile = "rest_enzymes.3"
    elif '5' in enzymeType:
        enzymeFile = "rest_enzymes.5"

    f = open(dataDir + enzymeFile, encoding="utf-8")
    for line in f:
        pieces = line.strip().split(' ')
        enzymeHash[pieces[0]] = enzymeType
    f.close()
    
def process_data(seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut):
    # --- DEBUG controls (no signature change) ---
    DEBUG = os.environ.get('RESTMAP_DEBUG', '0') == '1'
    CP = (os.environ.get('RESTMAP_CP') or '').upper()  # PARSE | WROTE_NO_CUT | EARLY_RETURN

    # Normalize once to be robust
    etype = (enzymetype or '').strip()
    etype_l = etype.lower()

    dataHash = {}
    offset = {}
    overhang = {}
    recognition_seq = {}
    notCutEnzyme = []
    # ---- PARSE OUTFILE ----
    with open(outfile, encoding="utf-8") as f:
        preLine = ''
        enzyme = ''
        for line in f:
            s = line.strip()
            if s.startswith('>>'):
                pieces = s.split(' ')
                enzyme = pieces[0].replace('>>', '').replace(':', '')
                # defensively guard pieces length
                if len(pieces) >= 4:
                    offset[enzyme] = pieces[1]
                    overhang[enzyme] = pieces[2]
                    recognition_seq[enzyme] = pieces[3]
                else:
                    # malformed header line
                    offset[enzyme] = offset.get(enzyme, '0')
                    overhang[enzyme] = overhang.get(enzyme, '0')
                    recognition_seq[enzyme] = recognition_seq.get(enzyme, 'N/A')

                # For "all" or "do not" modes, track enzymes that had no '>' lines
                if etype_l in ('all', '') or etype_l.startswith('enzymes that do not'):
                    if preLine.startswith('>>'):
                        prev = preLine.replace('>>', '').replace(':', '').split(' ')[0]
                        if prev not in notCutEnzyme:
                            notCutEnzyme.append(prev)

            elif s.startswith('>'):
                # expected format: "...:[beg,end]"
                parts = s.split(':', 1)
                if len(parts) == 2:
                    coords = parts[1].replace('[', '').replace(']', '')
                    if enzyme in dataHash and dataHash[enzyme]:
                        dataHash[enzyme] = dataHash[enzyme] + ':' + coords
                    else:
                        dataHash[enzyme] = coords
            preLine = s
    # handle last header possibly having no cuts
    if etype_l in ('all', '') or etype_l.startswith('enzymes that do not'):
        if preLine.startswith('>>'):
            last = preLine.replace('>>', '').replace(':', '').split(' ')[0]
            if last not in notCutEnzyme:
                notCutEnzyme.append(last)

    # ---- CHECKPOINT: after parse ----
    if DEBUG and CP == 'PARSE':
        # Return early so you can inspect what we parsed
        dbg = {
            "_debug": {
                "checkpoint": "PARSE",
                "etype": etype,
                "etype_lower": etype_l,
                "num_headers": len(offset),
                "num_with_cuts": sum(1 for k,v in dataHash.items() if v),
                "num_no_cut": len(notCutEnzyme),
                "sample_headers": list(offset.keys())[:10],
                "sample_no_cut": notCutEnzyme[:10],
                "has_data_for_first_header": bool(dataHash.get(next(iter(offset), ''), '')),
            }
        }
        # Keep the return shape the same: (data, notCutEnzyme)
        return (dbg, notCutEnzyme)

    # ---- WRITE "DO NOT CUT" LIST ----
    with open(downloadfile4notCut, 'w') as fw:
        for e in sorted(notCutEnzyme):
            fw.write(e + "\n")

    # ---- CHECKPOINT: after writing no-cut file ----
    if DEBUG and CP == 'WROTE_NO_CUT':
        dbg = {
            "_debug": {
                "checkpoint": "WROTE_NO_CUT",
                "etype": etype,
                "etype_lower": etype_l,
                "num_no_cut": len(notCutEnzyme),
                "first_5_no_cut": sorted(notCutEnzyme)[:5],
                "downloadfile4notCut_exists": os.path.isfile(downloadfile4notCut),
                "downloadfile4notCut": downloadfile4notCut,
            }
        }
        return (dbg, notCutEnzyme)

    # ---- EARLY RETURN for "do not" modes ----
    if etype_l.startswith('enzymes that do not'):
        if DEBUG and CP == 'EARLY_RETURN':
            dbg = {
                "_debug": {
                    "checkpoint": "EARLY_RETURN",
                    "etype": etype,
                    "etype_lower": etype_l,
                    "num_no_cut": len(notCutEnzyme),
                }
            }
            return (dbg, notCutEnzyme)
        # Normal early return
        return ({}, notCutEnzyme)
    # ---- FILTERS for "cut once/twice" ----
    if 'cut' in etype_l:
        cutLimit = 2 if 'twice' in etype_l else 1
        newDataHash = {}
        for key, coords_str in dataHash.items():
            if not coords_str:
                continue
            wCut = cCut = 0
            for coordPair in coords_str.split(':'):
                try:
                    beg, end = map(int, coordPair.split(','))
                except ValueError:
                    continue
                if beg < end:
                    wCut += 1
                else:
                    cCut += 1
            if (cCut == cutLimit and wCut <= cutLimit) or (wCut == cutLimit and cCut <= cutLimit):
                newDataHash[key] = coords_str
        dataHash = newDataHash

    # ---- ENZYME TYPE MAP ----
    enzyme_type = {}
    set_enzyme_types(enzyme_type, "3' overhang")
    set_enzyme_types(enzyme_type, "5' overhang")
    set_enzyme_types(enzyme_type, "blunt end")

    data = {}
    with open(downloadfile4cutSite, 'w') as fw:
        fw.write("Enzyme\toffset (bp)\toverhang (bp)\trecognition sequence\tenzyme type\tnumber of cuts\tordered fragment size\tsorted fragment size\tcut site on watson strand\tcut site on crick strand\n")
        for enzyme in sorted(dataHash):
            if not dataHash[enzyme]:
                continue

            et_label = enzyme_type.get(enzyme, 'unknown')
            # subset filters for overhang/blunt
            if (('overhang' in etype_l) or ('blunt' in etype_l)) and et_label.lower() != etype_l:
                continue
            cutW, cutC, cutAll = [], [], []
            for position in dataHash[enzyme].split(':'):
                try:
                    beg, end = map(int, position.split(','))
                except ValueError:
                    continue
                if beg < end:  # watson
                    cutSite = beg + int(offset[enzyme]) - 1
                    if cutSite not in cutW:
                        cutW.append(cutSite)
                else:  # crick
                    beg, end = end, beg
                    cutSite = beg + int(offset[enzyme]) + int(overhang[enzyme]) - 1
                    if cutSite not in cutC:
                        cutC.append(cutSite)
                if cutSite not in cutAll:
                    cutAll.append(cutSite)
            cutAll.append(seqLen)

            preCutSite = 0
            seen = set()
            cutFragments = []
            for cs in sorted(cutAll):
                cutSize = cs - preCutSite
                if cutSize and cutSize not in seen:
                    cutFragments.append(cutSize)
                    seen.add(cutSize)
                preCutSite = cs

            cutSiteW = ", ".join(str(x) for x in sorted(cutW))
            cutSiteC = ", ".join(str(x) for x in sorted(cutC))
            fragmentsReal = ", ".join(str(x) for x in cutFragments)
            fragments = ", ".join(str(x) for x in sorted(cutFragments, reverse=True))
            cutNum = len(cutFragments) - 1

            fw.write(f"{enzyme}\t{offset.get(enzyme,'?')}\t{overhang.get(enzyme,'?')}\t{recognition_seq.get(enzyme,'?')}\t{et_label}\t{cutNum}\t{fragmentsReal}\t{fragments}\t{cutSiteW}\t{cutSiteC}\n")

            data[enzyme] = {
                "cut_site_on_watson_strand": cutSiteW,
                "cut_site_on_crick_strand": cutSiteC,
                "fragment_size": fragments,
                "fragment_size_in_real_order": fragmentsReal,
                "offset": offset.get(enzyme, '?'),
                "overhang": overhang.get(enzyme, '?'),
                "recognition_seq": recognition_seq.get(enzyme, '?'),
                "enzyme_type": et_label,
            }

    return (data, notCutEnzyme)


def run_restriction_site_search(request, id):

    patfile = tmpDir + "patfile." + id + ".txt"
    outfile = tmpDir + "outfile." + id + ".txt"
    seqfile = tmpDir + "seqfile." + id + ".txt"

    cutSiteFile = "restrictionmapper." + id
    notCutFile = "restrictionmapper_not_cut_enzyme." + id

    downloadfile4cutSite = tmpDir + cutSiteFile
    downloadfile4notCut = tmpDir + notCutFile

    p = request.args
    f = request.form

    seq = f.get('seq') if f.get('seq') else p.get('seq')
    name = f.get('name') if f.get('name') else p.get('name')
    enzymetype = f.get('type') if f.get('type') else p.get('type', 'ALL')
    enzymetype = enzymetype.replace('+', ' ').replace("%27", "'")

    if enzymetype.startswith('3'):
        enzymetype = "3' overhang"
    elif enzymetype.startswith('5'):
        enzymetype = "5' overhang"

    defline = None
    if seq:
        defline = ">Unnamed sequence"
    else:
        (defline, seq) = get_sequence(name)

    (seqNm, chrCoords, seqLen) = write_seqfile(defline, seq, seqfile)
    
    enzymefile = set_enzyme_file(enzymetype)
    
    err = do_search(enzymefile, patfile, outfile, seqfile)
        
    if err == '':
        ## key is the enzyme
        (data, notCutEnzymeList) = process_data(seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut)
        (downloadUrl4cutSite, downloadUrl4notCut) = get_downloadURLs(cutSiteFile, notCutFile)
        
        return { "data": data,
                 "seqName": seqNm,
                 "chrCoords": chrCoords,
                 "seqLength": seqLen,
                 "notCutEnzyme": notCutEnzymeList,
                 "downloadUrl": downloadUrl4cutSite,
                 "downloadUrl4notCutEnzyme": downloadUrl4notCut }
    else:
        return { "ERROR": err,
                 "seqName": seqNm,
                 "chrCoords": chrCoords,
                 "seqLength": seqLen,
                 "notCutEnzyme": [],
                 "downloadUrl": '',
                 "downloadUrl4notCutEnzyme": '' }
    
