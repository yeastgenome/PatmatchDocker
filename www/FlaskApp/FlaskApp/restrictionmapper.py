import json
import os
import re
from patmatch import get_downloadUrl

binDir = '/var/www/bin/'
dataDir = '/data/restriction_mapper/'
tmpDir = "/var/www/tmp/"

scan4matches = binDir + "scan_for_matches"
fastafile = dataDir + "orf_genomic.seq"


def _normalize_type(t):
    """Map UI variants to server-friendly values."""
    if not t:
        return 'all'
    s = str(t).strip().lower().replace('+', ' ')
    s = s.replace('%27', "'")  # 3%27 -> 3'
    if s in {"3'overhang", "3' overhang", "3 overhang"}:
        return "3' overhang"
    if s in {"5'overhang", "5' overhang", "5 overhang"}:
        return "5' overhang"
    if s in {"six-base cutters", "six base cutters", "6-base cutters"}:
        return "Six-base cutters"
    if s in {"enzymes that do not cut", "do not cut", "no cut", "not cut"}:
        return "enzymes that do not cut"
    if s in {'blunt', 'blunt end'}:
        return 'blunt end'
    if s in {'cut once', 'cut 1', 'once'}:
        return 'cut once'
    if s in {'cut twice', 'cut 2', 'twice'}:
        return 'cut twice'
    if s == 'all':
        return 'all'
    return t  # fallback – your existing checks still handle this

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

def process_data(seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut, enzymefile):
    # 1) Parse which enzymes CUT (from outfile)
    dataHash = {}         # enzyme -> "beg,end:beg,end:..."
    offset = {}
    overhang = {}
    recognition_seq = {}

    with open(outfile, encoding="utf-8") as f:
        enzyme = ''
        for raw in f:
            line = raw.strip()
            if line.startswith('>>'):
                # ">>EcoRI: 1 4 GAATTC"
                parts = line.replace('>>', '').replace(':', '').split()
                # be defensive: some patterns could include spaces; keep indexes safe
                if len(parts) >= 4:
                    enzyme = parts[0]
                    offset[enzyme] = parts[1]
                    overhang[enzyme] = parts[2]
                    # everything after index 2 back together (handles rare spaces)
                    recognition_seq[enzyme] = ' '.join(parts[3:])
                else:
                    # malformed line; skip
                    enzyme = ''
            elif line.startswith('>') and enzyme:
                # ">{header}:{beg,end}"  (you were doing split(':')[1])
                try:
                    coords = line.split(':', 1)[1].replace('[', '').replace(']', '')
                    if enzyme in dataHash:
                        dataHash[enzyme] = dataHash[enzyme] + ':' + coords
                    else:
                        dataHash[enzyme] = coords
                except Exception:
                    pass

    # 2) Build the universe of enzymes for THIS selection (from the enzyme file you used)
    all_enzymes = []
    with open(enzymefile, encoding="utf-8") as ef:
        for raw in ef:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            # file format: "Enzyme offset pattern overhang"
            name = line.split()[0]
            all_enzymes.append(name)

    cut_enzymes = set(dataHash.keys())
    notCutEnzyme = sorted([e for e in all_enzymes if e not in cut_enzymes])

    # 3) Always write the "do not cut" download file
    with open(downloadfile4notCut, 'w') as fw:
        for e in notCutEnzyme:
            fw.write(e + "\n")

    # 4) If the user asked specifically for "enzymes that do not cut", return just that
    if enzymetype.lower() == 'enzymes that do not cut':
        return ({}, notCutEnzyme)

    # --- existing filtering logic below (unchanged) ---
    if "cut" in enzymetype.lower():
        cutLimit = 1 if 'once' in enzymetype.lower() else 2 if 'twice' in enzymetype.lower() else 1
        newDataHash = {}
        for key, coord_blob in dataHash.items():
            coords = coord_blob.split(':')
            wCut = cCut = 0
            for coordPair in coords:
                beg, end = map(int, coordPair.split(','))
                if beg < end:
                    wCut += 1
                else:
                    cCut += 1
            if (cCut == cutLimit and wCut <= cutLimit) or (wCut == cutLimit and cCut <= cutLimit):
                newDataHash[key] = coord_blob
        dataHash = newDataHash

    # Annotate enzyme types (unchanged)
    enzyme_type = {}
    set_enzyme_types(enzyme_type, "3' overhang")
    set_enzyme_types(enzyme_type, "5' overhang")
    set_enzyme_types(enzyme_type, "blunt end")

    data = {}

    with open(downloadfile4cutSite, 'w') as fw:
        fw.write("Enzyme\toffset (bp)\toverhang (bp)\trecognition sequence\tenzyme type\tnumber of cuts\tordered fragment size\tsorted fragment size\tcut site on watson strand\tcut site on crick strand\n")

        for enzyme in sorted(dataHash):
            # filter by overhang/blunt if requested
            if ("overhang" in enzymetype.lower() or "blunt" in enzymetype.lower()):
                if enzyme_type.get(enzyme) != enzymetype:
                    continue

            cutW, cutC, cutAll = [], [], []
            for position in dataHash[enzyme].split(':'):
                beg, end = map(int, position.split(','))
                if beg < end:  # watson
                    cutSite = beg + int(offset[enzyme]) - 1
                    if cutSite not in cutW:
                        cutW.append(cutSite)
                else:          # crick
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
                size = cs - preCutSite
                if size and size not in seen:
                    cutFragments.append(size)
                    seen.add(size)
                preCutSite = cs

            cutSiteW = ", ".join(map(str, sorted(cutW)))
            cutSiteC = ", ".join(map(str, sorted(cutC)))
            fragmentsReal = ", ".join(map(str, cutFragments))
            fragments = ", ".join(map(str, sorted(cutFragments, reverse=True)))
            cutNum = max(0, len(cutFragments) - 1)

            fw.write(f"{enzyme}\t{offset[enzyme]}\t{overhang[enzyme]}\t{recognition_seq[enzyme]}\t{enzyme_type.get(enzyme,'')}\t{cutNum}\t{fragmentsReal}\t{fragments}\t{cutSiteW}\t{cutSiteC}\n")

            data[enzyme] = {
                "cut_site_on_watson_strand": cutSiteW,
                "cut_site_on_crick_strand":  cutSiteC,
                "fragment_size":              fragments,
                "fragment_size_in_real_order":fragmentsReal,
                "offset":                     offset[enzyme],
                "overhang":                   overhang[enzyme],
                "recognition_seq":            recognition_seq[enzyme],
                "enzyme_type":                enzyme_type.get(enzyme, ''),
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

    enzymetype = _normalize_type(f.get('type') or p.get('type') or 'all')
    
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
        (data, notCutEnzymeList) = process_data(
            seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut, enzymefile
        )
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
    
