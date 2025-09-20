import json
import os
import re
from patmatch import get_downloadUrl

binDir = '/var/www/bin/'
dataDir = '/data/restriction_mapper/'
tmpDir = "/var/www/tmp/"

scan4matches = binDir + "scan_for_matches"
fastafile = dataDir + "orf_genomic.seq"

    
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

def _dump_debug(base_path: str, name: str, payload: dict):
    """Write debug JSON next to the cutSite file, never raise."""
    try:
        with open(f"{base_path}.{name}.json", "w", encoding="utf-8") as df:
            json.dump(payload, df, indent=2, ensure_ascii=False)
    except Exception:
        pass

def process_data(seqLen, enzymetype, outfile, downloadfile4cutSite, downloadfile4notCut):
    # --- Debug controls (set via env in ECS task def or Dockerfile) ---
    DEBUG = os.environ.get('RESTMAP_DEBUG', '0') == '1'
    CP    = (os.environ.get('RESTMAP_CP') or '').upper()  # PARSE | WROTE_NO_CUT | EARLY_RETURN

    etype   = (enzymetype or '').strip()
    etype_l = etype.lower()

    dataHash, offset, overhang, recognition_seq = {}, {}, {}, {}
    all_enzymes = []
    enzymes_with_cuts = set()

    # --- Parse scan output ---
    with open(outfile, encoding="utf-8") as f:
        current = None
        for raw in f:
            s = raw.strip()
            if s.startswith('>>'):
                parts = s.split(' ')
                enzyme = parts[0].replace('>>', '').replace(':', '')
                current = enzyme
                all_enzymes.append(enzyme)
                if len(parts) >= 4:
                    offset[enzyme], overhang[enzyme], recognition_seq[enzyme] = parts[1], parts[2], parts[3]
                else:
                    offset[enzyme] = offset.get(enzyme, '0')
                    overhang[enzyme] = overhang.get(enzyme, '0')
                    recognition_seq[enzyme] = recognition_seq.get(enzyme, 'N/A')
            elif s.startswith('>') and current:
                after = s.split(':', 1)
                if len(after) == 2:
                    coords = after[1].replace('[', '').replace(']', '')
                    dataHash[current] = (dataHash.get(current, '') + (':' if current in dataHash else '') + coords)
                    enzymes_with_cuts.add(current)

    notCutEnzyme = sorted([e for e in all_enzymes if e not in enzymes_with_cuts])
    # notCutEnzyme = sorted(set(all_enzymes) - enzymes_with_cuts)
    
    # --- Debug: after parse ---
    if DEBUG and CP == 'PARSE':
        _dump_debug(downloadfile4cutSite, "parse", {
            "etype": etype,
            "num_headers": len(all_enzymes),
            "num_with_cuts": len(enzymes_with_cuts),
            "num_no_cut": len(notCutEnzyme),
            "sample_headers": all_enzymes[:10],
            "sample_no_cut": notCutEnzyme[:10]
        })

    # Always write the “do not cut” list
    with open(downloadfile4notCut, 'w') as fw:
        for e in notCutEnzyme:
            fw.write(e + "\n")

    # --- Debug: after writing no-cut file ---
    if DEBUG and CP == 'WROTE_NO_CUT':
        _dump_debug(downloadfile4cutSite, "wrote_no_cut", {
            "no_cut_count": len(notCutEnzyme),
            "file": downloadfile4notCut,
            "exists": os.path.isfile(downloadfile4notCut)
        })
    # --- Early return for "enzymes that do not ..." ---
    if etype_l.startswith('enzymes that do not'):
        if DEBUG and CP == 'EARLY_RETURN':
            _dump_debug(downloadfile4cutSite, "early_return", {
                "etype": etype,
                "no_cut_count": len(notCutEnzyme)
            })
        return ({}, notCutEnzyme)

    # --- Cut once/twice filter ---
    if 'cut' in etype_l:
        cutLimit = 2 if 'twice' in etype_l else 1
        newDataHash = {}
        for enzyme, coords_str in dataHash.items():
            if not coords_str: continue
            wCut = cCut = 0
            for pair in coords_str.split(':'):
                beg, end = map(int, pair.split(','))
                if beg < end: wCut += 1
                else:         cCut += 1
            if (cCut == cutLimit and wCut <= cutLimit) or (wCut == cutLimit and cCut <= cutLimit):
                newDataHash[enzyme] = coords_str
        dataHash = newDataHash
    # --- Enzyme type mapping (unchanged) ---
    enzyme_type = {}
    set_enzyme_types(enzyme_type, "3' overhang")
    set_enzyme_types(enzyme_type, "5' overhang")
    set_enzyme_types(enzyme_type, "blunt end")

    data = {}
    with open(downloadfile4cutSite, 'w') as fw:
        fw.write("Enzyme\toffset (bp)\toverhang (bp)\trecognition sequence\tenzyme type\tnumber of cuts\tordered fragment size\tsorted fragment size\tcut site on watson strand\tcut site on crick strand\n")
        for enzyme in sorted(dataHash):
            if not dataHash[enzyme]: continue
            et_label = enzyme_type.get(enzyme, 'unknown')
            if (('overhang' in etype_l) or ('blunt' in etype_l)) and et_label.lower() != etype_l:
                continue

            cutW, cutC, cutAll = [], [], []
            for pos in dataHash[enzyme].split(':'):
                beg, end = map(int, pos.split(','))
                if beg < end:
                    cutSite = beg + int(offset[enzyme]) - 1
                    if cutSite not in cutW: cutW.append(cutSite)
                else:
                    beg, end = end, beg
                    cutSite = beg + int(offset[enzyme]) + int(overhang[enzyme]) - 1
                    if cutSite not in cutC: cutC.append(cutSite)
                if cutSite not in cutAll: cutAll.append(cutSite)
            cutAll.append(seqLen)

            pre, seen, frags = 0, set(), []
            for cs in sorted(cutAll):
                size = cs - pre
                if size and size not in seen:
                    frags.append(size); seen.add(size)
                pre = cs

            cutSiteW = ", ".join(str(x) for x in sorted(cutW))
            cutSiteC = ", ".join(str(x) for x in sorted(cutC))
            fragmentsReal = ", ".join(str(x) for x in frags)
            fragments     = ", ".join(str(x) for x in sorted(frags, reverse=True))
            cutNum        = len(frags) - 1

            fw.write(f"{enzyme}\t{offset[enzyme]}\t{overhang[enzyme]}\t{recognition_seq[enzyme]}\t{et_label}\t{cutNum}\t{fragmentsReal}\t{fragments}\t{cutSiteW}\t{cutSiteC}\n")

            data[enzyme] = {
                "cut_site_on_watson_strand": cutSiteW,
                "cut_site_on_crick_strand":  cutSiteC,
                "fragment_size":              fragments,
                "fragment_size_in_real_order":fragmentsReal,
                "offset":                     offset[enzyme],
                "overhang":                   overhang[enzyme],
                "recognition_seq":            recognition_seq[enzyme],
                "enzyme_type":                et_label,
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
        return {
            "data": data,
            "seqName": seqNm,
            "chrCoords": chrCoords,
            "seqLength": seqLen,
            "notCutEnzyme": notCutEnzymeList,
            "downloadUrl": get_downloadUrl(cutSiteFile),
            "downloadUrl4notCutEnzyme": get_downloadUrl(notCutFile)
        }
    else:
        return {
            "ERROR": err,
            "seqName": seqNm,
            "chrCoords": chrCoords,
            "seqLength": seqLen,
            "notCutEnzyme": [],
            "downloadUrl": '',
            "downloadUrl4notCutEnzyme": ''
        }
    
