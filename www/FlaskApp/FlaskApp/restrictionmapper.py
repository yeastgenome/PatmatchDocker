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
    t = (enzymetype or '').strip().lower()
    if 'six-base' in t:
        return dataDir + 'rest_enzymes.6base'
    if 'blunt' in t:
        return dataDir + 'rest_enzymes.blunt'
    if t.startswith('3') or "3'" in t:
        return dataDir + 'rest_enzymes.3'
    if t.startswith('5') or "5'" in t:
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
    # ---- normalize & mode flags ----
    t_raw = (enzymetype or '').replace('+', ' ').replace("%27", "'").strip()
    t = t_raw.lower()
    want_not_cut_only = (t == 'all' or t == '' or (t.startswith('enzymes') and 'not' in t))

    dataHash = {}
    offset = {}
    overhang = {}
    recognition_seq = {}
    notCutEnzyme = []

    current = None
    had_hit = False


    with open(outfile, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s.startswith('>>'):
                # finalize previous enzyme (no hits)
                if current is not None and not had_hit and want_not_cut_only:
                    notCutEnzyme.append(current)
                # parse header
                parts = s.split()
                current = parts[0].replace('>>', '').replace(':', '')
                offset[current] = parts[1]
                overhang[current] = parts[2]
                recognition_seq[current] = parts[3]
                had_hit = False
            elif s.startswith('>'):
                # hit line from scan_for_matches
                if current is None:
                    continue
                try:
                    coords = s.split(':', 1)[1].replace('[', '').replace(']', '')
                except IndexError:
                    continue
                dataHash[current] = f"{dataHash[current]}:{coords}" if current in dataHash else coords
                had_hit = True
            else:
                continue

    # tail finalize
    if current is not None and not had_hit and want_not_cut_only:
        notCutEnzyme.append(current)

    notCutEnzyme = sorted(set(notCutEnzyme))
    with open(downloadfile4notCut, 'w') as fw:
        for e in notCutEnzyme:
            fw.write(e + "\n")

    if want_not_cut_only:
        return ({}, notCutEnzyme)

    # ---- optional filters: 'cut once' / 'cut twice' ----
    if 'cut' in t:
        cutLimit = 2 if 'twice' in t else 1
        filtered = {}
        for key, coordstr in dataHash.items():
            wCut = cCut = 0
            for pair in coordstr.split(':'):
                beg, end = map(int, pair.split(','))
                if beg < end: wCut += 1
                else:         cCut += 1
            if ((cCut == cutLimit and wCut <= cutLimit) or
                (wCut == cutLimit and cCut <= cutLimit)):
                filtered[key] = coordstr
        dataHash = filtered

    # ---- enzyme types for display / filtering ----
    enzyme_type = {}
    set_enzyme_types(enzyme_type, "3' overhang")
    set_enzyme_types(enzyme_type, "5' overhang")
    set_enzyme_types(enzyme_type, "blunt end")

    data = {}
    with open(downloadfile4cutSite, 'w') as fw:
        fw.write("Enzyme\toffset (bp)\toverhang (bp)\trecognition sequence\tenzyme type\tnumber of cuts\tordered fragment size\tsorted fragment size\tcut site on watson strand\tcut site on crick strand\n")

        for enzyme in sorted(dataHash):
            # honor explicit overhang/blunt selections
            if (('overhang' in t or 'blunt' in t) and enzyme_type.get(enzyme) != t_raw):
                continue

            cutW, cutC, cutAll = [], [], []
            for position in dataHash[enzyme].split(':'):
                beg, end = map(int, position.split(','))
                if beg < end:
                    cutSite = beg + int(offset[enzyme]) - 1
                    if cutSite not in cutW: cutW.append(cutSite)
                else:
                    beg, end = end, beg
                    cutSite = beg + int(offset[enzyme]) + int(overhang[enzyme]) - 1
                    if cutSite not in cutC: cutC.append(cutSite)
                if cutSite not in cutAll: cutAll.append(cutSite)

            cutAll.append(seqLen)

            pre = 0
            seen = set()
            cutFragments = []
            for c in sorted(cutAll):
                size = c - pre
                if size and size not in seen:
                    cutFragments.append(size); seen.add(size)
                pre = c

            cutSiteW = ", ".join(str(x) for x in sorted(cutW))
            cutSiteC = ", ".join(str(x) for x in sorted(cutC))
            fragmentsReal = ", ".join(str(x) for x in cutFragments)
            fragments = ", ".join(str(x) for x in sorted(cutFragments, reverse=True))
            cutNum = max(len(cutFragments) - 1, 0)

            fw.write(f"{enzyme}\t{offset[enzyme]}\t{overhang[enzyme]}\t{recognition_seq[enzyme]}\t{enzyme_type.get(enzyme,'')}\t{cutNum}\t{fragmentsReal}\t{fragments}\t{cutSiteW}\t{cutSiteC}\n")

            data[enzyme] = {
                "cut_site_on_watson_strand": cutSiteW,
                "cut_site_on_crick_strand":  cutSiteC,
                "fragment_size":              fragments,
                "fragment_size_in_real_order": fragmentsReal,
                "offset":                      offset[enzyme],
                "overhang":                    overhang[enzyme],
                "recognition_seq":             recognition_seq[enzyme],
                "enzyme_type":                 enzyme_type.get(enzyme, "")
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

    enzymetype = (f.get('type') or p.get('type') or 'ALL')
    enzymetype = enzymetype.replace('+', ' ').replace("%27", "'").strip()
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
    
