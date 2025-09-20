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
    et = (enzymetype or '').strip().lower()
    if not et or et == 'all':
        return dataDir + 'rest_enzymes'
    if 'six-base' in et:
        return dataDir + 'rest_enzymes.6base'
    if 'blunt' in et:
        return dataDir + 'rest_enzymes.blunt'
    if et.startswith("3"):
        return dataDir + 'rest_enzymes.3'
    if et.startswith("5"):
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
    # normalize once
    etype = (enzymetype or '').strip().lower()

    dataHash = {}
    offset = {}
    overhang = {}
    recognition_seq = {}
    notCutEnzyme = []
    with open(outfile, encoding="utf-8") as f:
        preLine = ''
        enzyme = ''
        for line in f:
            if line.startswith('>>'):
                pieces = line.strip().split(' ')
                enzyme = pieces[0].replace('>>', '').replace(':', '')
                offset[enzyme] = pieces[1]
                overhang[enzyme] = pieces[2]
                recognition_seq[enzyme] = pieces[3]
                # mark previous enzyme as "no cut" if we saw a header followed by another header
                if etype in ('all', '') or etype.startswith('enzymes that do not'):
                    if preLine.startswith('>>'):
                        p2 = preLine.replace('>>', '').replace(':', '').split(' ')
                        if p2[0] not in notCutEnzyme:
                            notCutEnzyme.append(p2[0])
            elif line.startswith('>'):
                # capture a cut coordinate for the current enzyme
                coords = line.strip().split(':')[1].replace('[', '').replace(']', '')
                if enzyme in dataHash:
                    dataHash[enzyme] = dataHash[enzyme] + ':' + coords
                else:
                    dataHash[enzyme] = coords
            preLine = line.strip()

    # also check the last header in file (could be "no cut")
    if etype in ('all', '') or etype.startswith('enzymes that do not'):
        if preLine.startswith('>>'):
            p2 = preLine.replace('>>', '').replace(':', '').split(' ')
            if p2[0] not in notCutEnzyme:
                notCutEnzyme.append(p2[0])

    # write "do not cut" list (always)
    with open(downloadfile4notCut, 'w') as fw:
        for e in sorted(notCutEnzyme):
            fw.write(e + "\n")
    # if the request is specifically "enzymes that do not cut", return early
    if etype.startswith('enzymes that do not'):
        return ({}, notCutEnzyme)

    # handle "cut once/twice" filters
    if 'cut' in etype:
        cutLimit = 2 if 'twice' in etype else 1
        newDataHash = {}
        for key, coords_str in dataHash.items():
            if not coords_str:
                continue
            wCut = cCut = 0
            for coordPair in coords_str.split(':'):
                beg, end = coordPair.split(',')
                beg = int(beg); end = int(end)
                if beg < end:
                    wCut += 1
                else:
                    cCut += 1
            if (cCut == cutLimit and wCut <= cutLimit) or (wCut == cutLimit and cCut <= cutLimit):
                newDataHash[key] = coords_str
        dataHash = newDataHash

    # map enzyme -> type
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

            # subset filters like "3' overhang", "5' overhang", "blunt end"
            if (('overhang' in etype) or ('blunt' in etype)) and et_label.lower() != etype:
                continue

            cutW, cutC, cutAll = [], [], []
            for position in dataHash[enzyme].split(':'):
                beg, end = position.split(',')
                beg = int(beg); end = int(end)
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
            for cutSite in sorted(cutAll):
                cutSize = cutSite - preCutSite
                if cutSize and cutSize not in seen:
                    cutFragments.append(cutSize)
                    seen.add(cutSize)
                preCutSite = cutSite

            cutSiteW = ", ".join(str(x) for x in sorted(cutW))
            cutSiteC = ", ".join(str(x) for x in sorted(cutC))
            fragmentsReal = ", ".join(str(x) for x in cutFragments)
            fragments = ", ".join(str(x) for x in sorted(cutFragments, reverse=True))
            cutNum = len(cutFragments) - 1
            fw.write(f"{enzyme}\t{offset[enzyme]}\t{overhang[enzyme]}\t{recognition_seq[enzyme]}\t{et_label}\t{cutNum}\t{fragmentsReal}\t{fragments}\t{cutSiteW}\t{cutSiteC}\n")

            data[enzyme] = {
                "cut_site_on_watson_strand": cutSiteW,
                "cut_site_on_crick_strand": cutSiteC,
                "fragment_size": fragments,
                "fragment_size_in_real_order": fragmentsReal,
                "offset": offset[enzyme],
                "overhang": overhang[enzyme],
                "recognition_seq": recognition_seq[enzyme],
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
    
