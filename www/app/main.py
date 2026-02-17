"""
PatmatchDocker FastAPI Application

A FastAPI-based service for pattern matching in biological sequences.
"""
import os
import re
import random
import hashlib
import threading
from contextlib import asynccontextmanager
from typing import Optional

from fastapi import FastAPI, Request, Query, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse

from patmatch_service import (
    run_patmatch,
    get_config,
    get_download_url,
    DATA_DIR,
    TMP_DIR,
    CONF_DIR
)
from restriction_service import (
    RestrictionService,
    run_restriction_search,
    get_enzyme_types
)
from fasta_utils import get_sequence

# Try to import boto3
try:
    import boto3
    HAS_BOTO3 = True
except ImportError:
    boto3 = None
    HAS_BOTO3 = False


# Configuration loaded at startup
_config = None
RANDOM_MAX = 10000000

# Restriction mapper directories
RESTRICTION_DATA_DIR = os.environ.get('RESTRICTION_DATA_DIR', '/data/restriction_mapper/')


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler - load config at startup."""
    global _config
    try:
        _config = get_config('patmatch')
    except Exception as e:
        print(f"Warning: Could not load config: {e}")
        _config = {}
    yield


app = FastAPI(
    title="PatmatchDocker",
    description="Pattern matching service for biological sequences",
    version="2.0.0",
    lifespan=lifespan
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


def get_request_id() -> str:
    """Generate a unique request ID."""
    return str(random.randint(1, RANDOM_MAX))


def get_param(request: Request, name: str, default: Optional[str] = None) -> Optional[str]:
    """Get a parameter from query string or form data."""
    value = request.query_params.get(name)
    if value is None and hasattr(request, '_form'):
        value = request._form.get(name)
    return value if value is not None else default


@app.get("/")
async def hello():
    """Health check endpoint."""
    return "Hello, we all love SGD!!"


@app.get("/patmatch")
@app.post("/patmatch")
async def patmatch(
    request: Request,
    conf: Optional[str] = Query(None),
    file: Optional[str] = Query(None),
    seqname: Optional[str] = Query(None),
    dataset: Optional[str] = Query(None),
    pattern: Optional[str] = Query(None),
    seqtype: str = Query('pep'),
    strand: Optional[str] = Query(None),
    insertion: Optional[str] = Query(None),
    deletion: Optional[str] = Query(None),
    substitution: Optional[str] = Query(None),
    mismatch: int = Query(0),
    max_hits: int = Query(500)
):
    """
    Pattern matching endpoint.

    Supports both GET and POST methods with the following parameters:
    - conf: Configuration name to retrieve
    - file: File to download
    - seqname: Specific sequence name to retrieve
    - dataset: Dataset name
    - pattern: Search pattern
    - seqtype: Sequence type ('pep', 'protein', 'dna', 'nuc')
    - strand: Strand option for DNA
    - insertion: Allow insertions
    - deletion: Allow deletions
    - substitution: Allow substitutions
    - mismatch: Number of allowed mismatches
    - max_hits: Maximum number of hits
    """
    # Handle form data for POST
    if request.method == "POST":
        form = await request.form()
        conf = form.get('conf', conf)
        file = form.get('file', file)
        seqname = form.get('seqname', seqname)
        dataset = form.get('dataset', dataset)
        pattern = form.get('pattern', pattern)
        seqtype = form.get('seqtype', seqtype) or 'pep'
        strand = form.get('strand', strand)
        insertion = form.get('insertion', insertion)
        deletion = form.get('deletion', deletion)
        substitution = form.get('substitution', substitution)
        mismatch = int(form.get('mismatch', mismatch) or 0)
        max_hits = int(form.get('max_hits', max_hits) or 500)

    # Return configuration
    if conf:
        config_data = get_config(conf)
        return JSONResponse(content=config_data)

    # Return file download
    if file:
        file_path = os.path.join(TMP_DIR, file)
        if os.path.exists(file_path):
            return FileResponse(
                file_path,
                filename=file,
                media_type='application/text'
            )
        return JSONResponse(content={"error": "File not found"}, status_code=404)

    # Return sequence
    if seqname and dataset:
        data = get_sequence(dataset, seqname, DATA_DIR)
        return JSONResponse(content=data)

    # Run pattern search
    if not pattern:
        return JSONResponse(content={"error": "Pattern is required"}, status_code=400)

    result = run_patmatch(
        pattern=pattern,
        dataset=dataset,
        seqtype=seqtype,
        seqname=seqname,
        strand=strand,
        insertion=insertion,
        deletion=deletion,
        substitution=substitution,
        mismatch=mismatch,
        max_hits=max_hits,
        request_id=get_request_id()
    )

    return JSONResponse(content=result)


@app.get("/restrictionmapper")
@app.post("/restrictionmapper")
async def restrictionmapper(
    request: Request,
    file: Optional[str] = Query(None),
    seq: Optional[str] = Query(None),
    name: Optional[str] = Query(None),
    type: Optional[str] = Query('ALL')
):
    """
    Restriction enzyme mapping endpoint.

    Supports both GET and POST methods with the following parameters:
    - file: File to download
    - seq: DNA sequence to analyze
    - name: Gene/sequence name to look up
    - type: Enzyme type filter
    """
    # Handle form data for POST
    if request.method == "POST":
        form = await request.form()
        file = form.get('file', file)
        seq = form.get('seq', seq)
        name = form.get('name', name)
        type = form.get('type', type) or 'ALL'

    # Return file download
    if file:
        file_path = os.path.join(TMP_DIR, file)
        if os.path.exists(file_path):
            return FileResponse(
                file_path,
                filename=file,
                media_type='application/text'
            )
        return JSONResponse(content={"error": "File not found"}, status_code=404)

    request_id = get_request_id()
    pat_file = os.path.join(TMP_DIR, f"patfile.{request_id}.txt")
    out_file = os.path.join(TMP_DIR, f"outfile.{request_id}.txt")
    seq_file = os.path.join(TMP_DIR, f"seqfile.{request_id}.txt")

    cut_site_file = f"restrictionmapper.{request_id}"
    not_cut_file = f"restrictionmapper_not_cut_enzyme.{request_id}"

    download_file_cut_site = os.path.join(TMP_DIR, cut_site_file)
    download_file_not_cut = os.path.join(TMP_DIR, not_cut_file)

    # Clean up enzyme type parameter
    enzyme_type = type.replace('+', ' ').replace("%27", "'")

    if enzyme_type.startswith('3'):
        enzyme_type = "3' overhang"
    elif enzyme_type.startswith('5'):
        enzyme_type = "5' overhang"

    # Get sequence
    service = RestrictionService(RESTRICTION_DATA_DIR)
    defline = None
    sequence = seq

    if sequence:
        defline = ">Unnamed sequence"
    else:
        defline, sequence = service.get_sequence_from_fasta(name)

    # Write sequence file and get info
    seq_name, chr_coords, seq_len = service.write_seq_file(defline, sequence, seq_file)

    # Run search
    try:
        data, not_cut_enzymes, err = run_restriction_search(
            sequence=sequence,
            enzyme_type=enzyme_type,
            data_dir=RESTRICTION_DATA_DIR
        )

        # Write output files
        write_cut_site_file(data, download_file_cut_site)
        write_not_cut_file(not_cut_enzymes, download_file_not_cut)

        # Get download URLs
        download_url_cut_site = ''
        download_url_not_cut = ''

        try:
            if os.path.exists(download_file_cut_site):
                download_url_cut_site = get_download_url(cut_site_file)
            if os.path.exists(download_file_not_cut):
                download_url_not_cut = get_download_url(not_cut_file)
        except Exception as e:
            err = f"Error generating download URL: {e}"

        return JSONResponse(content={
            "data": data,
            "seqName": seq_name,
            "chrCoords": chr_coords,
            "seqLength": seq_len,
            "notCutEnzyme": not_cut_enzymes,
            "downloadUrl": download_url_cut_site,
            "downloadUrl4notCutEnzyme": download_url_not_cut
        })

    except Exception as e:
        return JSONResponse(content={
            "ERROR": str(e),
            "seqName": seq_name,
            "chrCoords": chr_coords,
            "seqLength": seq_len,
            "notCutEnzyme": [],
            "downloadUrl": '',
            "downloadUrl4notCutEnzyme": ''
        })


def write_cut_site_file(data: dict, file_path: str):
    """Write cut site data to file."""
    with open(file_path, 'w') as fw:
        header = ("Enzyme\toffset (bp)\toverhang (bp)\trecognition sequence\t"
                  "enzyme type\tnumber of cuts\tordered fragment size\t"
                  "sorted fragment size\tcut site on watson strand\t"
                  "cut site on crick strand\n")
        fw.write(header)

        for enzyme in sorted(data.keys()):
            info = data[enzyme]
            cut_num = len(info.get('fragment_size_in_real_order', '').split(', ')) - 1
            line = (f"{enzyme}\t{info.get('offset', '')}\t{info.get('overhang', '')}\t"
                    f"{info.get('recognition_seq', '')}\t{info.get('enzyme_type', '')}\t"
                    f"{cut_num}\t{info.get('fragment_size_in_real_order', '')}\t"
                    f"{info.get('fragment_size', '')}\t"
                    f"{info.get('cut_site_on_watson_strand', '')}\t"
                    f"{info.get('cut_site_on_crick_strand', '')}\n")
            fw.write(line)


def write_not_cut_file(enzymes: list, file_path: str):
    """Write non-cutting enzymes to file."""
    with open(file_path, 'w') as fw:
        for enzyme in sorted(enzymes):
            fw.write(enzyme + '\n')


@app.get("/download/{filename}")
async def download_file(filename: str):
    """
    Download a result file.

    Args:
        filename: Name of the file to download
    """
    file_path = os.path.join(TMP_DIR, filename)
    if os.path.exists(file_path):
        return FileResponse(
            file_path,
            filename=filename,
            media_type='application/text'
        )
    return JSONResponse(content={"error": "File not found"}, status_code=404)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
