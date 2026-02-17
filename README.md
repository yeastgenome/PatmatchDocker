# PatmatchDocker

A FastAPI-based microservice for pattern matching in biological sequences. Supports protein and nucleotide pattern searches with fuzzy matching capabilities, as well as restriction enzyme mapping.

## Features

- **Pattern Matching**: Search for patterns in protein and DNA sequences using Patmatch syntax
- **Fuzzy Matching**: Support for insertions, deletions, and substitutions with configurable error tolerance
- **IUPAC Codes**: Full support for IUPAC ambiguity codes for both nucleotides and amino acids
- **Restriction Mapping**: Find restriction enzyme cut sites and calculate fragment sizes
- **S3 Integration**: Optional upload of result files to AWS S3
- **Multiple Datasets**: Support for multiple yeast strain genomes

## Quick Start

### Using Docker

```bash
# Build the image
docker build -t patmatch-fastapi .

# Run with data volumes
docker run -d \
  -p 8000:8000 \
  -v /path/to/patmatch/data:/data/patmatch:ro \
  -v /path/to/restriction_mapper/data:/data/restriction_mapper:ro \
  patmatch-fastapi
```

### Using Docker Compose

```bash
# Set data paths and run
DATA_DIR=/path/to/patmatch/data \
RESTRICTION_DATA_DIR=/path/to/restriction_mapper/data \
docker-compose up -d
```

### Local Development

```bash
# Install dependencies
pip install -r requirements.txt

# Run the application
cd www/app
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

## API Endpoints

### Health Check

```
GET /
```

Returns: `"Hello, we all love SGD!!"`

### Pattern Matching

```
GET/POST /patmatch
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pattern` | string | required | Search pattern (Patmatch syntax) |
| `dataset` | string | auto | Dataset name (e.g., `orf_pep`, `orf_dna`) |
| `seqtype` | string | `pep` | Sequence type: `pep`, `protein`, `dna`, `nuc` |
| `strand` | string | null | Strand option: `Both`, `Watson`, `Crick` |
| `insertion` | string | null | Allow insertions: `insertion` |
| `deletion` | string | null | Allow deletions: `deletion` |
| `substitution` | string | null | Allow substitutions: `substitution` |
| `mismatch` | int | 0 | Maximum number of mismatches |
| `max_hits` | int | 500 | Maximum hits to return |
| `conf` | string | null | Return configuration instead of searching |
| `seqname` | string | null | Return specific sequence |

**Example:**

```bash
# Search for pattern in protein sequences
curl "http://localhost:8000/patmatch?pattern=MFVL&dataset=orf_pep&seqtype=pep"

# Search with fuzzy matching (1 mismatch allowed)
curl "http://localhost:8000/patmatch?pattern=ATGCATGC&dataset=orf_dna&seqtype=dna&mismatch=1"

# Get configuration
curl "http://localhost:8000/patmatch?conf=patmatch"
```

**Response:**

```json
{
  "hits": [
    {
      "seqname": "YAL001C",
      "beg": 1,
      "end": 4,
      "count": 1,
      "matchingPattern": "MFVL",
      "gene_name": "TFC3",
      "sgdid": "S000000001",
      "desc": "Subunit of RNA polymerase III..."
    }
  ],
  "uniqueHits": 1,
  "totalHits": 1,
  "downloadUrl": "https://bucket.s3.amazonaws.com/patmatch/abc123.txt",
  "error_message": null
}
```

### Restriction Mapping

```
GET/POST /restrictionmapper
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `seq` | string | null | DNA sequence to analyze |
| `name` | string | null | Gene name to look up sequence |
| `type` | string | `ALL` | Enzyme type filter |

**Example:**

```bash
# Analyze a sequence
curl "http://localhost:8000/restrictionmapper?seq=ATGCGAATTCATGC"

# Look up gene and analyze
curl "http://localhost:8000/restrictionmapper?name=ACT1"
```

**Response:**

```json
{
  "data": {
    "EcoRI": {
      "cut_site_on_watson_strand": "5",
      "cut_site_on_crick_strand": "9",
      "fragment_size": "14, 5",
      "fragment_size_in_real_order": "5, 9",
      "offset": "1",
      "overhang": "4",
      "recognition_seq": "GAATTC",
      "enzyme_type": "5' overhang"
    }
  },
  "seqName": "ACT1/YFL039C",
  "chrCoords": "Chr VI from 54696-53560",
  "seqLength": 1137,
  "notCutEnzyme": ["AatII", "AflII", ...],
  "downloadUrl": "https://...",
  "downloadUrl4notCutEnzyme": "https://..."
}
```

### File Download

```
GET /download/{filename}
```

Download result files from the server.

## Pattern Syntax

### Basic Patterns

- Single characters: `MFVL` matches exactly "MFVL"
- Wildcards: `X` or `.` matches any character
- Character classes: `[AILV]` matches any of A, I, L, or V
- Exclusions: `[^P]` matches any character except P

### Repetitions

- Exact: `A{3}` matches "AAA"
- Range: `A{2,4}` matches "AA", "AAA", or "AAAA"
- Minimum: `A{2,}` matches "AA" or more
- Maximum: `A{,3}` matches up to "AAA"

### IUPAC Codes (Nucleotides)

| Code | Bases | Description |
|------|-------|-------------|
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| S | G, C | Strong |
| W | A, T | Weak |
| M | A, C | Amino |
| K | G, T | Keto |
| V | A, C, G | Not T |
| H | A, C, T | Not G |
| D | A, G, T | Not C |
| B | C, G, T | Not A |
| N | A, C, G, T | Any |

### IUPAC Codes (Protein)

| Code | Amino Acids | Description |
|------|-------------|-------------|
| X | Any | Any amino acid |
| B | D, N | Aspartic acid or Asparagine |
| Z | E, Q | Glutamic acid or Glutamine |
| J | I, F, V, L, W, M, A, G, C, Y | Hydrophobic |
| O | T, S, H, E, D, Q, N, K, R | Hydrophilic |

### Anchors

- `<pattern` - Match at start of sequence
- `pattern>` - Match at end of sequence

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `DATA_DIR` | `/data/patmatch/` | Path to patmatch data files |
| `RESTRICTION_DATA_DIR` | `/data/restriction_mapper/` | Path to restriction enzyme data |
| `TMP_DIR` | `/var/www/tmp/` | Temporary file directory |
| `CONF_DIR` | `/var/www/conf/` | Configuration file directory |
| `S3_BUCKET` | (empty) | S3 bucket for result uploads |
| `AWS_ACCESS_KEY_ID` | (empty) | AWS credentials |
| `AWS_SECRET_ACCESS_KEY` | (empty) | AWS credentials |
| `AWS_DEFAULT_REGION` | `us-east-1` | AWS region |

## Data Files

### Patmatch Data (`/data/patmatch/`)

- `orf_pep.seq` - Protein sequences (FASTA format)
- `orf_dna.seq` - DNA coding sequences (FASTA format)
- `orf_genomic.seq` - Genomic DNA sequences (FASTA format)
- `locus.txt` - Gene annotation data (TSV format)
- Additional strain-specific datasets

### Restriction Mapper Data (`/data/restriction_mapper/`)

- `orf_genomic.seq` - Genomic sequences for lookup
- `rest_enzymes` - All restriction enzymes
- `rest_enzymes.6base` - Six-base cutters
- `rest_enzymes.blunt` - Blunt end cutters
- `rest_enzymes.3` - 3' overhang enzymes
- `rest_enzymes.5` - 5' overhang enzymes

## Project Structure

```
PatmatchDocker/
├── Dockerfile              # Docker image definition
├── docker-compose.yml      # Docker Compose configuration
├── requirements.txt        # Python dependencies
├── README.md               # This file
└── www/
    ├── app/                # FastAPI application
    │   ├── main.py         # Application entry point
    │   ├── schemas.py      # Pydantic models
    │   ├── patmatch_service.py    # Pattern matching service
    │   ├── pattern_converter.py   # Pattern syntax conversion
    │   ├── fuzzy_matcher.py       # Fuzzy matching engine
    │   ├── fasta_utils.py         # FASTA file utilities
    │   └── restriction_service.py # Restriction enzyme service
    └── conf/
        └── patmatch.json   # Dataset configuration
```

## Development

### Running Tests

```bash
# Install dev dependencies
pip install pytest httpx

# Run tests
pytest tests/
```

### Building for Production

```bash
# Build optimized image
docker build -t patmatch-fastapi:latest .

# Tag for registry
docker tag patmatch-fastapi:latest your-registry/patmatch-fastapi:latest

# Push to registry
docker push your-registry/patmatch-fastapi:latest
```

## API Compatibility

This service maintains backward compatibility with the legacy Flask-based API. All endpoints accept the same parameters and return the same response format.

## License

Copyright (c) Stanford University / SGD Project

## Related Projects

- [SGD (Saccharomyces Genome Database)](https://www.yeastgenome.org/)
- [Patmatch at TAIR](https://www.arabidopsis.org/cgi-bin/patmatch/nph-patmatch.pl)
