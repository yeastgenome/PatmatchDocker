"""Pydantic models for PatmatchDocker FastAPI application."""
from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field


class PatmatchRequest(BaseModel):
    """Request model for pattern matching."""
    pattern: str
    dataset: Optional[str] = None
    seqtype: str = "pep"
    strand: Optional[str] = None
    insertion: Optional[str] = None
    deletion: Optional[str] = None
    substitution: Optional[str] = None
    mismatch: int = 0
    max_hits: int = 500


class PatmatchHit(BaseModel):
    """A single pattern match hit."""
    seqname: str
    beg: int
    end: int
    count: int
    matchingPattern: str
    gene_name: Optional[str] = None
    sgdid: Optional[str] = None
    desc: Optional[str] = None
    # Additional fields for NotFeature dataset
    orfs: Optional[str] = None
    chr: Optional[str] = None


class PatmatchResponse(BaseModel):
    """Response model for pattern matching."""
    hits: List[Dict[str, Any]]
    uniqueHits: int
    totalHits: int
    downloadUrl: str = ""
    error_message: Optional[str] = None


class SequenceResponse(BaseModel):
    """Response model for sequence retrieval."""
    defline: str
    seq: str


class RestrictionEnzymeData(BaseModel):
    """Data for a single restriction enzyme."""
    cut_site_on_watson_strand: str
    cut_site_on_crick_strand: str
    fragment_size: str
    fragment_size_in_real_order: str
    offset: str
    overhang: str
    recognition_seq: str
    enzyme_type: str


class RestrictionMapperResponse(BaseModel):
    """Response model for restriction mapper."""
    data: Dict[str, RestrictionEnzymeData]
    seqName: str
    chrCoords: str
    seqLength: int
    notCutEnzyme: List[str]
    downloadUrl: str = ""
    downloadUrl4notCutEnzyme: str = ""
    ERROR: Optional[str] = None


class ConfigGenome(BaseModel):
    """Genome configuration entry."""
    strain: str
    label: str


class ConfigDataset(BaseModel):
    """Dataset configuration entry."""
    dataset_file_name: str
    seqtype: str
    label: str
    seqcount: str


class PatmatchConfig(BaseModel):
    """Full patmatch configuration."""
    genome: List[ConfigGenome]
    dataset: Dict[str, List[ConfigDataset]]
