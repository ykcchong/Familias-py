"""Pydantic v2 models for the Familias web REST API."""
from __future__ import annotations
from typing import Dict, List, Literal, Optional, Tuple

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
class PersonModel(BaseModel):
    id: str
    sex: Literal["male", "female"]


class RelationModel(BaseModel):
    parent: str
    child: str
    flag: Literal["fixed", "test"] = "fixed"


class GenotypeModel(BaseModel):
    """A single (person, locus) genotype: two alleles."""
    person: str
    a1: str
    a2: str


class LocusInputModel(BaseModel):
    """Frequencies + observed genotypes at a single locus."""
    name: str
    frequencies: Dict[str, float] = Field(
        default_factory=dict,
        description="{allele: freq}. May contain a 'Rest' bin.",
    )
    genotypes: List[GenotypeModel] = Field(default_factory=list)


class CaseBase(BaseModel):
    mutation_model: str = "Equal"
    mutation_rate: float = 0.001
    kinship: float = 0.0
    prior: float = Field(0.5, ge=0.0, le=1.0,
                         description="Prior P(H1). Used only for posterior.")
    loci: List[LocusInputModel]


class OneParentCase(CaseBase):
    """AF (alleged father) + CH (child)."""
    pass


class TwoParentCase(CaseBase):
    """MO (known mother) + AF (alleged father) + CH (child)."""
    pass


class ArbitraryCase(CaseBase):
    persons: List[PersonModel]
    relations: List[RelationModel]


# ---------------------------------------------------------------------------
# Single-locus quick LR
# ---------------------------------------------------------------------------
class SingleLocusRequest(BaseModel):
    mode: Literal["one-parent", "two-parent", "arbitrary"]
    locus: str
    frequencies: Dict[str, float]
    genotypes: List[GenotypeModel]
    mutation_model: str = "Equal"
    mutation_rate: float = 0.001
    kinship: float = 0.0
    # only used for arbitrary mode:
    persons: Optional[List[PersonModel]] = None
    relations: Optional[List[RelationModel]] = None


class SingleLocusResponse(BaseModel):
    LR: Optional[float]
    L1: Optional[float] = None
    L2: Optional[float] = None
    reason: Optional[str] = None


# ---------------------------------------------------------------------------
# Mendelian check
# ---------------------------------------------------------------------------
class MendelianRequest(BaseModel):
    mode: Literal["one-parent", "two-parent"]
    parent: Optional[Tuple[str, str]] = None       # one-parent
    known: Optional[Tuple[str, str]] = None        # two-parent
    alleged: Optional[Tuple[str, str]] = None      # two-parent
    child: Optional[Tuple[str, str]] = None
    frequencies: Optional[Dict[str, float]] = None  # for PE


class MendelianResponse(BaseModel):
    mismatch: bool
    power_of_exclusion: Optional[float] = None


# ---------------------------------------------------------------------------
# Compute response
# ---------------------------------------------------------------------------
class PerLocusResult(BaseModel):
    locus: str
    L1: float
    L2: float
    LR: float


class ComputeResponse(BaseModel):
    LR: float
    log10_LR: float
    likelihoods: Tuple[float, float]
    posterior: Tuple[float, float]
    posterior_h1: float = Field(
        ...,
        description="Posterior P(H1) using the supplied prior (overrides "
                    "the engine's uniform prior).",
    )
    per_locus: List[PerLocusResult]
    warnings: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Frequency DB listing
# ---------------------------------------------------------------------------
class DatabaseListResponse(BaseModel):
    databases: List[str]
    default: str


class DatabaseLociResponse(BaseModel):
    database: str
    loci: List[str]


class DatabaseFullResponse(BaseModel):
    database: str
    loci: Dict[str, Dict[str, float]]
