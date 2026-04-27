"""Glue between TUI inputs and the :mod:`familias` engine.

The functions in this module each take the user's textual inputs (DNA
table as TSV string, frequency database, mutation parameters, etc.) and
return a uniform ``Result`` dict:

    {"LR": float,                       # H1 / H2
     "log10_LR": float,
     "posterior": (p1, p2),
     "likelihoods": (L1, L2),
     "per_locus": [(locus, LR_i), ...],
     "warnings": [str, ...]}
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple
import math

import numpy as np
import pandas as pd

from .. import FamiliasLocus, FamiliasPedigree, FamiliasPosterior
from .defaults import canonical_locus_name


FreqDB = Dict[str, Dict[str, float]]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class Person:
    id: str
    sex: str  # "male" | "female"


@dataclass
class Relation:
    """A parent->child relationship.

    ``flag`` is either ``"fixed"`` (present in both pedigrees) or ``"test"``
    (present only in H1; absent in H2). For modes (1) and (2) the engine
    builds these automatically.
    """
    parent: str
    child: str
    flag: str = "fixed"

    def __post_init__(self):
        if self.flag not in ("fixed", "test"):
            raise ValueError(f"Bad flag {self.flag!r}: expected fixed|test.")


@dataclass
class TypedAllele:
    person: str
    locus: str
    a1: str
    a2: str


@dataclass
class Result:
    LR: float
    log10_LR: float
    posterior: Tuple[float, float]
    likelihoods: Tuple[float, float]
    per_locus: List[Tuple[str, float, float, float]]   # (locus, L1, L2, LR)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Parsers for TSV input fields
# ---------------------------------------------------------------------------
def parse_dna_table(text: str) -> List[TypedAllele]:
    """Parse a TSV / CSV / whitespace table of ``Person Locus A1 A2``.

    Lines starting with ``#`` and the literal header row are skipped.
    """
    out: List[TypedAllele] = []
    for raw in text.splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        parts = [p.strip() for p in s.replace(",", "\t").split("\t") if p.strip() != ""]
        if len(parts) < 4:
            parts = s.split()
        if len(parts) < 4:
            continue
        if parts[0].lower() == "person" and parts[1].lower() == "locus":
            continue
        out.append(TypedAllele(parts[0], parts[1], parts[2], parts[3]))
    return out


def parse_persons(text: str) -> List[Person]:
    out: List[Person] = []
    for raw in text.splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        parts = [p.strip() for p in s.replace(",", "\t").split("\t") if p.strip() != ""]
        if len(parts) < 2:
            parts = s.split()
        if len(parts) < 2:
            raise ValueError(f"Bad persons row: {raw!r}")
        if parts[0].lower() == "id":
            continue
        sex = parts[1].lower()
        if sex.startswith("m"):
            sex = "male"
        elif sex.startswith("f"):
            sex = "female"
        else:
            raise ValueError(f"Bad sex value {parts[1]!r} for {parts[0]!r}")
        out.append(Person(parts[0], sex))
    return out


def parse_relations(text: str) -> List[Relation]:
    out: List[Relation] = []
    for raw in text.splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        parts = [p.strip() for p in s.replace(",", "\t").split("\t") if p.strip() != ""]
        if len(parts) < 2:
            parts = s.split()
        if len(parts) < 2:
            raise ValueError(f"Bad relations row: {raw!r}")
        if parts[0].lower() == "parent":
            continue
        flag = parts[2].lower() if len(parts) >= 3 else "fixed"
        if flag.startswith("t"):
            flag = "test"
        else:
            flag = "fixed"
        out.append(Relation(parts[0], parts[1], flag))
    return out


# ---------------------------------------------------------------------------
# Locus / DNA assembly
# ---------------------------------------------------------------------------
def _resolve_freqs(locus_name: str, freq_db: FreqDB) -> Tuple[str, Dict[str, float]]:
    """Look a locus up in ``freq_db`` (with alias-tolerance)."""
    if locus_name in freq_db:
        return locus_name, freq_db[locus_name]
    canon = canonical_locus_name(locus_name)
    if canon in freq_db:
        return canon, freq_db[canon]
    # case-insensitive fallback
    lower = {k.lower(): k for k in freq_db}
    for cand in (locus_name, canon, locus_name.upper(), locus_name.lower()):
        if cand.lower() in lower:
            real = lower[cand.lower()]
            return real, freq_db[real]
    raise KeyError(f"Locus {locus_name!r} not found in frequency database.")


def _ensure_alleles(freqs: Dict[str, float],
                    needed: Iterable[str]) -> Tuple[Dict[str, float], List[str]]:
    """Return a copy of ``freqs`` extended so every allele in ``needed`` is
    present. If some alleles are missing, they are added by borrowing from a
    ``Rest`` bin (when present) or otherwise assigned a small floor (1e-4)
    and the dict is renormalised. Returns the new dict + a warnings list.
    """
    warnings: List[str] = []
    f = dict(freqs)
    needed = [a for a in needed if a]
    missing = [a for a in needed if a not in f]
    if not missing:
        return f, warnings
    rest = f.pop("Rest", None)
    if rest is not None and rest > 0:
        share = rest / len(missing)
        for a in missing:
            f[a] = share
        warnings.append(
            f"Distributed CAP 'Rest' bin ({rest:.4f}) evenly across "
            f"unspecified alleles {missing}."
        )
    else:
        floor = 1e-4
        for a in missing:
            f[a] = floor
        warnings.append(
            f"Allele(s) {missing} not in frequency database; "
            f"assigned floor frequency {floor}."
        )
    s = sum(f.values())
    if s > 0:
        f = {k: v / s for k, v in f.items()}
    return f, warnings


def build_loci(locus_names: List[str],
               freq_db: FreqDB,
               dna: List[TypedAllele],
               *,
               mutation_model: str = "Equal",
               mutation_rate: float = 0.001) -> Tuple[List, pd.DataFrame, List[str]]:
    """Build :class:`FamiliasLocus` objects + a typed datamatrix DataFrame."""
    warnings: List[str] = []
    persons = sorted({d.person for d in dna})
    cols: List[str] = []
    rows = {p: [] for p in persons}
    loci_objs = []
    for locus in locus_names:
        try:
            real, freqs = _resolve_freqs(locus, freq_db)
        except KeyError as e:
            warnings.append(str(e) + " — skipped.")
            continue
        # Collect alleles observed at this locus
        observed: List[str] = []
        per_person: Dict[str, Tuple[str, str]] = {}
        for d in dna:
            if d.locus == locus or canonical_locus_name(d.locus) == canonical_locus_name(locus):
                observed += [d.a1, d.a2]
                per_person[d.person] = (d.a1, d.a2)
        freqs2, w = _ensure_alleles(freqs, observed)
        warnings.extend(f"[{locus}] {m}" for m in w)
        loc = FamiliasLocus(
            frequencies=list(freqs2.values()),
            allelenames=list(freqs2.keys()),
            name=locus,
            MutationModel=mutation_model,
            MutationRate=float(mutation_rate),
        )
        loci_objs.append(loc)
        cols += [f"{locus}.1", f"{locus}.2"]
        for p in persons:
            a1, a2 = per_person.get(p, ("NA", "NA"))
            rows[p].extend([a1, a2])
    if not loci_objs:
        raise ValueError("No usable loci. Check locus names against the frequency database.")
    dm = pd.DataFrame([rows[p] for p in persons], index=persons, columns=cols)
    return loci_objs, dm, warnings


# ---------------------------------------------------------------------------
# Mode 1: paternity duo (AF + CH)
# ---------------------------------------------------------------------------
def make_duo_pedigrees() -> Tuple[FamiliasPedigree, FamiliasPedigree]:
    h1 = FamiliasPedigree(
        id=["AF", "CH"], dadid=["NA", "AF"], momid=["NA", "NA"],
        sex=["male", "male"],
    )
    h2 = FamiliasPedigree(
        id=["AF", "CH"], dadid=["NA", "NA"], momid=["NA", "NA"],
        sex=["male", "male"],
    )
    return h1, h2


# Mode 2: trio (AF + MO + CH)
def make_trio_pedigrees() -> Tuple[FamiliasPedigree, FamiliasPedigree]:
    h1 = FamiliasPedigree(
        id=["AF", "MO", "CH"], dadid=["NA", "NA", "AF"],
        momid=["NA", "NA", "MO"], sex=["male", "female", "male"],
    )
    h2 = FamiliasPedigree(
        id=["AF", "MO", "CH"], dadid=["NA", "NA", "NA"],
        momid=["NA", "NA", "MO"], sex=["male", "female", "male"],
    )
    return h1, h2


# ---------------------------------------------------------------------------
# Mode 3: arbitrary
# ---------------------------------------------------------------------------
def make_pedigrees_from_relations(persons: List[Person],
                                  relations: List[Relation],
                                  ) -> Tuple[FamiliasPedigree, FamiliasPedigree]:
    """Build (H1, H2) pedigrees: with vs without the ``test`` relation(s)."""
    if not any(r.flag == "test" for r in relations):
        raise ValueError(
            "Mode 3 requires at least one relation flagged 'test' (defines "
            "what is being tested)."
        )
    sex_of = {p.id: p.sex for p in persons}
    ids = [p.id for p in persons]

    def _build(include_test: bool) -> FamiliasPedigree:
        dad = ["NA"] * len(ids)
        mom = ["NA"] * len(ids)
        idx = {pid: i for i, pid in enumerate(ids)}
        for r in relations:
            if r.flag == "test" and not include_test:
                continue
            if r.parent not in idx:
                raise ValueError(f"Unknown parent id {r.parent!r}.")
            if r.child not in idx:
                raise ValueError(f"Unknown child id {r.child!r}.")
            i = idx[r.child]
            if sex_of[r.parent] == "male":
                if dad[i] != "NA":
                    raise ValueError(f"{r.child!r} already has a father.")
                dad[i] = r.parent
            else:
                if mom[i] != "NA":
                    raise ValueError(f"{r.child!r} already has a mother.")
                mom[i] = r.parent
        return FamiliasPedigree(id=ids, dadid=dad, momid=mom,
                                sex=[sex_of[i] for i in ids])
    return _build(True), _build(False)


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------
def run_posterior(h1: FamiliasPedigree,
                  h2: FamiliasPedigree,
                  loci: List,
                  dm: pd.DataFrame,
                  *, kinship: float = 0.0) -> Result:
    res = FamiliasPosterior([h1, h2], loci=loci, datamatrix=dm,
                            ref=2, kinship=kinship)
    L = res["likelihoods"]
    LRtot = float(res["LR"][0])
    per = res["likelihoodsPerSystem"]   # shape (nloci, npeds)
    LRm = res["LRperMarker"]
    per_locus = []
    for i, name in enumerate(res["locusnames"]):
        per_locus.append(
            (str(name), float(per[i, 0]), float(per[i, 1]), float(LRm[i, 0]))
        )
    log10 = math.log10(LRtot) if LRtot > 0 else float("-inf")
    return Result(
        LR=LRtot,
        log10_LR=log10,
        posterior=(float(res["posterior"][0]), float(res["posterior"][1])),
        likelihoods=(float(L[0]), float(L[1])),
        per_locus=per_locus,
    )


# ---------------------------------------------------------------------------
# Single-locus quick LR (used by the live TUI table)
# ---------------------------------------------------------------------------
def single_locus_lr(
    *,
    locus_name: str,
    frequencies: Dict[str, float],
    genotypes: Dict[str, Tuple[str, str]],
    h1: FamiliasPedigree,
    h2: FamiliasPedigree,
    mutation_model: str = "Equal",
    mutation_rate: float = 0.001,
    kinship: float = 0.0,
) -> Optional[Tuple[float, float, float]]:
    """Return ``(L_h1, L_h2, LR)`` for a single locus, or ``None`` if the
    inputs are insufficient.

    ``frequencies`` must already contain every allele referenced by
    ``genotypes`` (use :func:`_ensure_alleles` first if necessary).
    """
    if not genotypes:
        return None
    s = sum(frequencies.values())
    if s <= 0:
        return None
    # Preserve the user-supplied (or DB) frequencies verbatim. If they don't
    # sum to 1, top up with a synthetic "Rest" allele so that p_a stays
    # exactly p_a (no spurious renormalisation that would otherwise inflate
    # the rare alleles and bias the LR). If they sum to slightly more than
    # 1 due to rounding, only then renormalise.
    f = {a: float(v) for a, v in frequencies.items() if v and v > 0}
    if not f:
        return None
    s = sum(f.values())
    if s > 1.0 + 1e-9:
        f = {a: v / s for a, v in f.items()}
    elif s < 1.0 - 1e-9:
        rest = 1.0 - s
        # Use a unique allele name unlikely to collide with real STR alleles.
        f["__rest__"] = rest
    fnorm = f
    loc = FamiliasLocus(
        frequencies=list(fnorm.values()),
        allelenames=list(fnorm.keys()),
        name=locus_name,
        MutationModel=mutation_model,
        MutationRate=float(mutation_rate),
    )
    persons = list(genotypes.keys())
    rows = []
    for p in persons:
        a, b = genotypes[p]
        rows.append([a, b])
    dm = pd.DataFrame(rows, index=persons,
                      columns=[f"{locus_name}.1", f"{locus_name}.2"])
    res = FamiliasPosterior([h1, h2], loci=[loc], datamatrix=dm, ref=2,
                            kinship=float(kinship))
    L1 = float(res["likelihoods"][0])
    L2 = float(res["likelihoods"][1])
    if L2 <= 0:
        return (L1, L2, float("nan"))
    return (L1, L2, L1 / L2)


# ---------------------------------------------------------------------------
# Mendelian mismatch detection
# ---------------------------------------------------------------------------
def is_mismatch_one_parent(
    parent: Optional[Tuple[str, str]],
    child: Optional[Tuple[str, str]],
) -> bool:
    """True iff parent and child genotypes share no allele."""
    if parent is None or child is None:
        return False
    return not (set(parent) & set(child))


def is_mismatch_two_parents(
    known: Optional[Tuple[str, str]],
    alleged: Optional[Tuple[str, str]],
    child: Optional[Tuple[str, str]],
) -> bool:
    """True iff no allele assignment makes the trio Mendelian-consistent."""
    if known is None or alleged is None or child is None:
        return False
    ca, cb = child
    for k_alle in known:
        for a_alle in alleged:
            if {k_alle, a_alle} == {ca, cb} or (
                k_alle == ca and a_alle == cb
            ) or (k_alle == cb and a_alle == ca):
                return False
    return True


# Power-of-exclusion (one-parent, mother unknown) — Brenner 1992 form.
def power_of_exclusion_one_parent(child: Tuple[str, str],
                                  freqs: Dict[str, float]) -> Optional[float]:
    """PE = probability a random man would be excluded as the father.

    For child genotype (a,b): a man is *not* excluded iff he transmits
    one of {a,b}. Using random-transmission probability ½(p_a + p_b) for
    each haplotype and HWE for the man's genotype yields::

        P(not excluded) = (p_a + p_b)·(2 - p_a - p_b) - p_a·p_b
        PE              = 1 − P(not excluded)

    Returns ``None`` if a frequency is missing.
    """
    a, b = child
    pa = freqs.get(a)
    pb = freqs.get(b)
    if pa is None or pb is None:
        return None
    s = pa + pb
    p_not = s * (2 - s) - pa * pb
    return max(0.0, min(1.0, 1.0 - p_not))
