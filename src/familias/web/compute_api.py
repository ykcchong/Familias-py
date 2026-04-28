"""Adapter from web request models to :mod:`familias.tui.compute`."""
from __future__ import annotations
from typing import Dict, List, Tuple

import pandas as pd

from .. import FamiliasLocus
from ..tui import compute as tc
from .models import (
    ArbitraryCase,
    CaseBase,
    ComputeResponse,
    GenotypeModel,
    LocusInputModel,
    OneParentCase,
    PerLocusResult,
    PersonModel,
    RelationModel,
    SingleLocusRequest,
    SingleLocusResponse,
    TwoParentCase,
)


# Floor frequency assigned to observed alleles missing from the user-supplied
# frequency dict (mirrors the TUI behaviour in ``compute._ensure_alleles``).
_FLOOR = 1e-4
# Synthetic allele used to make frequencies sum to exactly 1 without
# renormalising user-supplied values (preserves rounded DB inputs verbatim,
# mirroring the trick used in ``compute.single_locus_lr``).
_REST_KEY = "__rest__"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_freq_db(loci: List[LocusInputModel]) -> tc.FreqDB:
    return {l.name: {a: float(v) for a, v in l.frequencies.items()
                     if v and v > 0}
            for l in loci}


def _build_dna(loci: List[LocusInputModel]) -> List[tc.TypedAllele]:
    out: List[tc.TypedAllele] = []
    for l in loci:
        for g in l.genotypes:
            if not g.a1 or not g.a2:
                continue
            out.append(tc.TypedAllele(g.person, l.name, g.a1, g.a2))
    return out


def _posterior_with_prior(L1: float, L2: float, prior_h1: float
                          ) -> Tuple[float, float]:
    """Return (P(H1|D), P(H2|D)) for the supplied prior P(H1)."""
    p1 = float(prior_h1)
    p2 = 1.0 - p1
    num1 = L1 * p1
    num2 = L2 * p2
    s = num1 + num2
    if s <= 0:
        return (p1, p2)
    return (num1 / s, num2 / s)


def _to_response(result: tc.Result, prior_h1: float) -> ComputeResponse:
    L1, L2 = result.likelihoods
    post = _posterior_with_prior(L1, L2, prior_h1)
    return ComputeResponse(
        LR=result.LR,
        log10_LR=result.log10_LR,
        likelihoods=result.likelihoods,
        posterior=result.posterior,
        posterior_h1=post[0],
        per_locus=[PerLocusResult(locus=n, L1=l1, L2=l2, LR=lr)
                   for (n, l1, l2, lr) in result.per_locus],
        warnings=list(result.warnings),
    )


def _build_case_loci(case: CaseBase):
    """Build ``(loci_objs, datamatrix, locus_names, warnings)``.

    Unlike :func:`familias.tui.compute.build_loci`, this preserves the
    user's frequency values verbatim and tops the per-locus sum to exactly
    1 with a synthetic ``__rest__`` allele when the supplied values fall
    short. This matches what the TUI's live single-locus path does and
    avoids inflating rare alleles by renormalising rounded DB inputs.
    """
    warnings: List[str] = []
    locus_names = [l.name for l in case.loci]

    # Collect every distinct person across loci (preserve first-seen order).
    persons: List[str] = []
    seen: set[str] = set()
    for l in case.loci:
        for g in l.genotypes:
            if g.a1 and g.a2 and g.person not in seen:
                seen.add(g.person)
                persons.append(g.person)
    if not persons:
        raise ValueError("No genotypes provided.")

    cols: List[str] = []
    rows: Dict[str, List[str]] = {p: [] for p in persons}
    loci_objs = []
    used_locus_names: List[str] = []

    for l in case.loci:
        per_person: Dict[str, Tuple[str, str]] = {
            g.person: (g.a1, g.a2) for g in l.genotypes if g.a1 and g.a2
        }
        if not per_person:
            continue  # locus has no observations; skip silently.

        observed: List[str] = []
        for a1, a2 in per_person.values():
            for a in (a1, a2):
                if a not in observed:
                    observed.append(a)

        # Start with positive user frequencies (verbatim).
        f: Dict[str, float] = {a: float(v) for a, v in l.frequencies.items()
                               if v and v > 0}
        # Floor any observed allele missing from the user dict.
        missing = [a for a in observed if a not in f]
        if missing:
            for a in missing:
                f[a] = _FLOOR
            warnings.append(
                f"[{l.name}] allele(s) {missing} missing from frequencies; "
                f"assigned floor {_FLOOR}."
            )

        s = sum(f.values())
        if s > 1.0 + 1e-9:
            f = {a: v / s for a, v in f.items()}
        elif s < 1.0 - 1e-9:
            f[_REST_KEY] = 1.0 - s

        loc = FamiliasLocus(
            frequencies=list(f.values()),
            allelenames=list(f.keys()),
            name=l.name,
            MutationModel=case.mutation_model,
            MutationRate=float(case.mutation_rate),
        )
        loci_objs.append(loc)
        used_locus_names.append(l.name)
        cols += [f"{l.name}.1", f"{l.name}.2"]
        for p in persons:
            a1, a2 = per_person.get(p, ("NA", "NA"))
            rows[p].extend([a1, a2])

    if not loci_objs:
        raise ValueError(
            "No usable loci. Each locus must have at least one complete "
            "genotype.")
    dm = pd.DataFrame([rows[p] for p in persons], index=persons, columns=cols)
    return loci_objs, dm, used_locus_names, warnings


# ---------------------------------------------------------------------------
# Mode entry points
# ---------------------------------------------------------------------------
def compute_one_parent(case: OneParentCase) -> ComputeResponse:
    h1, h2 = tc.make_duo_pedigrees()
    loci_objs, dm, _, warns = _build_case_loci(case)
    res = tc.run_posterior(h1, h2, loci_objs, dm, kinship=case.kinship)
    res.warnings.extend(warns)
    return _to_response(res, case.prior)


def compute_two_parent(case: TwoParentCase) -> ComputeResponse:
    h1, h2 = tc.make_trio_pedigrees()
    loci_objs, dm, _, warns = _build_case_loci(case)
    res = tc.run_posterior(h1, h2, loci_objs, dm, kinship=case.kinship)
    res.warnings.extend(warns)
    return _to_response(res, case.prior)


def compute_arbitrary(case: ArbitraryCase) -> ComputeResponse:
    persons = [tc.Person(p.id, p.sex) for p in case.persons]
    relations = [tc.Relation(r.parent, r.child, r.flag) for r in case.relations]
    h1, h2 = tc.make_pedigrees_from_relations(persons, relations)
    loci_objs, dm, _, warns = _build_case_loci(case)
    res = tc.run_posterior(h1, h2, loci_objs, dm, kinship=case.kinship)
    res.warnings.extend(warns)
    return _to_response(res, case.prior)


# ---------------------------------------------------------------------------
# Single-locus quick LR (debounced live updates)
# ---------------------------------------------------------------------------
def compute_single_locus(req: SingleLocusRequest) -> SingleLocusResponse:
    if not req.genotypes:
        return SingleLocusResponse(LR=None, reason="No genotypes.")
    if not req.frequencies:
        return SingleLocusResponse(LR=None, reason="No frequencies.")

    if req.mode == "one-parent":
        h1, h2 = tc.make_duo_pedigrees()
    elif req.mode == "two-parent":
        h1, h2 = tc.make_trio_pedigrees()
    else:
        if not req.persons or not req.relations:
            return SingleLocusResponse(
                LR=None,
                reason="Arbitrary mode requires persons + relations.",
            )
        persons = [tc.Person(p.id, p.sex) for p in req.persons]
        relations = [tc.Relation(r.parent, r.child, r.flag)
                     for r in req.relations]
        try:
            h1, h2 = tc.make_pedigrees_from_relations(persons, relations)
        except ValueError as e:
            return SingleLocusResponse(LR=None, reason=str(e))

    genos: Dict[str, tuple] = {g.person: (g.a1, g.a2) for g in req.genotypes
                               if g.a1 and g.a2}
    if not genos:
        return SingleLocusResponse(LR=None, reason="No complete genotypes.")

    # Preserve user-supplied (DB) frequencies verbatim. Only floor observed
    # alleles missing from the dict; ``single_locus_lr`` will add a __rest__
    # bin to top the sum to 1 without renormalising rare-allele values.
    freqs: Dict[str, float] = {a: float(v) for a, v in req.frequencies.items()
                               if v and v > 0}
    needed: List[str] = []
    for a1, a2 in genos.values():
        for a in (a1, a2):
            if a not in needed:
                needed.append(a)
    for a in needed:
        if a not in freqs:
            freqs[a] = _FLOOR

    out = tc.single_locus_lr(
        locus_name=req.locus,
        frequencies=freqs,
        genotypes=genos,
        h1=h1,
        h2=h2,
        mutation_model=req.mutation_model,
        mutation_rate=req.mutation_rate,
        kinship=req.kinship,
    )
    if out is None:
        return SingleLocusResponse(LR=None, reason="Insufficient inputs.")
    L1, L2, LR = out
    return SingleLocusResponse(LR=LR, L1=L1, L2=L2)
