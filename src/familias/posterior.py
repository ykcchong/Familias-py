"""Port of ``R/FamiliasPosterior.R``."""
from __future__ import annotations
from typing import List, Optional, Sequence, Union
import numpy as np
import pandas as pd  # only used for nicer return labels (optional)

from ._engine.interface import FamInterface
from ._register import normalise_pedigrees, register
from .locus import _FamiliasLocusObj


def _is_na(x) -> bool:
    if x is None:
        return True
    try:
        return bool(np.isnan(x))
    except (TypeError, ValueError):
        return isinstance(x, str) and x in ("NA", "na", "<NA>")


def FamiliasPosterior(
    pedigrees,
    loci,
    datamatrix,                    # 2D array-like with row labels (person ids)
    prior: Optional[Sequence[float]] = None,
    ref: int = 1,
    kinship: float = 0.0,
    simplifyMutations: bool = False,
    *,
    person_ids: Optional[Sequence[str]] = None,
):
    """Compute posterior probabilities, LRs and likelihoods.

    ``datamatrix`` should be a 2-D array of shape ``(n_persons, 2 * n_loci)``.
    Row labels are the person ids (use the optional ``person_ids`` parameter
    if you pass a plain numpy array instead of a pandas DataFrame).
    """
    pedigrees = normalise_pedigrees(pedigrees)
    npeds = len(pedigrees)
    if ref < 1 or ref > npeds:
        raise ValueError("Impossible reference pedigree index.")
    if prior is None:
        prior = np.full(npeds, 1.0 / npeds)
    prior = np.asarray(prior, dtype=float)
    if prior.size != npeds:
        raise ValueError("Length of prior must match number of pedigrees.")
    if (prior < 0).any() or round(float(prior.sum()), 6) != 1:
        raise ValueError("Prior must be non-negative and sum to 1.")

    if isinstance(loci, _FamiliasLocusObj):
        loci = [loci]
    loci = list(loci)
    nloci = len(loci)
    if nloci < 1:
        raise ValueError("At least one locus required.")
    if len({lo["locusname"] for lo in loci}) != nloci:
        raise ValueError("Locus names must be unique.")

    # Materialise datamatrix and row labels
    if isinstance(datamatrix, pd.DataFrame):
        row_labels = list(datamatrix.index.astype(str))
        dm = datamatrix.to_numpy(dtype=object)
    else:
        dm = np.asarray(datamatrix, dtype=object)
        if person_ids is None:
            raise ValueError("person_ids required when datamatrix is not a DataFrame.")
        row_labels = [str(s) for s in person_ids]
    if dm.ndim != 2 or dm.shape[1] != 2 * nloci:
        raise ValueError("datamatrix must be 2-D with 2 * n_loci columns.")
    if dm.shape[0] != len(row_labels):
        raise ValueError("Number of row labels must match datamatrix rows.")

    # Optionally extend with blank rows for untyped pedigree members (kinship>0)
    if kinship > 0:
        all_ped_ids = []
        seen = set()
        for p in pedigrees:
            for pid in p.id:
                if pid not in seen:
                    seen.add(pid)
                    all_ped_ids.append(pid)
        to_add = [pid for pid in all_ped_ids if pid not in row_labels]
        if to_add:
            blanks = np.full((len(to_add), dm.shape[1]), np.nan, dtype=object)
            dm = np.vstack([dm, blanks])
            row_labels = row_labels + to_add

    # Validate that every typed person occurs in the (intersection-of-)pedigrees,
    # and find the common-person ordering used by the engine.
    npers = len(row_labels)
    # The engine treats the persons in row order as person 0..npers-1
    # Each person must appear in *some* pedigree (matching the R semantics).
    for j, name in enumerate(row_labels):
        in_some = any(name in p.id for p in pedigrees)
        typed = not all(_is_na(dm[j, k]) for k in range(dm.shape[1]))
        if typed and not in_some:
            raise ValueError(
                f"Person {name!r} is typed but does not occur in any pedigree."
            )

    # Build engine
    fi = FamInterface()
    register(pedigrees, fi, row_labels, require_consistent_sex=True)

    # Register loci
    for lo in loci:
        n_alleles = len(lo["alleles"])
        simplify = bool(lo.get("simpleMutationMatrices", False) or simplifyMutations)
        has_silent = bool(lo.get("hasSilentAllele", False))
        fi.add_allele_system(
            n_alleles=n_alleles,
            mutation_matrix_female=np.asarray(lo["femaleMutationMatrix"], dtype=float),
            mutation_matrix_male=np.asarray(lo["maleMutationMatrix"], dtype=float),
            simplify_mutation_matrix=simplify,
            frequencies=np.fromiter(lo["alleles"].values(), dtype=float),
            has_silent_allele=has_silent,
        )

    # Register DNA observations
    for i, lo in enumerate(loci):
        names_list = list(lo["alleles"].keys())
        name_to_idx = {n: k for k, n in enumerate(names_list)}
        for j in range(npers):
            a1 = dm[j, 2 * i]
            a2 = dm[j, 2 * i + 1]
            if _is_na(a1) and _is_na(a2):
                continue
            if _is_na(a1):
                a1 = a2
            if _is_na(a2):
                a2 = a1
            try:
                m1 = name_to_idx[str(a1)]
            except KeyError:
                raise ValueError(f"Allele {a1!r} not in locus {lo['locusname']!r}")
            try:
                m2 = name_to_idx[str(a2)]
            except KeyError:
                raise ValueError(f"Allele {a2!r} not in locus {lo['locusname']!r}")
            fi.add_dna_observation(j, i, m1, m2)

    if kinship < 0:
        raise ValueError("Kinship cannot be negative.")
    redundant, _probs, likelihoods = fi.get_probabilities(
        generations_parameter=1.0,
        max_generations=-1,
        inbreeding_parameter=1.0,
        promiscuity_parameter=1.0,
        using_dna_observations=True,
        kinship=float(kinship),
    )
    if any(redundant):
        dups = [k + 1 for k, r in enumerate(redundant) if r]
        raise ValueError(f"Some pedigrees are duplicates: {dups}")

    # likelihoods is shape (npeds, nloci) — match R's (nloci, npeds) layout
    likelihoods_per_system = np.asarray(likelihoods, dtype=float).T  # (nloci, npeds)
    likelihoods_total = likelihoods_per_system.prod(axis=0)  # (npeds,)
    posterior = prior * likelihoods_total
    s = posterior.sum()
    if s > 0:
        posterior = posterior / s
    ref0 = ref - 1
    LR = likelihoods_total / likelihoods_total[ref0]
    LRperMarker = likelihoods_per_system / likelihoods_per_system[:, ref0:ref0 + 1]
    locus_names = [lo["locusname"] for lo in loci]

    return {
        "posterior": posterior,
        "prior": prior,
        "LR": LR,
        "LRperMarker": LRperMarker,
        "likelihoods": likelihoods_total,
        "likelihoodsPerSystem": likelihoods_per_system,
        "locusnames": locus_names,
    }
