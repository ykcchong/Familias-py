"""Port of ``R/FamiliasPrior.R``."""
from __future__ import annotations
from typing import List, Optional
import numpy as np
from ._engine.interface import FamInterface
from ._register import normalise_pedigrees, _common_persons, register


def FamiliasPrior(
    pedigrees,
    generationsParameter: float = 1.0,
    inbreedingParameter: float = 1.0,
    partnerParameter: float = 1.0,
    maxGenerations: Optional[int] = None,
) -> np.ndarray:
    """Return the prior distribution over ``pedigrees``."""
    pedigrees = normalise_pedigrees(pedigrees)
    if generationsParameter < 0 or inbreedingParameter < 0 or partnerParameter < 0:
        raise ValueError("Parameters cannot be negative.")
    persons = _common_persons(pedigrees)
    if len(persons) < 2:
        raise ValueError("At least two persons must be common to all pedigrees.")
    fi = FamInterface()
    register(pedigrees, fi, persons)
    redundant, probs, _ = fi.get_probabilities(
        generationsParameter,
        -1 if maxGenerations is None else int(maxGenerations),
        inbreedingParameter,
        partnerParameter,
        using_dna_observations=False,
        kinship=0.0,
    )
    if any(redundant):
        raise ValueError("Some pedigrees are duplicates; remove duplicates.")
    return np.asarray(probs, dtype=float)
