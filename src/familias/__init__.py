"""familias - Python port of the R package Familias (v2.6.4).

Forensic-genetics pedigree probability calculations.
"""
from __future__ import annotations

from .locus import FamiliasLocus
from .pedigree import FamiliasPedigree
from .prior import FamiliasPrior
from .posterior import FamiliasPosterior
from ._data import NorwegianFrequencies

__all__ = [
    "FamiliasLocus",
    "FamiliasPedigree",
    "FamiliasPrior",
    "FamiliasPosterior",
    "NorwegianFrequencies",
    "plot_familias_pedigree",
]

__version__ = "2.6.4"


def plot_familias_pedigree(ped, **kwargs):  # pragma: no cover - thin wrapper
    """Convenience wrapper; imports matplotlib lazily."""
    from .plot import plot_familias_pedigree as _impl
    return _impl(ped, **kwargs)
