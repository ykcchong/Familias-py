"""Embedded data: NorwegianFrequencies.

A dict mapping locus name -> dict of {allele_name: frequency},
extracted from the original ``data/NorwegianFrequencies.rda`` shipped with
the R package (35 STR markers; allele names are repeat-count strings).
"""
from __future__ import annotations
import json
from importlib import resources

_DATA_FILE = "norwegian_frequencies.json"


def _load() -> dict:
    with resources.files(__package__).joinpath("data", _DATA_FILE).open("r") as f:
        raw = json.load(f)
    # Preserve allele insertion order (Python 3.7+ dicts do this) so user sees
    # the same allele ordering as in R.
    return {locus: dict(alleles.items()) for locus, alleles in raw.items()}


NorwegianFrequencies: dict[str, dict[str, float]] = _load()

__all__ = ["NorwegianFrequencies"]
