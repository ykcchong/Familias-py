"""Allele-frequency loaders.

Two formats are recognised, dispatched on file extension:

* ``.json`` — ``{locus: {allele: freq}}``.
* ``.rda`` — R data file (requires the optional ``pyreadr`` dependency).

In addition, ``builtin_databases()`` exposes the bundled
``NorwegianFrequencies`` plus any ``.json`` or ``.rda`` files found in the
``data/`` directory shipped inside the ``familias`` package.
"""
from __future__ import annotations
from importlib import resources
from pathlib import Path
from typing import Dict
import json


FreqDB = Dict[str, Dict[str, float]]


# ---------------------------------------------------------------------------
# JSON
# ---------------------------------------------------------------------------
def load_json(path: str | Path) -> FreqDB:
    with open(path, "r") as fh:
        raw = json.load(fh)
    return {str(k): {str(a): float(v) for a, v in d.items()} for k, d in raw.items()}


# ---------------------------------------------------------------------------
# Auto-dispatch
# ---------------------------------------------------------------------------
def load_freq_file(path: str | Path) -> FreqDB:
    p = Path(path)
    suf = p.suffix.lower()
    if suf == ".json":
        return load_json(p)
    if suf == ".rda":
        return load_rda(p)
    raise ValueError(f"Unrecognised frequency-file extension: {suf!r}")


def load_rda(path: str | Path) -> FreqDB:
    """Load an R ``.rda`` data file via :mod:`pyreadr` (optional dep).

    The expected R object is a list-of-populations-of-locus-frequency-vectors
    such as ``FBI2015freqs``; if a single such population is found, it is
    returned. Otherwise the first population is used and a hint is printed.
    """
    try:
        import pyreadr  # type: ignore
    except ImportError as e:                                  # pragma: no cover
        raise RuntimeError(
            "Reading .rda files requires the optional dependency "
            "'pyreadr' (pip install pyreadr)."
        ) from e
    res = pyreadr.read_r(str(path))
    # pyreadr returns OrderedDict[name -> DataFrame]; population databases
    # are nested lists which pyreadr flattens to multiple DataFrames keyed
    # by '<top>$<pop>'. We accept either flat or nested.
    out: FreqDB = {}
    for _name, df in res.items():
        # heuristic: dataframe with columns like ['locus','allele','frequency']
        cols_lower = {c.lower(): c for c in df.columns}
        if {"locus", "allele"}.issubset(cols_lower):
            fcol = cols_lower.get("frequency") or cols_lower.get("freq")
            if fcol is None:
                continue
            for _, row in df.iterrows():
                loc = str(row[cols_lower["locus"]])
                al = str(row[cols_lower["allele"]])
                v = float(row[fcol])
                if v > 0:
                    out.setdefault(loc, {})[al] = v
    if not out:
        raise ValueError(
            f"{path}: could not extract a (locus, allele, frequency) "
            "table from this .rda file."
        )
    # normalise
    for k, d in out.items():
        s = sum(d.values())
        if s > 0:
            out[k] = {a: v / s for a, v in d.items()}
    return out


# ---------------------------------------------------------------------------
# Built-in / repo-discovered databases
# ---------------------------------------------------------------------------
def builtin_databases() -> Dict[str, FreqDB]:
    """Return ``{name: freq_db}`` for each readily-available database.

    Always includes ``NorwegianFrequencies``. Also scans the ``data/``
    directory shipped inside the ``familias`` package for any ``.json``
    or ``.rda`` files.
    """
    from .. import NorwegianFrequencies
    out: Dict[str, FreqDB] = {"NorwegianFrequencies (built-in)": dict(NorwegianFrequencies)}

    data_dir = resources.files("familias") / "data"
    if data_dir.is_dir():
        for f in sorted(data_dir.iterdir()):
            if not f.is_file():
                continue
            suf = f.suffix.lower()
            if suf not in (".json", ".rda"):
                continue
            # Skip the built-in Norwegian frequencies to avoid duplication.
            if f.stem == "norwegian_frequencies":
                continue
            try:
                out[f.stem] = load_freq_file(f)
            except Exception:
                pass
    return out
