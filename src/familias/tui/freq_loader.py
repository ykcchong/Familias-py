"""Allele-frequency loaders.

Three formats are recognised, dispatched on file extension:

* ``.csv`` — FSIgen format used by the ``forensicpopdata`` package
  (one column per locus, first column = allele, last row = ``N``).
* ``.txt`` — CAP-style: tab-separated blocks of ``allele\\tfreq`` lines
  with locus headers and an optional ``Rest`` row for the unspecified-allele
  bin.
* ``.json`` — ``{locus: {allele: freq}}``.

In addition, ``builtin_databases()`` exposes the bundled
``NorwegianFrequencies`` plus any ``.json``, ``.csv``, ``.txt``,
``.tsv`` or ``.rda`` files found in the ``data/`` directory shipped
inside the ``familias`` package.
"""
from __future__ import annotations
from importlib import resources
from pathlib import Path
from typing import Dict
import csv
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
# FSIgen-style CSV (used by forensicpopdata)
# ---------------------------------------------------------------------------
def load_fsigen_csv(path: str | Path,
                    *, remove_zeroes: bool = True,
                    normalise: bool = True) -> FreqDB:
    rows = []
    with open(path, newline="") as fh:
        rdr = csv.reader(fh)
        for row in rdr:
            if any(c.strip() for c in row):
                rows.append(row)
    if len(rows) < 3:
        raise ValueError(f"{path}: not enough rows for FSIgen format.")
    header = rows[0]
    body = rows[1:-1]                       # last row is N
    out: FreqDB = {}
    for j in range(1, len(header)):
        locus = header[j].strip()
        if not locus:
            continue
        d: Dict[str, float] = {}
        for r in body:
            if j >= len(r):
                continue
            allele = r[0].strip()
            cell = r[j].strip()
            if not allele or not cell:
                continue
            try:
                v = float(cell)
            except ValueError:
                continue
            if remove_zeroes and v == 0.0:
                continue
            d[allele] = v
        if not d:
            continue
        if normalise:
            s = sum(d.values())
            if s > 0:
                d = {k: v / s for k, v in d.items()}
        out[locus] = d
    return out


# ---------------------------------------------------------------------------
# CAP "Dry challenge" tab-separated TXT
# ---------------------------------------------------------------------------
def load_cap_txt(path: str | Path,
                 *, normalise: bool = True) -> FreqDB:
    """Parse the CAP allele-frequency database TXT.

    Layout::

        CSF1PO\\t
        11\\t0.2797
        12\\t0.375
        Rest\\t0.3453
        \\t              <- blank separator
        D2S441\\t
        ...

    The ``Rest`` bin is preserved as a pseudo-allele named ``Rest`` so that
    case alleles which fall outside the explicitly-listed ones can be
    mapped onto it.
    """
    out: FreqDB = {}
    cur: str | None = None
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            parts = [p.strip() for p in line.split("\t")]
            if not any(parts):                       # blank separator
                cur = None
                continue
            head, *rest = parts
            tail = rest[0] if rest else ""
            if not tail:                              # locus header
                cur = head
                out.setdefault(cur, {})
                continue
            if cur is None:
                continue
            try:
                v = float(tail)
            except ValueError:
                continue
            if v > 0:
                out[cur][head] = v
    if normalise:
        for k, d in list(out.items()):
            s = sum(d.values())
            if s > 0:
                out[k] = {a: v / s for a, v in d.items()}
    # drop empty loci
    return {k: v for k, v in out.items() if v}


# ---------------------------------------------------------------------------
# Auto-dispatch
# ---------------------------------------------------------------------------
def load_freq_file(path: str | Path) -> FreqDB:
    p = Path(path)
    suf = p.suffix.lower()
    if suf == ".json":
        return load_json(p)
    if suf == ".csv":
        return load_fsigen_csv(p)
    if suf in (".txt", ".tsv"):
        return load_cap_txt(p)
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
    directory shipped inside the ``familias`` package for any ``.json``,
    ``.csv``, ``.txt``, ``.tsv`` or ``.rda`` files.
    """
    from .. import NorwegianFrequencies
    out: Dict[str, FreqDB] = {"NorwegianFrequencies (built-in)": dict(NorwegianFrequencies)}

    data_dir = resources.files("familias") / "data"
    if data_dir.is_dir():
        for f in sorted(data_dir.iterdir()):
            if not f.is_file():
                continue
            suf = f.suffix.lower()
            if suf not in (".json", ".csv", ".txt", ".tsv", ".rda"):
                continue
            # Skip the built-in Norwegian frequencies to avoid duplication.
            if f.stem == "norwegian_frequencies":
                continue
            try:
                out[f.stem] = load_freq_file(f)
            except Exception:
                pass
    return out
