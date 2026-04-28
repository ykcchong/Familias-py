# familias (Python port)

A pure-Python re-implementation of the R package
[**Familias**](https://cran.r-project.org/package=Familias) v2.6.4 by
P. Mostad, T. Egeland, F. Marsico and I. Simonsson.

`Familias` computes probabilities for pedigrees given DNA marker data, using
the Elston-Stewart peeling algorithm with cut-set based dynamic programming.
It supports several mutation models (Equal, Proportional, Stepwise, Custom),
silent alleles, theta-correction (kinship), and prior weighting based on
generations / inbreeding / promiscuity.

The original C++ engine has been ported directly to Python. NumPy / SciPy are
used where convenient (mutation matrices, optimisation, root-finding) but the
peeling recursion itself is plain Python, mirroring the C++ structure.

## Install

```bash
pip install -e .
```

Optional extras:

```bash
pip install -e ".[plot]"        # matplotlib + networkx pedigree drawing
pip install -e ".[tui]"         # Textual terminal UI (familias-tui)
pip install -e ".[web]"         # FastAPI/uvicorn web UI (familias-web)
pip install -e ".[rda]"         # load R .rda allele-frequency files
pip install -e ".[test]"        # pytest + httpx for the test suite
```

## Quick start

```python
import pandas as pd
from familias import (
    FamiliasLocus,
    FamiliasPedigree,
    FamiliasPosterior,
    NorwegianFrequencies,
)

# Hypothesis 1: AF is the alleged father of CH; M is CH's mother.
ped1 = FamiliasPedigree(
    id=["AF", "M", "CH"],
    dadid=["NA", "NA", "AF"],
    momid=["NA", "NA", "M"],
    sex=["male", "female", "male"],
)
# Hypothesis 2: AF is unrelated to CH.
ped2 = FamiliasPedigree(
    id=["AF", "M", "CH"],
    dadid=["NA", "NA", "NA"],
    momid=["NA", "NA", "M"],
    sex=["male", "female", "male"],
)

# Build a locus from the bundled Norwegian allele frequencies.
freqs = NorwegianFrequencies["TPOX"]
loc = FamiliasLocus(
    list(freqs.values()),
    allelenames=list(freqs.keys()),
    name="TPOX",
    MutationModel="Equal",
    MutationRate=0.001,
)

# datamatrix: rows = persons, columns = (allele1, allele2) per locus
dm = pd.DataFrame(
    [["8", "9"],     # AF
     ["8", "10"],    # M
     ["8", "10"]],   # CH
    index=["AF", "M", "CH"],
    columns=["TPOX.1", "TPOX.2"],
)
res = FamiliasPosterior([ped1, ped2], loci=[loc], datamatrix=dm, ref=2)
print("LR(ped1 vs ped2):", res["LR"][0])
print("Posterior:       ", res["posterior"])
```

## Verified correctness

The test suite checks against closed-form formulas, e.g. for the canonical
paternity case with AF=A/B, M=A/C, CH=A/A and no mutation:
``LR = 1 / (2·p_A)``; for two full siblings both A/A:
``LR = (1 + p_A)² / (4·p_A²)``.

## Public API

* `FamiliasLocus(frequencies, name=..., allele_names=..., ...)`
* `FamiliasPedigree(id, dadid, momid, sex)`
* `FamiliasPrior(pedigrees, ...)`
* `FamiliasPosterior(pedigrees, loci, datamatrix, ...)`
* `plot_familias_pedigree(ped)`
* `NorwegianFrequencies` (a dict of locus -> {allele: freq})

### Bundled allele-frequency datasets

The package ships two locus-frequency databases under `familias/data/`:

* `norwegian_frequencies.json` — the legacy R-package dataset, exposed
  programmatically as `familias.NorwegianFrequencies`. Not exposed in the
  TUI / web picker.
* `fsi.json` — Hong Kong Chinese frequencies, exposed in the TUI / web
  picker as **Hong Kong Chinese** (the default selection).

## Low-level engine API

The `familias._engine.FamInterface` class mirrors the C++ `FamInterface`
1:1 (the same underlying engine the original R package binds to). It is
useful for building pedigrees programmatically, enumerating all valid
pedigrees from extras, applying year-of-birth (YOB) / is-child fertility
constraints, or attaching fixed relations:

```python
from familias._engine.interface import FamInterface

fi = FamInterface()
af  = fi.add_person(male=True,  yob=1970)
mo  = fi.add_person(male=False, yob=1972)
ch  = fi.add_person(male=True,  yob=2000, is_child=True)
fi.add_pedigree(0, 0)
fi.add_relation(af, ch, 0)
fi.add_relation(mo, ch, 0)
fi.add_fixed_relation(af, ch)        # propagate to all pedigrees
fi.generate_pedigrees(0, 1)          # enumerate alternates with one extra male
print(fi.get_number_of_pedigrees())
```

YOB sentinels follow the C++/R convention: `0` or `-1` mean "unknown".

## Status

This is a faithful port of the algorithm in v2.6.4. Numerical results are
intended to agree with the R version. Please report any discrepancies.

## License

GPL-3.0-or-later (matching the upstream R package).

## Terminal UI

A [Textual](https://textual.textualize.io/) front-end is bundled. It exposes
three modes:

1. **Parent-child** — alleged parent + child (1 trio member missing)
2. **Trio** — alleged father + mother + child
3. **Arbitrary** — define persons, relations and tag one or more relations as
   `test`. The app builds two pedigrees (with vs. without the test relations)
   and reports the LR / posterior.

Modes 1 and 2 default to the 15-locus panel:

```
D8S1179, D21S11, D7S820, CSF1PO, D3S1358, THO1, D13S317, D16S539,
D2S1338, D19S433, vWA, TPOX, D18S51, D5S818, FGA
```

Allele frequencies can be loaded from:

* the bundled **Hong Kong Chinese** database (shown as
  `Hong Kong Chinese` in the picker; sourced from `src/familias/data/fsi.json`),
* `NorwegianFrequencies` (importable as a Python dict, but not exposed in the
  TUI/web picker by default),
* a JSON file (`{locus: {allele: freq}}`),
* an R `.rda` file (requires the optional `pyreadr` extra).

Install + run:

```bash
pip install -e ".[tui]"          # add textual + pandas
familias-tui                      # or: python -m familias.tui
```

## Web UI

A FastAPI + React/AG Grid web front-end is also bundled. It mirrors the TUI
features and is served by a single `uvicorn` process that hosts both the
REST API (`/api/*`) and the pre-built single-page application from
`familias/web/dist/`. Three workflows are exposed (matching the TUI):

1. **One-parent** — Tested Parent vs. unrelated. Two persons
   (Tested Parent, Child).
2. **Two-parent** — Known Parent + Tested Parent vs. Known Parent only.
   Three persons (Tested Parent, Known Parent, Child).
3. **Arbitrary** — up to 8 persons with custom parent → child relations;
   tag one or more as `test` to define H₁ vs H₂.

Person column headers use the friendly labels above; the underlying wire
identifiers remain `AF`/`MO`/`CH` so the bundled pedigree templates apply
unchanged.

Highlights:

* **Live per-locus LR** — debounced single-locus calls update the LR column
  as you type genotypes / frequencies.
* **Plain-decimal output** — combined and per-locus LR are shown as
  `toFixed(4)` (no scientific notation); posterior P(H₁) is shown as a
  percentage (`xx.xxxx%`).
* **Frequency database picker** — defaults to the bundled
  **Hong Kong Chinese** database. A **Manual Input** entry disables
  database-driven auto-fill (used automatically after loading a JSON case
  file so user-supplied frequencies are preserved verbatim). Choosing a
  real database and clicking **Apply DB** overwrites observed-allele
  frequencies for the loaded data.
* **Verbatim frequencies + Rest top-up** — frequencies are used as-is; if
  the per-locus sum is < 1, the missing mass is placed in a synthetic
  `__rest__` allele rather than renormalising rare alleles.
* **HTTPS** — pass `--ssl-certfile` / `--ssl-keyfile` (and optionally
  `--ssl-keyfile-password`) to enable TLS.

Install + run:

```bash
pip install -e ".[web]"           # add fastapi + uvicorn + pandas
familias-web                       # default: http://127.0.0.1:8765/
# Common options:
familias-web --host 0.0.0.0 --port 8443 \
    --ssl-certfile cert.pem --ssl-keyfile key.pem
familias-web --no-browser --reload     # development
```

The web UI ships its compiled bundle inside the wheel, so the JavaScript
toolchain is **not** required at install time. Rebuilding the bundle
(only needed when changing `frontend/src/`) requires Node 18+:

```bash
cd frontend
npm install
npm run build      # writes ../src/familias/web/dist/
```
