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

Optional plotting:

```bash
pip install -e ".[plot]"
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

* the bundled `NorwegianFrequencies`,
* a JSON file (`{locus: {allele: freq}}`),
* an R `.rda` file (requires the optional `pyreadr` extra).

Install + run:

```bash
pip install -e ".[tui]"          # add textual + pandas
familias-tui                      # or: python -m familias.tui
```
