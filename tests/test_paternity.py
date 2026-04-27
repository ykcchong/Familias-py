"""Basic paternity-test smoke tests.

For a single-locus paternity case with allele frequencies (A=0.2, B=0.3,
C=0.5) and AF=A/B, MO=A/C, CH=A/A, the textbook LR for "AF is the father"
vs "AF unrelated" is 1 / (2·p_A) = 2.5.
"""
import numpy as np
import pandas as pd
from familias import (
    FamiliasLocus, FamiliasPedigree, FamiliasPosterior, FamiliasPrior,
    NorwegianFrequencies,
)


def _make_peds():
    p1 = FamiliasPedigree(
        id=["AF", "MO", "CH"],
        dadid=["NA", "NA", "AF"],
        momid=["NA", "NA", "MO"],
        sex=["male", "female", "male"],
    )
    p2 = FamiliasPedigree(
        id=["AF", "MO", "CH"],
        dadid=["NA", "NA", "NA"],
        momid=["NA", "NA", "MO"],
        sex=["male", "female", "male"],
    )
    return p1, p2


def test_paternity_textbook_LR_no_mutation():
    locus = FamiliasLocus(
        frequencies=[0.2, 0.3, 0.5],
        allelenames=["A", "B", "C"],
        name="L1",
        MutationRate=0.0,
    )
    p1, p2 = _make_peds()
    dm = pd.DataFrame(
        [["A", "B"], ["A", "C"], ["A", "A"]],
        index=["AF", "MO", "CH"],
        columns=["L1.1", "L1.2"],
    )
    res = FamiliasPosterior([p1, p2], [locus], dm, ref=2)
    assert np.isclose(res["LR"][0], 1 / (2 * 0.2))
    assert np.isclose(res["LR"][1], 1.0)
    # Likelihoods (closed form):
    # L1 = 2·pA·pB · 2·pA·pC · pA = 0.006
    # L2 = 2·pA·pB · 2·pA·pC · pA^2 = 0.0024
    assert np.isclose(res["likelihoods"][0], 0.006)
    assert np.isclose(res["likelihoods"][1], 0.0024)


def test_paternity_with_mutation_runs():
    locus = FamiliasLocus(
        frequencies=[0.2, 0.3, 0.5],
        allelenames=["A", "B", "C"],
        name="L1",
        MutationModel="Equal",
        MutationRate=0.005,
    )
    p1, p2 = _make_peds()
    dm = pd.DataFrame(
        [["A", "B"], ["A", "C"], ["A", "A"]],
        index=["AF", "MO", "CH"],
        columns=["L1.1", "L1.2"],
    )
    res = FamiliasPosterior([p1, p2], [locus], dm, ref=2)
    assert res["LR"][0] > 1.0
    assert np.isclose(res["LR"][1], 1.0)
    assert np.isclose(res["posterior"].sum(), 1.0)


def test_paternity_norwegian_locus():
    freqs = NorwegianFrequencies["D3S1358"]
    names = list(freqs.keys())
    vals = list(freqs.values())
    locus = FamiliasLocus(frequencies=vals, allelenames=names, name="D3S1358",
                          MutationRate=0.001, MutationModel="Equal")
    p1, p2 = _make_peds()
    a = names[0]
    b = names[1]
    c = names[2]
    dm = pd.DataFrame(
        [[a, b], [a, c], [a, a]],
        index=["AF", "MO", "CH"],
        columns=["D3S1358.1", "D3S1358.2"],
    )
    res = FamiliasPosterior([p1, p2], [locus], dm, ref=2)
    # AF carries an A; LR for paternity should exceed 1.
    assert res["LR"][0] > 1.0


def test_prior_uniform():
    p1, p2 = _make_peds()
    prior = FamiliasPrior([p1, p2])
    assert prior.shape == (2,)
    assert np.isclose(prior.sum(), 1.0)
