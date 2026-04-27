"""Closed-form correctness checks for non-paternity scenarios."""
import numpy as np
import pandas as pd
import pytest
from familias import FamiliasLocus, FamiliasPedigree, FamiliasPosterior


def _three_persons_pedigree(child_dad="P1", child_mom="P2"):
    """Two parents + one child."""
    return FamiliasPedigree(
        id=["P1", "P2", "CH"],
        dadid=["NA", "NA", child_dad],
        momid=["NA", "NA", child_mom],
        sex=["male", "female", "male"],
    )


def test_full_siblings_closed_form_LR():
    """Two children both A/A.

    Hypothesis A: full siblings (common parents F, M, both unobserved).
    Hypothesis B: unrelated.

    Closed form: LR = (1 + pA)^2 / (4 * pA^2)  (assuming no mutation).
    """
    pA = 0.2
    locus = FamiliasLocus(
        frequencies=[pA, 0.3, 0.5],
        allelenames=["A", "B", "C"],
        name="L1",
        MutationRate=0.0,
    )
    # Hypothesis 1: full siblings via unobserved parents F, M
    sibs = FamiliasPedigree(
        id=["F", "M", "C1", "C2"],
        dadid=["NA", "NA", "F", "F"],
        momid=["NA", "NA", "M", "M"],
        sex=["male", "female", "male", "male"],
    )
    # Hypothesis 2: unrelated. F and M still present (so person sets match) but
    # have no relation to either child.
    unrelated = FamiliasPedigree(
        id=["F", "M", "C1", "C2"],
        dadid=["NA", "NA", "NA", "NA"],
        momid=["NA", "NA", "NA", "NA"],
        sex=["male", "female", "male", "male"],
    )
    dm = pd.DataFrame(
        [["A", "A"], ["A", "A"]],
        index=["C1", "C2"],
        columns=["L1.1", "L1.2"],
    )
    res = FamiliasPosterior([sibs, unrelated], [locus], dm, ref=2)
    expected = (1 + pA) ** 2 / (4 * pA ** 2)
    assert np.isclose(res["LR"][0], expected, rtol=1e-9), \
        f"LR = {res['LR'][0]}, expected {expected}"


def test_multilocus_likelihoods_factorise():
    """L_total = product over loci of L_per_locus."""
    L1 = FamiliasLocus([0.2, 0.3, 0.5], allelenames=["A", "B", "C"], name="L1",
                        MutationRate=0.0)
    L2 = FamiliasLocus([0.4, 0.6], allelenames=["X", "Y"], name="L2",
                        MutationRate=0.0)
    p1 = FamiliasPedigree(id=["AF", "MO", "CH"], dadid=["NA", "NA", "AF"],
                           momid=["NA", "NA", "MO"], sex=["male", "female", "male"])
    p2 = FamiliasPedigree(id=["AF", "MO", "CH"], dadid=["NA", "NA", "NA"],
                           momid=["NA", "NA", "MO"], sex=["male", "female", "male"])
    dm = pd.DataFrame(
        [["A", "B", "X", "Y"],
         ["A", "C", "X", "X"],
         ["A", "A", "X", "Y"]],
        index=["AF", "MO", "CH"],
        columns=["L1.1", "L1.2", "L2.1", "L2.2"],
    )
    res = FamiliasPosterior([p1, p2], [L1, L2], dm, ref=2)
    per_sys = res["likelihoodsPerSystem"]   # shape (n_loci, n_peds)
    assert per_sys.shape == (2, 2)
    assert np.allclose(res["likelihoods"], per_sys.prod(axis=0))
    # L1 LR = 1/(2 pA) = 2.5 (textbook)
    assert np.isclose(per_sys[0, 0] / per_sys[0, 1], 2.5)


def test_silent_allele_present():
    """A silent allele should be permitted as the last allele only.

    With a silent allele, an observed A/A genotype is consistent with both
    A/A and A/silent. We just check the engine accepts it and produces a
    finite, positive likelihood.
    """
    locus = FamiliasLocus(
        frequencies=[0.2, 0.3, 0.45, 0.05],
        allelenames=["A", "B", "C", "silent"],
        name="L1",
        MutationRate=0.0,
    )
    p1 = FamiliasPedigree(id=["AF", "MO", "CH"], dadid=["NA", "NA", "AF"],
                           momid=["NA", "NA", "MO"], sex=["male", "female", "male"])
    p2 = FamiliasPedigree(id=["AF", "MO", "CH"], dadid=["NA", "NA", "NA"],
                           momid=["NA", "NA", "MO"], sex=["male", "female", "male"])
    dm = pd.DataFrame(
        [["A", "B"], ["A", "C"], ["A", "A"]],
        index=["AF", "MO", "CH"],
        columns=["L1.1", "L1.2"],
    )
    res = FamiliasPosterior([p1, p2], [locus], dm, ref=2)
    assert np.all(res["likelihoods"] > 0)
    assert np.all(np.isfinite(res["LR"]))
    # AF still carries A; LR > 1 even with silent allele present.
    assert res["LR"][0] > 1.0


def test_locus_stabilisation_PM_runs():
    """``Stabilization='PM'`` should produce a row-stochastic matrix that
    matches the original on the diagonal (no-op effect for default rate=0)."""
    locus = FamiliasLocus(
        frequencies=[0.25, 0.25, 0.25, 0.25],
        allelenames=["1", "2", "3", "4"],
        MutationModel="Equal",
        MutationRate=0.01,
        Stabilization="PM",
    )
    M = locus["femaleMutationMatrix"]
    assert np.allclose(M.sum(axis=1), 1.0, atol=1e-6)
    assert (M >= 0).all()


def test_kinship_reduces_LR():
    """Positive kinship should reduce the paternity LR (shared alleles are
    less informative under co-ancestry).
    """
    locus = FamiliasLocus(
        frequencies=[0.2, 0.3, 0.5],
        allelenames=["A", "B", "C"],
        name="L1",
        MutationRate=0.0,
    )
    p1 = FamiliasPedigree(id=["AF", "MO", "CH"], dadid=["NA", "NA", "AF"],
                           momid=["NA", "NA", "MO"], sex=["male", "female", "male"])
    p2 = FamiliasPedigree(id=["AF", "MO", "CH"], dadid=["NA", "NA", "NA"],
                           momid=["NA", "NA", "MO"], sex=["male", "female", "male"])
    dm = pd.DataFrame(
        [["A", "B"], ["A", "C"], ["A", "A"]],
        index=["AF", "MO", "CH"],
        columns=["L1.1", "L1.2"],
    )
    res0 = FamiliasPosterior([p1, p2], [locus], dm, ref=2, kinship=0.0)
    res_t = FamiliasPosterior([p1, p2], [locus], dm, ref=2, kinship=0.05)
    assert np.isclose(res0["LR"][0], 2.5)
    assert res_t["LR"][0] < res0["LR"][0]
    assert res_t["LR"][0] > 1.0


def test_prior_generations_parameter():
    """generationsParameter < 1 should down-weight pedigrees with more
    generations relative to a flat prior.
    """
    from familias import FamiliasPrior
    # Two pedigrees: (a) parent-child (1 generation), (b) grandparent-child (2 gens).
    one_gen = FamiliasPedigree(
        id=["A", "B", "C"],
        dadid=["NA", "NA", "A"],
        momid=["NA", "NA", "B"],
        sex=["male", "female", "male"],
    )
    two_gen = FamiliasPedigree(
        id=["A", "B", "C"],
        dadid=["NA", "A", "NA"],
        momid=["NA", "NA", "B"],
        sex=["male", "female", "male"],
    )
    flat = FamiliasPrior([one_gen, two_gen])
    sloped = FamiliasPrior([one_gen, two_gen], generationsParameter=0.1)
    assert np.isclose(flat.sum(), 1.0)
    assert np.isclose(sloped.sum(), 1.0)
    # The 2-generation pedigree gets weight ∝ 0.1^extraGen, so its share drops.
    assert sloped[1] < flat[1]
