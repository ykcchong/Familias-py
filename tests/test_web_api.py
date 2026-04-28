"""Tests for the Familias web REST API.

Exercises every endpoint via FastAPI's TestClient and asserts numerical
parity with :mod:`familias.tui.compute` so the web layer cannot diverge
from the TUI / engine.
"""
from __future__ import annotations
import math
import pytest

pytest.importorskip("fastapi")
pytest.importorskip("pandas")

from fastapi.testclient import TestClient

from familias.tui import compute as tc
from familias.web.app import create_app


@pytest.fixture(scope="module")
def client() -> TestClient:
    return TestClient(create_app())


# ---------------------------------------------------------------------------
# Meta / frequencies
# ---------------------------------------------------------------------------
def test_health(client: TestClient) -> None:
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.json() == {"ok": True}


def test_databases_excludes_norwegian(client: TestClient) -> None:
    r = client.get("/api/frequencies/databases")
    assert r.status_code == 200
    body = r.json()
    assert isinstance(body["databases"], list) and body["databases"]
    # Norwegian DB is no longer exposed via the picker.
    assert not any("norwegian" in d.lower() for d in body["databases"])
    assert body["default"] in body["databases"]


def test_database_round_trip(client: TestClient) -> None:
    dbs = client.get("/api/frequencies/databases").json()["databases"]
    name = dbs[0]
    r = client.get(f"/api/frequencies/{name}")
    assert r.status_code == 200
    body = r.json()
    assert body["database"] == name
    assert isinstance(body["loci"], dict) and body["loci"]
    # sanity: every freq dict sums to ~1.0.
    for _loc, freqs in body["loci"].items():
        s = sum(freqs.values())
        assert 0.99 <= s <= 1.01


def test_database_loci_endpoint(client: TestClient) -> None:
    dbs = client.get("/api/frequencies/databases").json()["databases"]
    name = dbs[0]
    r = client.get(f"/api/frequencies/{name}/loci")
    assert r.status_code == 200
    loci = r.json()["loci"]
    assert isinstance(loci, list) and len(loci) > 5


def test_default_loci(client: TestClient) -> None:
    r = client.get("/api/defaults/loci")
    assert r.status_code == 200
    loci = r.json()["loci"]
    assert "D8S1179" in loci and len(loci) == 15


# ---------------------------------------------------------------------------
# Single-locus quick LR
# ---------------------------------------------------------------------------
def test_single_locus_two_parent_textbook(client: TestClient) -> None:
    r = client.post("/api/compute/single-locus", json={
        "mode": "two-parent",
        "locus": "L1",
        "frequencies": {"A": 0.2, "B": 0.3, "C": 0.5},
        "genotypes": [
            {"person": "AF", "a1": "A", "a2": "B"},
            {"person": "MO", "a1": "A", "a2": "C"},
            {"person": "CH", "a1": "A", "a2": "A"},
        ],
        "mutation_model": "Equal",
        "mutation_rate": 0.0,
        "kinship": 0.0,
    })
    assert r.status_code == 200
    body = r.json()
    # Textbook: LR = 1 / (2 p_A) = 2.5
    assert body["LR"] == pytest.approx(1 / (2 * 0.2), rel=1e-6)


def test_single_locus_one_parent_match(client: TestClient) -> None:
    r = client.post("/api/compute/single-locus", json={
        "mode": "one-parent",
        "locus": "L1",
        "frequencies": {"A": 0.2, "B": 0.3, "C": 0.5},
        "genotypes": [
            {"person": "AF", "a1": "A", "a2": "B"},
            {"person": "CH", "a1": "A", "a2": "A"},
        ],
        "mutation_model": "Equal",
        "mutation_rate": 0.0,
        "kinship": 0.0,
    })
    assert r.status_code == 200
    body = r.json()
    assert body["LR"] is not None and body["LR"] > 1.0


def test_single_locus_missing_genotypes(client: TestClient) -> None:
    r = client.post("/api/compute/single-locus", json={
        "mode": "one-parent",
        "locus": "L1",
        "frequencies": {"A": 0.5, "B": 0.5},
        "genotypes": [],
        "mutation_model": "Equal",
        "mutation_rate": 0.001,
        "kinship": 0.0,
    })
    assert r.status_code == 200
    body = r.json()
    assert body["LR"] is None and body["reason"]


# ---------------------------------------------------------------------------
# Full-case compute parity vs. tc.run_posterior
# ---------------------------------------------------------------------------
def _ref_one_parent(case: dict) -> tc.Result:
    h1, h2 = tc.make_duo_pedigrees()
    freq_db = {l["name"]: dict(l["frequencies"]) for l in case["loci"]}
    dna = [tc.TypedAllele(g["person"], l["name"], g["a1"], g["a2"])
           for l in case["loci"] for g in l["genotypes"]]
    loci_objs, dm, _w = tc.build_loci(
        [l["name"] for l in case["loci"]], freq_db, dna,
        mutation_model=case["mutation_model"],
        mutation_rate=case["mutation_rate"],
    )
    return tc.run_posterior(h1, h2, loci_objs, dm, kinship=case["kinship"])


def test_one_parent_matches_reference(client: TestClient) -> None:
    case = {
        "mutation_model": "Equal",
        "mutation_rate": 0.001,
        "kinship": 0.0,
        "prior": 0.5,
        "loci": [
            {
                "name": "L1",
                "frequencies": {"A": 0.2, "B": 0.3, "C": 0.5},
                "genotypes": [
                    {"person": "AF", "a1": "A", "a2": "B"},
                    {"person": "CH", "a1": "A", "a2": "C"},
                ],
            },
            {
                "name": "L2",
                "frequencies": {"X": 0.4, "Y": 0.6},
                "genotypes": [
                    {"person": "AF", "a1": "X", "a2": "Y"},
                    {"person": "CH", "a1": "X", "a2": "X"},
                ],
            },
        ],
    }
    ref = _ref_one_parent(case)
    r = client.post("/api/compute/one-parent", json=case)
    assert r.status_code == 200
    body = r.json()
    assert body["LR"] == pytest.approx(ref.LR, rel=1e-9)
    assert body["log10_LR"] == pytest.approx(ref.log10_LR, rel=1e-9)
    # Posterior with prior=0.5 must equal LR/(LR+1).
    assert body["posterior_h1"] == pytest.approx(
        ref.LR / (ref.LR + 1.0), rel=1e-9)
    assert len(body["per_locus"]) == 2


def test_two_parent_matches_reference(client: TestClient) -> None:
    case = {
        "mutation_model": "Equal",
        "mutation_rate": 0.001,
        "kinship": 0.0,
        "prior": 0.5,
        "loci": [
            {
                "name": "L1",
                "frequencies": {"A": 0.2, "B": 0.3, "C": 0.5},
                "genotypes": [
                    {"person": "AF", "a1": "A", "a2": "B"},
                    {"person": "MO", "a1": "A", "a2": "C"},
                    {"person": "CH", "a1": "A", "a2": "A"},
                ],
            },
        ],
    }
    h1, h2 = tc.make_trio_pedigrees()
    freq_db = {"L1": {"A": 0.2, "B": 0.3, "C": 0.5}}
    dna = [tc.TypedAllele(g["person"], "L1", g["a1"], g["a2"])
           for g in case["loci"][0]["genotypes"]]
    loci_objs, dm, _ = tc.build_loci(["L1"], freq_db, dna,
                                     mutation_rate=0.001)
    ref = tc.run_posterior(h1, h2, loci_objs, dm)

    r = client.post("/api/compute/two-parent", json=case)
    assert r.status_code == 200
    body = r.json()
    assert body["LR"] == pytest.approx(ref.LR, rel=1e-9)


def test_arbitrary_matches_two_parent(client: TestClient) -> None:
    """Express a trio paternity case through the arbitrary endpoint and
    verify it produces the same LR as the two-parent endpoint."""
    loci = [
        {
            "name": "L1",
            "frequencies": {"A": 0.2, "B": 0.3, "C": 0.5},
            "genotypes": [
                {"person": "AF", "a1": "A", "a2": "B"},
                {"person": "MO", "a1": "A", "a2": "C"},
                {"person": "CH", "a1": "A", "a2": "A"},
            ],
        },
    ]
    base = {"mutation_model": "Equal", "mutation_rate": 0.001,
            "kinship": 0.0, "prior": 0.5, "loci": loci}
    r2 = client.post("/api/compute/two-parent", json=base).json()
    arb = {
        **base,
        "persons": [
            {"id": "AF", "sex": "male"},
            {"id": "MO", "sex": "female"},
            {"id": "CH", "sex": "male"},
        ],
        "relations": [
            {"parent": "MO", "child": "CH", "flag": "fixed"},
            {"parent": "AF", "child": "CH", "flag": "test"},
        ],
    }
    r3 = client.post("/api/compute/arbitrary", json=arb).json()
    assert r3["LR"] == pytest.approx(r2["LR"], rel=1e-9)


def test_arbitrary_requires_test_relation(client: TestClient) -> None:
    body = {
        "mutation_model": "Equal", "mutation_rate": 0.001,
        "kinship": 0.0, "prior": 0.5,
        "loci": [{
            "name": "L1",
            "frequencies": {"A": 0.5, "B": 0.5},
            "genotypes": [
                {"person": "P1", "a1": "A", "a2": "B"},
                {"person": "P2", "a1": "A", "a2": "A"},
            ],
        }],
        "persons": [
            {"id": "P1", "sex": "male"},
            {"id": "P2", "sex": "male"},
        ],
        "relations": [
            {"parent": "P1", "child": "P2", "flag": "fixed"},
        ],
    }
    r = client.post("/api/compute/arbitrary", json=body)
    assert r.status_code == 400


# ---------------------------------------------------------------------------
# Frequency handling: preserve verbatim user/DB values + Rest top-up.
# ---------------------------------------------------------------------------
def test_frequencies_not_renormalised_when_sum_below_one(
        client: TestClient) -> None:
    """When user-supplied frequencies sum to < 1 (e.g. due to DB rounding),
    rare-allele values must be preserved verbatim and the missing mass put
    in a synthetic Rest bin — NOT smeared by renormalisation."""
    # Two cases: (a) verbatim frequencies summing to 0.9; (b) the same
    # frequencies post-renormalisation. They must give DIFFERENT LRs.
    case_verbatim = {
        "mutation_model": "Equal", "mutation_rate": 0.0,
        "kinship": 0.0, "prior": 0.5,
        "loci": [{
            "name": "L1",
            "frequencies": {"A": 0.18, "B": 0.27, "C": 0.45},  # sums to 0.9
            "genotypes": [
                {"person": "AF", "a1": "A", "a2": "B"},
                {"person": "MO", "a1": "A", "a2": "C"},
                {"person": "CH", "a1": "A", "a2": "A"},
            ],
        }],
    }
    case_renorm = {
        **case_verbatim,
        "loci": [{
            **case_verbatim["loci"][0],
            "frequencies": {"A": 0.2, "B": 0.3, "C": 0.5},  # 0.18/0.9 etc.
        }],
    }
    lr_v = client.post("/api/compute/two-parent", json=case_verbatim).json()["LR"]
    lr_r = client.post("/api/compute/two-parent", json=case_renorm).json()["LR"]
    # Textbook for renormalised case: 1 / (2 * 0.2) = 2.5
    assert lr_r == pytest.approx(2.5, rel=1e-9)
    # Verbatim case must use p_A = 0.18 → 1 / (2 * 0.18) ≈ 2.7778
    assert lr_v == pytest.approx(1 / (2 * 0.18), rel=1e-9)
    assert lr_v != pytest.approx(lr_r)


def test_single_locus_preserves_verbatim_frequencies(client: TestClient
                                                     ) -> None:
    r = client.post("/api/compute/single-locus", json={
        "mode": "two-parent",
        "locus": "L1",
        "frequencies": {"A": 0.18, "B": 0.27, "C": 0.45},  # sums to 0.9
        "genotypes": [
            {"person": "AF", "a1": "A", "a2": "B"},
            {"person": "MO", "a1": "A", "a2": "C"},
            {"person": "CH", "a1": "A", "a2": "A"},
        ],
        "mutation_model": "Equal",
        "mutation_rate": 0.0,
        "kinship": 0.0,
    }).json()
    assert r["LR"] == pytest.approx(1 / (2 * 0.18), rel=1e-9)


def test_observed_alleles_floored_when_missing_from_freqs(
        client: TestClient) -> None:
    """An observed allele not present in the frequency dict must be
    floored (not crash), and the rest of the dict preserved verbatim."""
    r = client.post("/api/compute/two-parent", json={
        "mutation_model": "Equal", "mutation_rate": 0.0,
        "kinship": 0.0, "prior": 0.5,
        "loci": [{
            "name": "L1",
            "frequencies": {"A": 0.5, "B": 0.4},  # C missing, sum 0.9
            "genotypes": [
                {"person": "AF", "a1": "A", "a2": "B"},
                {"person": "MO", "a1": "A", "a2": "C"},
                {"person": "CH", "a1": "A", "a2": "A"},
            ],
        }],
    })
    assert r.status_code == 200
    body = r.json()
    assert body["LR"] == pytest.approx(1 / (2 * 0.5), rel=1e-9)
    assert any("missing" in w.lower() for w in body["warnings"])


# ---------------------------------------------------------------------------
# Mendelian check
# ---------------------------------------------------------------------------
def test_mendelian_one_parent_match(client: TestClient) -> None:
    r = client.post("/api/check/mendelian", json={
        "mode": "one-parent",
        "parent": ["13", "14"],
        "child": ["13", "16"],
        "frequencies": {"13": 0.1, "16": 0.2},
    })
    assert r.status_code == 200
    body = r.json()
    assert body["mismatch"] is False
    assert body["power_of_exclusion"] is not None
    assert 0.0 < body["power_of_exclusion"] < 1.0


def test_mendelian_one_parent_mismatch(client: TestClient) -> None:
    r = client.post("/api/check/mendelian", json={
        "mode": "one-parent",
        "parent": ["13", "14"],
        "child": ["15", "16"],
    })
    assert r.json()["mismatch"] is True


def test_mendelian_two_parent(client: TestClient) -> None:
    r = client.post("/api/check/mendelian", json={
        "mode": "two-parent",
        "known": ["12", "14"],
        "alleged": ["16", "18"],
        "child": ["12", "20"],
    })
    assert r.json()["mismatch"] is True
