"""Smoke tests for the TUI compute layer + headless app run."""
import asyncio
from familias.tui.compute import (
    is_mismatch_one_parent, is_mismatch_two_parents,
    power_of_exclusion_one_parent, single_locus_lr,
    make_duo_pedigrees, make_trio_pedigrees,
)
from familias.tui.freq_loader import (
    builtin_databases, load_cap_txt, load_fsigen_csv,
)
from familias.tui.widgets import parse_genotype, fmt_genotype


def test_builtin_databases_includes_norwegian():
    dbs = builtin_databases()
    assert any(k.startswith("NorwegianFrequencies") for k in dbs)


def test_cap_loader_parses_blocks():
    db = load_cap_txt(
        "/Users/calvinchong/projects/forensicpopdata/CAPdata/2024b/"
        "Dry challenge allele frequency database.txt"
    )
    assert "CSF1PO" in db
    assert 0.999 < sum(db["CSF1PO"].values()) < 1.001


def test_fsigen_loader():
    db = load_fsigen_csv(
        "/Users/calvinchong/projects/forensicpopdata/inst/extdata/"
        "FBI_extended_Cauc_022024.csv"
    )
    assert "TPOX" in db


def test_parse_genotype_variants():
    assert parse_genotype("13") == ("13", "13")
    assert parse_genotype("13,15") == ("13", "15")
    assert parse_genotype("13/15") == ("13", "15")
    assert parse_genotype("13 15") == ("13", "15")
    assert parse_genotype("?") is None
    assert parse_genotype("") is None


def test_fmt_genotype():
    assert fmt_genotype(("13", "13")) == "13"
    assert fmt_genotype(("13", "15")) == "13,15"
    assert fmt_genotype(None) == ""


def test_one_parent_mismatch():
    assert is_mismatch_one_parent(("13", "14"), ("15", "16")) is True
    assert is_mismatch_one_parent(("13", "14"), ("13", "16")) is False
    assert is_mismatch_one_parent(None, ("15", "16")) is False


def test_two_parent_mismatch():
    assert is_mismatch_two_parents(("12", "14"), ("16", "18"),
                                   ("12", "16")) is False
    assert is_mismatch_two_parents(("12", "14"), ("16", "18"),
                                   ("12", "20")) is True


def test_power_of_exclusion_in_unit_range():
    pe = power_of_exclusion_one_parent(
        ("13", "15"), {"13": 0.10, "15": 0.20}
    )
    assert pe is not None and 0.0 < pe < 1.0


def test_single_locus_lr_textbook_paternity():
    h1, h2 = make_trio_pedigrees()
    res = single_locus_lr(
        locus_name="L1",
        frequencies={"A": 0.2, "B": 0.3, "C": 0.5},
        genotypes={"AF": ("A", "B"), "MO": ("A", "C"), "CH": ("A", "A")},
        h1=h1, h2=h2,
        mutation_rate=0.0,
    )
    assert res is not None
    _l1, _l2, lr = res
    assert abs(lr - 1 / (2 * 0.2)) < 1e-6


def test_single_locus_lr_duo_with_match():
    h1, h2 = make_duo_pedigrees()
    res = single_locus_lr(
        locus_name="L1",
        frequencies={"A": 0.2, "B": 0.3, "C": 0.5},
        genotypes={"AF": ("A", "B"), "CH": ("A", "A")},
        h1=h1, h2=h2,
        mutation_rate=0.0,
    )
    _l1, _l2, lr = res
    assert lr > 1.0


def test_tui_screens_mount_and_pop():
    from familias.tui.app import FamiliasApp

    async def go():
        app = FamiliasApp()
        async with app.run_test() as pilot:
            await pilot.pause()
            for sid in ("mode_one", "mode_two", "mode_arb"):
                await pilot.click("#" + sid)
                await pilot.pause()
                # escape is captured by focused Inputs; pop programmatically
                app.pop_screen()
                await pilot.pause()
    asyncio.run(go())


def test_tui_one_parent_live_lr():
    """Enter genotypes in the one-parent screen and verify a finite LR appears."""
    from familias.tui.app import FamiliasApp, OneParentScreen
    from familias.tui.widgets import LocusRow

    async def go():
        app = FamiliasApp()
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.click("#mode_one")
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, OneParentScreen)
            # Pick a built-in DB
            screen.topbar.db_select.value = "NorwegianFrequencies (built-in)"
            await pilot.pause()
            row = screen.rows["TPOX"]
            row._inputs["AF"].value = "8,9"
            row._inputs["AF"].post_message(
                row._inputs["AF"].Submitted(row._inputs["AF"], "8,9", None))
            await pilot.pause()
            row._inputs["CH"].value = "8,10"
            row._inputs["CH"].post_message(
                row._inputs["CH"].Submitted(row._inputs["CH"], "8,10", None))
            await pilot.pause()
            # The LR widget should now show a numeric value.
            assert row.lr_widget is not None
            text = str(row.lr_widget.render())
            assert text not in ("?", "—")
            float(text.replace("[red]err[/]", "nan"))
    asyncio.run(go())
