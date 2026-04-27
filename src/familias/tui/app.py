"""Textual TUI application — table-style design.

Run with::

    python -m familias.tui          # or `familias-tui`

Three modes:

1. **One-parent** — alleged parent + child.
2. **Two-parent** — known parent + alleged parent + child.
3. **Arbitrary** — N persons + an editable relationships table; flag one or
   more relations as ``test`` to define ``H1`` (with) vs ``H2`` (without).
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import math

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import (
    Horizontal, ScrollableContainer, Vertical, VerticalScroll,
)
from textual.message import Message
from textual.screen import ModalScreen, Screen
from textual.widget import Widget
from textual.widgets import (
    Button, Checkbox, Footer, Header, Input, Label, Select, Static,
)

from .. import FamiliasPedigree
from .compute import (
    _ensure_alleles, _resolve_freqs, is_mismatch_one_parent,
    is_mismatch_two_parents, make_duo_pedigrees, make_trio_pedigrees,
    make_pedigrees_from_relations, power_of_exclusion_one_parent,
    single_locus_lr, Person, Relation,
)
from .defaults import DEFAULT_LOCI
from .freq_loader import FreqDB, builtin_databases, load_freq_file
from .widgets import AmeloRow, CommitInput, FreqCell, LocusRow, ParamCell


# ---------------------------------------------------------------------------
# Path-prompt modal (used by Save / Load buttons)
# ---------------------------------------------------------------------------
class _ActionLabel(Static):
    """A clickable Static rendered like a markup label, e.g. ``[Save]``.

    Used in place of :class:`Button` because the default Button background /
    foreground combination renders as solid colour with invisible text under
    some Textual color schemes. Emits :class:`Pressed` when clicked or when
    Enter is pressed while it has focus.
    """

    DEFAULT_CSS = """
    _ActionLabel { width: auto; height: 1; padding: 0 1; color: $accent;
                   text-style: bold; }
    _ActionLabel:hover { background: $boost; }
    _ActionLabel:focus { background: $boost; text-style: bold reverse; }
    """

    can_focus = True

    class Pressed(Message):
        def __init__(self, sender_id: Optional[str]) -> None:
            self.sender_id = sender_id
            super().__init__()

    def __init__(self, label: str, *, id: Optional[str] = None) -> None:
        # Escape the leading ``[`` so Textual markup doesn't try to parse it.
        super().__init__(f"\\[{label}]", id=id)
        self._label = label

    def on_click(self) -> None:
        self.post_message(self.Pressed(self.id))

    def on_key(self, event) -> None:                           # noqa: ANN001
        if event.key in ("enter", "space"):
            event.stop()
            self.post_message(self.Pressed(self.id))


class _PathDialog(ModalScreen[Optional[str]]):
    DEFAULT_CSS = """
    _PathDialog { align: center middle; }
    _PathDialog #box { width: 100; height: auto; padding: 3 2;
                       border: round $accent; background: $surface; }
    _PathDialog #title { text-style: bold; padding-bottom: 1; }
    _PathDialog Input { width: 1fr; height: 3; border: none;
                        background: $boost; color: $foreground; padding: 0 1; }
    _PathDialog #buttons { height: 1; padding-top: 1; }
    _PathDialog Button { min-width: 10; margin-left: 1; height: 1;
                         border: none; }
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self, title: str, default: str = "") -> None:
        super().__init__()
        self._title = title
        self._default = default

    def compose(self) -> ComposeResult:
        with Vertical(id="box"):
            yield Static(self._title, id="title")
            self.input = Input(value=self._default, placeholder="path…",
                               id="path_input")
            yield self.input
            with Horizontal(id="buttons"):
                yield Static("", classes="grow")
                yield _ActionLabel("Cancel", id="cancel")
                yield _ActionLabel("OK", id="ok")

    def on_mount(self) -> None:
        self.input.focus()

    def on__action_label_pressed(self, event: "_ActionLabel.Pressed") -> None:  # noqa: N802
        if event.sender_id == "ok":
            self.dismiss(self.input.value.strip() or None)
        else:
            self.dismiss(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "ok":
            self.dismiss(self.input.value.strip() or None)
        else:
            self.dismiss(None)

    def on_input_submitted(self, event: Input.Submitted) -> None:
        self.dismiss(event.value.strip() or None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ---------------------------------------------------------------------------
# Top bar (DB picker + mutation parameters + file loader)
# ---------------------------------------------------------------------------
class TopBar(Vertical):
    DEFAULT_CSS = """
    TopBar { height: auto; padding: 0 1; }
    TopBar > Horizontal { height: 1; }
    TopBar > #db_row { height: 4; padding-bottom: 2; }
    TopBar Label { padding: 0 1; }
    TopBar .grow { width: 1fr; }
    TopBar #freq_path { width: 1fr; height: 1; border: none;
                        background: $boost; color: $foreground;
                        padding: 0 1; }
    TopBar #freq_path:focus { background: $surface; color: $foreground;
                              text-style: bold; }
    TopBar Select { width: 40; }
    TopBar SelectCurrent { color: $foreground; background: $boost; }
    TopBar SelectCurrent Static { color: $foreground; }
    TopBar #file_row.hidden { display: none; }
    TopBar .action_lbl { width: auto; height: 1; padding: 0 1;
                         color: $accent; text-style: bold; }
    TopBar .action_lbl:hover { background: $boost; }
    TopBar Button { height: 1; min-width: 8; border: none;
                    background: $primary; color: $foreground; }
    """

    class Changed(Message):
        """``kind`` is ``"db"`` when the database selection changed,
        ``"param"`` for mut.rate / theta edits.
        """
        def __init__(self, kind: str = "param") -> None:
            self.kind = kind
            super().__init__()

    class SaveRequested(Message):
        def __init__(self, path: str) -> None:
            self.path = path
            super().__init__()

    class LoadRequested(Message):
        def __init__(self, path: str) -> None:
            self.path = path
            super().__init__()

    def __init__(self, *, mode_label: str,
                 extras: Optional[List[Widget]] = None) -> None:
        super().__init__()
        self._mode_label = mode_label
        self._extras = extras or []

    def compose(self) -> ComposeResult:
        with Horizontal():
            yield Label("[b]Familias-py TUI[/]")
            yield Static("", classes="grow")
            yield Label(f"[b]{self._mode_label}[/]")
        with Horizontal(id="db_row"):
            yield Label("Allele Frequency Database:")
            self.databases = builtin_databases()
            opts: List[Tuple[str, str]] = [("Manual Input", "__manual__")]
            opts += [(name, name) for name in self.databases]
            opts.append(("Load from file\u2026", "__file__"))
            self.db_select = Select(opts, value="__manual__",
                                    allow_blank=False, id="db_select")
            yield self.db_select
            yield Static("", classes="grow")
            yield _ActionLabel("Save", id="save_btn")
            yield _ActionLabel("Load", id="load_btn")
            for w in self._extras:
                yield w
        with Horizontal():
            yield Label("Mut.Rate:")
            self.mut_rate = ParamCell("", "0.001", input_id="mut_rate")
            yield self.mut_rate
            yield Label("\u03b8:")
            self.kin = ParamCell("", "0.0", input_id="kin")
            yield self.kin
            yield Static("", classes="grow")
        with Horizontal(id="file_row", classes="hidden"):
            yield Label("Load file:")
            self.path_input = Input(
                placeholder="path to .csv / .txt / .json / .rda",
                id="freq_path",
            )
            yield self.path_input
            yield Button("Load", id="freq_load", variant="primary")
            self.status = Static("", id="db_status")
            yield self.status

    # events ---------------------------------------------------------------
    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "freq_load":
            event.stop()
            self._load_db_from_path()

    def on__action_label_pressed(self, event: "_ActionLabel.Pressed") -> None:  # noqa: N802
        if event.sender_id == "save_btn":
            event.stop()
            self.app.push_screen(
                _PathDialog("Save case to JSON", ""),
                callback=lambda p: p and self.post_message(
                    self.SaveRequested(p)),
            )
        elif event.sender_id == "load_btn":
            event.stop()
            self.app.push_screen(
                _PathDialog("Load case from JSON", ""),
                callback=lambda p: p and self.post_message(
                    self.LoadRequested(p)),
            )

    def _load_db_from_path(self) -> None:
        path = self.path_input.value.strip()
        if not path:
            self.status.update("[red]Enter a path first.[/]")
            return
        try:
            db = load_freq_file(path)
        except Exception as e:                           # noqa: BLE001
            self.status.update(f"[red]{e}[/]")
            return
        name = f"file:{Path(path).name}"
        self.databases[name] = db
        # Rebuild option list, keeping the "Load from file…" entry last.
        self.db_select.set_options(
            [("Manual Input", "__manual__")]
            + [(n, n) for n in self.databases]
            + [("Load from file…", "__file__")]
        )
        self.db_select.value = name
        self.status.update(f"[green]✓ {name} ({len(db)} loci)[/]")
        self.post_message(self.Changed("db"))

    def on_select_changed(self, event: Select.Changed) -> None:
        if event.select.id != "db_select":
            return
        event.stop()
        # Show the file-loader row only when the user picks "Load from file…".
        try:
            file_row = self.query_one("#file_row")
            file_row.set_class(event.value != "__file__", "hidden")
        except Exception:                                       # noqa: BLE001
            pass
        # Don't fire a recompute for the synthetic "__file__" sentinel —
        # nothing has changed yet on the data side.
        if event.value != "__file__":
            self.post_message(self.Changed("db"))

    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input.id in ("mut_rate", "kin"):
            event.stop()
            self.post_message(self.Changed("param"))

    # accessors ------------------------------------------------------------
    def db(self) -> Optional[FreqDB]:
        v = self.db_select.value
        if v in ("__manual__", "__file__") or v is Select.BLANK:
            return None
        return self.databases.get(str(v))

    def mut(self) -> float:
        try:
            return float(self.mut_rate.value)
        except ValueError:
            return 0.001

    def theta(self) -> float:
        try:
            return float(self.kin.value)
        except ValueError:
            return 0.0


# ---------------------------------------------------------------------------
# Combined-stats footer
# ---------------------------------------------------------------------------
class StatsFooter(Vertical):
    DEFAULT_CSS = """
    StatsFooter { height: auto; padding: 1 2; border-top: tall $accent; }
    StatsFooter Static { height: 1; }
    """

    def __init__(self, prior: float = 0.5) -> None:
        super().__init__()
        self.prior = prior
        self.lr_line = Static(
            "Combined Likelihood Ratio                       ?")
        self.post_line = Static(
            f"Post-test Probability (Pre-test={int(prior*100)}%)             ?"
        )

    def compose(self) -> ComposeResult:
        yield self.lr_line
        yield self.post_line

    def update_lr(self, lr: Optional[float]) -> None:
        if lr is None or not math.isfinite(lr) or lr <= 0:
            self.lr_line.update(
                "Combined Likelihood Ratio                       ?")
            self.post_line.update(
                f"Post-test Probability (Pre-test={int(self.prior*100)}%)             ?"
            )
            return
        post = (lr * self.prior) / (lr * self.prior + (1 - self.prior))
        self.lr_line.update(
            f"Combined Likelihood Ratio                       {lr:9.2f}")
        self.post_line.update(
            f"Post-test Probability (Pre-test={int(self.prior*100)}%)             "
            f"{post*100:.4f} %"
        )


# ---------------------------------------------------------------------------
# Base for paternity tables (modes 1 and 2)
# ---------------------------------------------------------------------------
class _PaternityBase(Screen):
    PERSONS: List[str] = []
    ALLEGED: str = "AF"
    KNOWN: Optional[str] = None
    CHILD: str = "CH"
    MODE_LABEL = ""
    AMELO_DEFAULTS: Dict[str, str] = {}

    BINDINGS = [
        Binding("escape", "app.pop_screen", "Back"),
        Binding("ctrl+r", "recompute_all", "Recompute"),
    ]

    DEFAULT_CSS = """
    Screen { layout: vertical; }
    #col_header { height: 1; padding: 0 1; background: $boost; }
    #col_header Static { content-align: left middle; }
    .h_tag { width: 7; }
    .h_loci { width: 10; }
    .h_geno { width: 12; }
    .h_lr { width: 10; content-align: right middle; }
    .h_freqs { width: 1fr; }
    #table_area { height: 1fr; padding: 0 1; }
    """

    def compose(self) -> ComposeResult:
        yield Header()
        self.topbar = TopBar(mode_label=self.MODE_LABEL)
        yield self.topbar
        with Horizontal(id="col_header"):
            yield Static("", classes="h_tag")
            yield Static("[b]Loci[/]", classes="h_loci")
            for header_text in self._geno_headers():
                yield Static(header_text, classes="h_geno")
            yield Static("[b]LR[/]", classes="h_lr")
            yield Static("[b]Frequencies[/]", classes="h_freqs")
        with VerticalScroll(id="table_area"):
            self.rows: Dict[str, LocusRow] = {}
            for locus in DEFAULT_LOCI:
                row = LocusRow(locus, self.PERSONS)
                self.rows[locus] = row
                yield row
            self.amelo = AmeloRow(self.PERSONS, defaults=self.AMELO_DEFAULTS)
            yield self.amelo
        self.stats = StatsFooter(prior=0.5)
        yield self.stats
        yield Footer()

    # subclasses fill in -------------------------------------------------
    def _geno_headers(self) -> List[str]:
        return [f"[b]{p}[/]" for p in self.PERSONS]

    def _pedigrees(self) -> Tuple[FamiliasPedigree, FamiliasPedigree]:
        raise NotImplementedError

    def _is_mismatch(self, gtypes: Dict[str, Optional[Tuple[str, str]]]) -> bool:
        return False

    # event plumbing -----------------------------------------------------
    def on_mount(self) -> None:
        self._lrs: Dict[str, Optional[float]] = {L: None for L in DEFAULT_LOCI}
        self._force_defaults = False
        self._update_combined()

    def on_locus_row_changed(self, event: LocusRow.Changed) -> None:
        self._recompute_locus(event.row)
        self._update_combined()

    def on_top_bar_changed(self, event: TopBar.Changed) -> None:
        # When the DB changed, overwrite any user-entered freqs with the
        # newly-selected database's values for this round of recompute.
        self._force_defaults = (event.kind == "db")
        try:
            self.action_recompute_all()
        finally:
            self._force_defaults = False

    def on_top_bar_save_requested(self, event: TopBar.SaveRequested) -> None:
        self._save_case(event.path)

    def on_top_bar_load_requested(self, event: TopBar.LoadRequested) -> None:
        self._load_case(event.path)

    def action_recompute_all(self) -> None:
        for row in self.rows.values():
            self._recompute_locus(row)
        self._update_combined()

    # core ---------------------------------------------------------------
    def _db_freqs(self, locus: str) -> Dict[str, float]:
        db = self.topbar.db()
        if not db:
            return {}
        try:
            _real, freqs = _resolve_freqs(locus, db)
            return dict(freqs)
        except KeyError:
            return {}

    def _recompute_locus(self, row: LocusRow) -> None:
        gtypes = row.get_genotypes()
        alleles_present: List[str] = []
        for g in gtypes.values():
            if g is None:
                continue
            for a in g:
                if a not in alleles_present:
                    alleles_present.append(a)

        if not alleles_present:
            row.freq.set_alleles([])
            row.set_lr("?")
            self._lrs[row.locus] = None
            return

        child_set = set(gtypes.get(self.CHILD) or ())
        # NR rule: an allele is "Not Relevant" if it does not appear in the
        # *child*'s genotype. The parent-genotype priors cancel between H1
        # and H2 in both the duo and trio LR, so only frequencies of the
        # child's observed alleles enter the result.
        def _nr(a: str) -> bool:
            return bool(child_set) and a not in child_set

        mismatch = self._is_mismatch(gtypes)
        db_freqs = self._db_freqs(row.locus)

        items: List[Tuple[str, bool, Optional[str]]] = []
        if mismatch:
            items.append(("MR", False, f"{self.topbar.mut():.4f}"))
            child = gtypes.get(self.CHILD)
            pe_val = (power_of_exclusion_one_parent(child, db_freqs)
                      if child is not None else None)
            items.append(("PE", False,
                          f"{pe_val:.4f}" if pe_val is not None else "?"))
        for a in alleles_present:
            items.append((a, _nr(a) and not mismatch, None))
        defaults = {a: db_freqs[a] for a in alleles_present if a in db_freqs}
        row.freq.set_alleles(items, defaults=defaults,
                             force_defaults=self._force_defaults)

        # Build engine-side frequency dict.
        user_freqs = row.freq.get_freqs()
        engine_freqs: Dict[str, float] = {}
        missing: List[str] = []
        for a in alleles_present:
            v = user_freqs.get(a)
            if v is None:
                v = db_freqs.get(a)
            if v is None:
                if _nr(a):
                    v = 0.01            # NR: irrelevant for LR
                else:
                    missing.append(a)
                    v = 0.01
            engine_freqs[a] = float(v)

        if missing:
            row.set_lr("?")
            self._lrs[row.locus] = None
            return

        try:
            full_freqs, _ = _ensure_alleles(engine_freqs, alleles_present)
            h1, h2 = self._pedigrees()
            typed = {p: g for p, g in gtypes.items() if g is not None}
            res = single_locus_lr(
                locus_name=row.locus,
                frequencies=full_freqs,
                genotypes=typed,
                h1=h1, h2=h2,
                mutation_rate=self.topbar.mut(),
                kinship=self.topbar.theta(),
            )
        except Exception:                                          # noqa: BLE001
            row.set_lr("[red]err[/]")
            self._lrs[row.locus] = None
            return
        if res is None:
            row.set_lr("?")
            self._lrs[row.locus] = None
            return
        l1, l2, lr = res
        if not math.isfinite(lr):
            row.set_lr("∞" if l1 > 0 else "0")
            self._lrs[row.locus] = None
            return
        row.set_lr(f"{lr:.4f}")
        self._lrs[row.locus] = lr

    def _update_combined(self) -> None:
        vals = [v for v in self._lrs.values()
                if v is not None and v > 0 and math.isfinite(v)]
        if not vals:
            self.stats.update_lr(None)
            return
        prod = 1.0
        for v in vals:
            prod *= v
        self.stats.update_lr(prod)

    # case I/O -----------------------------------------------------------
    def _save_case(self, path: str) -> None:
        data = {
            "schema": "familias-py-tui-case/1",
            "mode": self.MODE_LABEL,
            "persons": list(self.PERSONS),
            "mut_rate": self.topbar.mut(),
            "kinship": self.topbar.theta(),
            "database": str(self.topbar.db_select.value),
            "loci": {},
            "amelo": {p: self.amelo._inputs[p].value
                      for p in self.PERSONS},
        }
        for locus, row in self.rows.items():
            entry = {
                "genotypes": {p: row._inputs[p].value
                              for p in self.PERSONS},
                "freqs": {a: v for a, v in row.freq._values.items() if v},
            }
            data["loci"][locus] = entry
        try:
            with open(path, "w") as fh:
                json.dump(data, fh, indent=2)
            self.topbar.status.update(f"[green]✓ saved {path}[/]")
        except Exception as e:                                  # noqa: BLE001
            self.topbar.status.update(f"[red]save failed: {e}[/]")

    def _load_case(self, path: str) -> None:
        try:
            with open(path, "r") as fh:
                data = json.load(fh)
        except Exception as e:                                  # noqa: BLE001
            self.topbar.status.update(f"[red]load failed: {e}[/]")
            return
        try:
            self.topbar.mut_rate.value = str(data.get("mut_rate", 0.001))
            self.topbar.kin.value = str(data.get("kinship", 0.0))
            for p, v in (data.get("amelo") or {}).items():
                if p in self.amelo._inputs:
                    self.amelo._inputs[p].value = v
            for locus, entry in (data.get("loci") or {}).items():
                row = self.rows.get(locus)
                if row is None:
                    continue
                for p, g in (entry.get("genotypes") or {}).items():
                    if p in row._inputs:
                        row._inputs[p].value = g
                # stash freqs into the cell so that the next set_alleles
                # rebuild picks them up as the user value.
                for a, v in (entry.get("freqs") or {}).items():
                    row.freq._values[a] = str(v)
        except Exception as e:                                  # noqa: BLE001
            self.topbar.status.update(f"[red]load failed: {e}[/]")
            return
        self.topbar.status.update(f"[green]✓ loaded {path}[/]")
        self.action_recompute_all()


# ---------------------------------------------------------------------------
# Mode 1
# ---------------------------------------------------------------------------
class OneParentScreen(_PaternityBase):
    PERSONS = ["AF", "CH"]
    ALLEGED = "AF"
    CHILD = "CH"
    MODE_LABEL = "One-parent mode"
    AMELO_DEFAULTS = {"AF": "XY", "CH": ""}

    def _geno_headers(self) -> List[str]:
        return ["[b]Parent[/]", "[b]Child[/]"]

    def _pedigrees(self):
        return make_duo_pedigrees()

    def _is_mismatch(self, gtypes):
        return is_mismatch_one_parent(gtypes.get("AF"), gtypes.get("CH"))


# ---------------------------------------------------------------------------
# Mode 2
# ---------------------------------------------------------------------------
class TwoParentScreen(_PaternityBase):
    PERSONS = ["MO", "AF", "CH"]
    ALLEGED = "AF"
    KNOWN = "MO"
    CHILD = "CH"
    MODE_LABEL = "Two-parent mode"
    AMELO_DEFAULTS = {"MO": "XX", "AF": "XY", "CH": ""}

    def _geno_headers(self) -> List[str]:
        return ["[b]Known\nParent[/]", "[b]Alleged\nParent[/]", "[b]Child[/]"]

    def _pedigrees(self):
        return make_trio_pedigrees()

    def _is_mismatch(self, gtypes):
        return is_mismatch_two_parents(
            gtypes.get("MO"), gtypes.get("AF"), gtypes.get("CH"),
        )


# ---------------------------------------------------------------------------
# Arbitrary mode: N persons + relationships table
# ---------------------------------------------------------------------------
class _RelRow(Horizontal):
    DEFAULT_CSS = """
    _RelRow { height: 3; padding: 0 1; }
    _RelRow Label { width: 4; height: 3; content-align: right middle;
                    padding-right: 1; }
    _RelRow Select { width: 14; height: 3; }
    _RelRow Checkbox { width: 12; height: 3; }
    """

    class Changed(Message):
        pass

    def __init__(self, idx: int, persons: List[str]) -> None:
        super().__init__()
        self.idx = idx
        self.persons = persons

    def compose(self) -> ComposeResult:
        opts = [(p, p) for p in self.persons]
        self.parent_sel = Select(opts, prompt="?", id=f"rel_p_{self.idx}")
        self.child_sel = Select(opts, prompt="?", id=f"rel_c_{self.idx}")
        self.test_chk = Checkbox(value=False, id=f"rel_t_{self.idx}")
        yield Label(f"{self.idx}.")
        yield self.parent_sel
        yield self.child_sel
        yield self.test_chk

    def get(self) -> Optional[Tuple[str, str, bool]]:
        p = self.parent_sel.value
        c = self.child_sel.value
        if p is Select.BLANK or c is Select.BLANK or p == c:
            return None
        return str(p), str(c), bool(self.test_chk.value)

    def on_select_changed(self, event: Select.Changed) -> None:
        event.stop()
        self.post_message(self.Changed())

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        event.stop()
        self.post_message(self.Changed())


class ArbitraryScreen(Screen):
    MODE_LABEL = "Arbitrary Relationships mode"
    BINDINGS = [
        Binding("escape", "app.pop_screen", "Back"),
        Binding("ctrl+r", "recompute_all", "Recompute"),
        Binding("ctrl+a", "add_relation", "Add relation"),
    ]

    DEFAULT_CSS = _PaternityBase.DEFAULT_CSS + """
    #npersons { width: 6; height: 1; border: none; background: $boost;
                padding: 0 1; }
    #rel_area { height: 14; padding: 0 2; border-top: tall $accent; }
    #rel_header { height: 1; padding: 0 1; }
    """

    def compose(self) -> ComposeResult:
        yield Header()
        self.npersons_input = CommitInput(value="4", id="npersons")
        extras: List[Widget] = [
            Label("Number of individuals:"), self.npersons_input
        ]
        self.topbar = TopBar(mode_label=self.MODE_LABEL, extras=extras)
        yield self.topbar
        self.col_header = Horizontal(id="col_header")
        yield self.col_header
        self.table = VerticalScroll(id="table_area")
        yield self.table
        self.rel_area = VerticalScroll(id="rel_area")
        yield self.rel_area
        self.stats = StatsFooter(prior=0.5)
        yield self.stats
        yield Footer()

    # build ---------------------------------------------------------------
    def on_mount(self) -> None:
        self._lrs: Dict[str, Optional[float]] = {}
        self._force_defaults = False
        self._build_persons(4)

    def _build_persons(self, n: int) -> None:
        n = max(2, min(8, n))
        if getattr(self, "persons", None) == [f"P{i+1}" for i in range(n)]:
            return
        self.persons = [f"P{i+1}" for i in range(n)]
        # column header
        for c in list(self.col_header.children):
            c.remove()
        self.col_header.mount(Static("", classes="h_tag"))
        self.col_header.mount(Static("[b]Loci[/]", classes="h_loci"))
        for p in self.persons:
            self.col_header.mount(Static(f"[b]{p}[/]", classes="h_geno"))
        self.col_header.mount(Static("[b]LR[/]", classes="h_lr"))
        self.col_header.mount(Static("[b]Frequencies[/]", classes="h_freqs"))
        # table rows
        for c in list(self.table.children):
            c.remove()
        self.rows: Dict[str, LocusRow] = {}
        for locus in DEFAULT_LOCI:
            row = LocusRow(locus, self.persons)
            self.rows[locus] = row
            self.table.mount(row)
        self.amelo = AmeloRow(self.persons)
        self.table.mount(self.amelo)
        # relationships
        for c in list(self.rel_area.children):
            c.remove()
        self.rel_area.mount(Static(
            "[b]Relationships[/]   "
            "(Parent → Child; tick the rows defining the *test* hypothesis)"
        ))
        hdr = Horizontal()
        self.rel_area.mount(hdr)
        hdr.mount(Label(""))
        hdr.mount(Static("[b]Parent[/]"))
        hdr.mount(Static("[b]Child[/]"))
        hdr.mount(Static("[b]Test[/]"))
        self.rel_rows: List[_RelRow] = []
        for _ in range(3):
            self._add_rel_row()
        self._lrs = {L: None for L in DEFAULT_LOCI}
        self._update_combined()

    def _add_rel_row(self) -> None:
        idx = len(self.rel_rows) + 1
        r = _RelRow(idx, self.persons)
        self.rel_rows.append(r)
        self.rel_area.mount(r)

    def action_add_relation(self) -> None:
        self._add_rel_row()

    # events --------------------------------------------------------------
    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input.id == "npersons":
            event.stop()
            try:
                n = int(event.input.value)
            except ValueError:
                return
            self._build_persons(n)

    def on_top_bar_changed(self, event: TopBar.Changed) -> None:
        self._force_defaults = (event.kind == "db")
        try:
            self.action_recompute_all()
        finally:
            self._force_defaults = False

    def on_top_bar_save_requested(self, event: TopBar.SaveRequested) -> None:
        self._save_case(event.path)

    def on_top_bar_load_requested(self, event: TopBar.LoadRequested) -> None:
        self._load_case(event.path)

    def on_locus_row_changed(self, event: LocusRow.Changed) -> None:
        self._recompute_locus(event.row)
        self._update_combined()

    def on__rel_row_changed(self, event: _RelRow.Changed) -> None:  # noqa: N802
        self.action_recompute_all()

    def action_recompute_all(self) -> None:
        for row in self.rows.values():
            self._recompute_locus(row)
        self._update_combined()

    # compute -------------------------------------------------------------
    def _collect_relations(self) -> List[Relation]:
        out: List[Relation] = []
        for r in self.rel_rows:
            v = r.get()
            if v is None:
                continue
            p, c, test = v
            out.append(Relation(p, c, "test" if test else "fixed"))
        return out

    def _build_pedigrees(self) -> Optional[Tuple[FamiliasPedigree, FamiliasPedigree]]:
        rels = self._collect_relations()
        if not any(r.flag == "test" for r in rels):
            return None
        sexes = self.amelo.get_sexes()
        persons = [Person(p, sexes.get(p, "male")) for p in self.persons]
        try:
            return make_pedigrees_from_relations(persons, rels)
        except Exception:
            return None

    def _db_freqs(self, locus: str) -> Dict[str, float]:
        db = self.topbar.db()
        if not db:
            return {}
        try:
            _real, freqs = _resolve_freqs(locus, db)
            return dict(freqs)
        except KeyError:
            return {}

    def _recompute_locus(self, row: LocusRow) -> None:
        gtypes = row.get_genotypes()
        alleles_present: List[str] = []
        for g in gtypes.values():
            if g is None:
                continue
            for a in g:
                if a not in alleles_present:
                    alleles_present.append(a)
        if not alleles_present:
            row.freq.set_alleles([])
            row.set_lr("?")
            self._lrs[row.locus] = None
            return
        db_freqs = self._db_freqs(row.locus)
        row.freq.set_alleles(
            [(a, False, None) for a in alleles_present],
            defaults={a: db_freqs[a] for a in alleles_present if a in db_freqs},
            force_defaults=self._force_defaults,
        )
        peds = self._build_pedigrees()
        if peds is None:
            row.set_lr("—")
            self._lrs[row.locus] = None
            return
        h1, h2 = peds
        user_freqs = row.freq.get_freqs()
        engine_freqs: Dict[str, float] = {}
        missing: List[str] = []
        for a in alleles_present:
            v = user_freqs.get(a) or db_freqs.get(a)
            if v is None:
                missing.append(a)
                v = 0.01
            engine_freqs[a] = float(v)
        if missing:
            row.set_lr("?")
            self._lrs[row.locus] = None
            return
        typed = {p: g for p, g in gtypes.items() if g is not None}
        try:
            full_freqs, _ = _ensure_alleles(engine_freqs, alleles_present)
            res = single_locus_lr(
                locus_name=row.locus,
                frequencies=full_freqs,
                genotypes=typed,
                h1=h1, h2=h2,
                mutation_rate=self.topbar.mut(),
                kinship=self.topbar.theta(),
            )
        except Exception:                                       # noqa: BLE001
            row.set_lr("[red]err[/]")
            self._lrs[row.locus] = None
            return
        if res is None:
            row.set_lr("?")
            self._lrs[row.locus] = None
            return
        _, _, lr = res
        if not math.isfinite(lr):
            row.set_lr("—")
            self._lrs[row.locus] = None
            return
        row.set_lr(f"{lr:.4f}")
        self._lrs[row.locus] = lr

    def _update_combined(self) -> None:
        vals = [v for v in self._lrs.values()
                if v is not None and v > 0 and math.isfinite(v)]
        if not vals:
            self.stats.update_lr(None)
            return
        prod = 1.0
        for v in vals:
            prod *= v
        self.stats.update_lr(prod)

    # case I/O -----------------------------------------------------------
    def _save_case(self, path: str) -> None:
        data = {
            "schema": "familias-py-tui-case/1",
            "mode": self.MODE_LABEL,
            "persons": list(self.persons),
            "mut_rate": self.topbar.mut(),
            "kinship": self.topbar.theta(),
            "database": str(self.topbar.db_select.value),
            "amelo": {p: self.amelo._inputs[p].value for p in self.persons},
            "relations": [
                {"parent": str(r.parent_sel.value)
                          if r.parent_sel.value is not Select.BLANK else "",
                 "child": str(r.child_sel.value)
                          if r.child_sel.value is not Select.BLANK else "",
                 "test": bool(r.test_chk.value)}
                for r in self.rel_rows
            ],
            "loci": {},
        }
        for locus, row in self.rows.items():
            data["loci"][locus] = {
                "genotypes": {p: row._inputs[p].value
                              for p in self.persons},
                "freqs": {a: v for a, v in row.freq._values.items() if v},
            }
        try:
            with open(path, "w") as fh:
                json.dump(data, fh, indent=2)
            self.topbar.status.update(f"[green]✓ saved {path}[/]")
        except Exception as e:                                  # noqa: BLE001
            self.topbar.status.update(f"[red]save failed: {e}[/]")

    def _load_case(self, path: str) -> None:
        try:
            with open(path, "r") as fh:
                data = json.load(fh)
        except Exception as e:                                  # noqa: BLE001
            self.topbar.status.update(f"[red]load failed: {e}[/]")
            return
        try:
            self._build_persons(len(data.get("persons") or self.persons))
            self.topbar.mut_rate.value = str(data.get("mut_rate", 0.001))
            self.topbar.kin.value = str(data.get("kinship", 0.0))
            for p, v in (data.get("amelo") or {}).items():
                if p in self.amelo._inputs:
                    self.amelo._inputs[p].value = v
            for locus, entry in (data.get("loci") or {}).items():
                row = self.rows.get(locus)
                if row is None:
                    continue
                for p, g in (entry.get("genotypes") or {}).items():
                    if p in row._inputs:
                        row._inputs[p].value = g
                for a, v in (entry.get("freqs") or {}).items():
                    row.freq._values[a] = str(v)
        except Exception as e:                                  # noqa: BLE001
            self.topbar.status.update(f"[red]load failed: {e}[/]")
            return
        self.topbar.status.update(f"[green]✓ loaded {path}[/]")
        self.action_recompute_all()


# ---------------------------------------------------------------------------
# Home + App
# ---------------------------------------------------------------------------
class HomeScreen(Screen):
    DEFAULT_CSS = """
    HomeScreen { align: center middle; }
    #menu { width: 70; height: auto; border: round $accent; padding: 1 2; }
    #menu Label.title { text-style: bold; color: $accent; padding-bottom: 1; }
    #menu Button { width: 1fr; margin: 0 0 1 0; }
    """
    BINDINGS = [Binding("q", "app.quit", "Quit")]

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        with Vertical(id="menu"):
            yield Label("Familias-py — choose a mode", classes="title")
            yield Button("1)  One-parent  (alleged parent + child)",
                         id="mode_one", variant="primary")
            yield Button("2)  Two-parent  (known + alleged + child)",
                         id="mode_two", variant="primary")
            yield Button("3)  Arbitrary relationships hypothesis test",
                         id="mode_arb", variant="primary")
            yield Button("Quit", id="quit_btn")
        yield Footer()

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "mode_one":
            self.app.push_screen(OneParentScreen())
        elif event.button.id == "mode_two":
            self.app.push_screen(TwoParentScreen())
        elif event.button.id == "mode_arb":
            self.app.push_screen(ArbitraryScreen())
        elif event.button.id == "quit_btn":
            self.app.exit()


class FamiliasApp(App):
    TITLE = "Familias-py TUI"
    SUB_TITLE = "DNA pedigree LR / posterior calculator"
    BINDINGS = [Binding("ctrl+c", "quit", "Quit", priority=True)]

    def on_mount(self) -> None:
        self.push_screen(HomeScreen())


def main() -> None:
    FamiliasApp().run()


if __name__ == "__main__":  # pragma: no cover
    main()
