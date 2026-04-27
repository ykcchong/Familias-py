"""Custom widgets for the Familias TUI.

The screens are built around a *table* metaphor: each row corresponds to one
locus and is composed of fixed-width cells (colour tag, locus name, one
``Input`` per typed individual, an LR readout, and a dynamic frequencies
cell). The frequencies cell is rebuilt whenever the genotype inputs change.
"""
from __future__ import annotations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.message import Message
from textual.widget import Widget
from textual.widgets import Input, Label, Static

from .defaults import is_first_in_group


# ---------------------------------------------------------------------------
# Genotype parsing
# ---------------------------------------------------------------------------
def parse_genotype(text: str) -> Optional[Tuple[str, str]]:
    """Parse '13', '13,15', '13/15', '13 15' → ``(a, b)``.

    Returns ``None`` for an empty / placeholder string.
    """
    s = text.strip()
    if not s or s == "?":
        return None
    for sep in (",", "/", "|", ";"):
        if sep in s:
            parts = [p.strip() for p in s.split(sep) if p.strip()]
            if len(parts) == 2:
                return parts[0], parts[1]
            if len(parts) == 1:
                return parts[0], parts[0]
            return None
    parts = s.split()
    if len(parts) == 2:
        return parts[0], parts[1]
    if len(parts) == 1:
        return parts[0], parts[0]
    return None


def fmt_genotype(g: Optional[Tuple[str, str]]) -> str:
    if g is None:
        return ""
    if g[0] == g[1]:
        return g[0]
    return f"{g[0]},{g[1]}"


def _safe_id(s: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in s)


# ---------------------------------------------------------------------------
# Input that commits on Enter *and* on focus loss (Tab / click away)
# ---------------------------------------------------------------------------
class CommitInput(Input):
    """An ``Input`` that re-emits :class:`Input.Submitted` whenever the user
    leaves the field — pressing Enter, Tab, or clicking elsewhere.

    Use this for fields that should trigger an expensive recompute only
    when the user is *done* typing, not on every keystroke.
    """

    def _on_blur(self, event) -> None:                       # noqa: ANN001
        super()._on_blur(event)
        self.post_message(self.Submitted(self, self.value, None))


# ---------------------------------------------------------------------------
# ParamCell: single labelled CommitInput, reuses FreqCell's `Input.freq`
# styling so the field is rendered identically to an allele frequency cell.
# ---------------------------------------------------------------------------
class ParamCell(Horizontal):
    DEFAULT_CSS = """
    ParamCell { height: 1; width: auto; }
    ParamCell .lab { width: 10; height: 1; content-align: right middle;
                     padding-right: 1; }
    ParamCell Input.freq { width: 12; height: 1; border: none; padding: 0 1;
                           background: $boost; color: $foreground; }
    ParamCell Input.freq:focus { background: $surface; color: $foreground;
                                 text-style: bold; }
    """

    def __init__(self, label: str, value: str, *, input_id: str) -> None:
        super().__init__()
        self._label = label
        self._value = value
        self._input_id = input_id

    def compose(self) -> ComposeResult:
        yield Label(self._label, classes="lab")
        self.input = CommitInput(value=self._value, classes="freq",
                                 id=self._input_id)
        yield self.input

    @property
    def value(self) -> str:
        return self.input.value

    @value.setter
    def value(self, v: str) -> None:
        self.input.value = v


# ---------------------------------------------------------------------------
# Frequencies cell
# ---------------------------------------------------------------------------
class FreqCell(Vertical):
    """Rebuildable cell of ``allele=[freq]`` editors (max ``PER_ROW`` per line).

    Use :meth:`set_alleles` to (re)populate. A *Not Relevant* allele is
    rendered as ``13=NR``; otherwise ``13=[Input]``. Pre-existing user
    edits are preserved across rebuilds when the allele still appears.

    Special "MR" / "PE" pseudo-alleles are used for mismatch display.
    """

    PER_ROW = 2

    DEFAULT_CSS = """
    FreqCell { height: auto; layout: vertical; }
    FreqCell > Horizontal { height: auto; width: auto; }
    FreqCell .lab { width: 6; height: 1; content-align: right middle; }
    FreqCell Input.freq { width: 12; height: 1; border: none; padding: 0 1;
                          background: $boost; color: $foreground; }
    FreqCell Input.freq:focus { background: $surface; color: $foreground;
                                text-style: bold; }
    FreqCell .nr { width: 12; height: 1; color: $text-muted; content-align: left middle; }
    FreqCell .ro { width: 12; height: 1; content-align: left middle; color: $warning; }
    """

    class Changed(Message):
        def __init__(self, sender: "FreqCell") -> None:
            self.sender = sender
            super().__init__()

    def __init__(self) -> None:
        super().__init__()
        self._values: Dict[str, str] = {}
        self._items: List[Tuple[str, bool, Optional[str]]] = []
        # (allele, is_NR, fixed_display) — fixed_display non-None ⇒ read-only

    # public API ------------------------------------------------------------
    def set_alleles(
        self,
        items: Sequence[Tuple[str, bool, Optional[str]]],
        defaults: Optional[Dict[str, float]] = None,
        *,
        force_defaults: bool = False,
    ) -> None:
        """Rebuild the cell with ``items``.

        Each item is ``(allele, is_NR, fixed_display)``. If ``fixed_display``
        is non-None it's shown as a non-editable Static (used for ``MR=...``
        on mismatches). ``defaults`` may pre-fill freq inputs.

        When ``force_defaults`` is True, ``defaults`` *overwrites* any
        previously-entered user value. Use this when the upstream database
        selection changes.
        """
        # snapshot current editable values before tearing down
        for a in [it[0] for it in self._items if not it[1] and it[2] is None]:
            v = self._read_input(a)
            if v is not None:
                self._values[a] = v
        if defaults:
            for a, v in defaults.items():
                if force_defaults:
                    self._values[a] = f"{v:.4f}"
                else:
                    self._values.setdefault(a, f"{v:.4f}")
        new_items = list(items)
        # If the cell already shows the same alleles in the same order with
        # the same NR/fixed flags, rebuild only when ``force_defaults`` asked
        # us to refresh values; otherwise keep the live Input untouched so
        # mid-edit focus / cursor are preserved.
        if new_items == self._items and not force_defaults:
            return
        self._items = new_items
        if self.is_mounted:
            self._rebuild()

    def set_user_value(self, allele: str, value: str) -> None:
        """Force a user value into the cell (used by case-file load)."""
        self._values[allele] = value
        if self.is_mounted:
            try:
                inp = self.query_one(f"#fq_{_safe_id(allele)}", Input)
                inp.value = value
            except Exception:                                  # noqa: BLE001
                pass

    def get_freqs(self) -> Dict[str, Optional[float]]:
        out: Dict[str, Optional[float]] = {}
        for a, nr, fixed in self._items:
            if nr or fixed is not None:
                continue
            txt = self._read_input(a)
            try:
                out[a] = float(txt) if txt else None
            except ValueError:
                out[a] = None
        return out

    # internals -------------------------------------------------------------
    def on_mount(self) -> None:
        self._rebuild()

    def _read_input(self, allele: str) -> Optional[str]:
        try:
            inp = self.query_one(f"#fq_{_safe_id(allele)}", Input)
        except Exception:
            return None
        return inp.value.strip() or None

    def _rebuild(self) -> None:
        for child in list(self.children):
            child.remove()
        if not self._items:
            return
        for chunk_start in range(0, len(self._items), self.PER_ROW):
            chunk = self._items[chunk_start:chunk_start + self.PER_ROW]
            row = Horizontal()
            self.mount(row)
            for allele, is_nr, fixed in chunk:
                row.mount(Label(f"{allele}=", classes="lab"))
                if fixed is not None:
                    row.mount(Static(fixed, classes="ro"))
                elif is_nr:
                    row.mount(Static("NR", classes="nr"))
                else:
                    row.mount(CommitInput(
                        value=self._values.get(allele, ""),
                        placeholder="?",
                        classes="freq",
                        id=f"fq_{_safe_id(allele)}",
                    ))

    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input.id and event.input.id.startswith("fq_"):
            event.stop()
            self.post_message(self.Changed(self))


# ---------------------------------------------------------------------------
# Locus row (one for each STR locus in the table)
# ---------------------------------------------------------------------------
class LocusRow(Horizontal):
    """A single locus row in the duo / trio / arbitrary tables.

    The row owns one ``Input`` per typed individual (``persons`` parameter),
    a read-only LR ``Static``, and a dynamic :class:`FreqCell`. It emits a
    :class:`LocusRow.Changed` message whenever any of its inputs change so
    the parent screen can recompute combined statistics.
    """

    DEFAULT_CSS = """
    LocusRow { height: auto; layout: horizontal; }
    LocusRow .tag { width: 7; content-align: left middle; text-style: bold; }
    LocusRow .name { width: 10; content-align: left middle; }
    LocusRow .geno { width: 12; height: 1; border: none; padding: 0 1;
                     background: $boost; color: $foreground; }
    LocusRow .geno:focus { background: $surface; color: $foreground;
                           text-style: bold; }
    LocusRow .lr { width: 10; content-align: right middle; padding-right: 1; }
    LocusRow FreqCell { width: 1fr; }
    """

    class Changed(Message):
        def __init__(self, row: "LocusRow") -> None:
            self.row = row
            super().__init__()

    def __init__(
        self,
        locus: str,
        persons: Sequence[str],
        *,
        show_lr: bool = True,
        show_color_tag: bool = True,
    ) -> None:
        super().__init__()
        self.locus = locus
        self.persons = list(persons)
        self.show_lr = show_lr
        self.show_color_tag = show_color_tag
        self._inputs: Dict[str, Input] = {}
        self.lr_widget: Optional[Static] = None
        self.freq: FreqCell = FreqCell()

    def compose(self) -> ComposeResult:
        if self.show_color_tag:
            first, color, label = is_first_in_group(self.locus)
            text = f"[{color}]{label}[/]" if (first and label) else ""
            yield Static(text, classes="tag")
        else:
            yield Static("", classes="tag")
        yield Static(self.locus, classes="name")
        for p in self.persons:
            inp = CommitInput(placeholder="?", classes="geno",
                              id=f"g_{_safe_id(self.locus)}_{_safe_id(p)}")
            self._inputs[p] = inp
            yield inp
        if self.show_lr:
            self.lr_widget = Static("?", classes="lr",
                                    id=f"lr_{_safe_id(self.locus)}")
            yield self.lr_widget
        else:
            yield Static("", classes="lr")
        yield self.freq

    # ---- public API ------------------------------------------------------
    def get_genotypes(self) -> Dict[str, Optional[Tuple[str, str]]]:
        return {p: parse_genotype(self._inputs[p].value) for p in self.persons}

    def get_freqs(self) -> Dict[str, Optional[float]]:
        return self.freq.get_freqs()

    def set_lr(self, text: str) -> None:
        if self.lr_widget is not None:
            self.lr_widget.update(text)

    def set_genotype(self, person: str, value: str) -> None:
        if person in self._inputs:
            self._inputs[person].value = value

    # ---- events ----------------------------------------------------------
    def on_input_submitted(self, event: Input.Submitted) -> None:
        if event.input in self._inputs.values():
            event.stop()
            self.post_message(self.Changed(self))

    def on_freq_cell_changed(self, event: FreqCell.Changed) -> None:
        event.stop()
        self.post_message(self.Changed(self))


# ---------------------------------------------------------------------------
# AMELO row — sex marker, no LR, no frequencies
# ---------------------------------------------------------------------------
class AmeloRow(Horizontal):
    DEFAULT_CSS = """
    AmeloRow { height: auto; layout: horizontal; }
    AmeloRow .tag { width: 7; }
    AmeloRow .name { width: 10; content-align: left middle; }
    AmeloRow .geno { width: 12; height: 1; border: none; padding: 0 1;
                     background: $boost; }
    AmeloRow .geno:focus { background: $surface; text-style: bold; }
    AmeloRow .pad { width: 1fr; }
    """

    def __init__(self, persons: Sequence[str], defaults: Optional[Dict[str, str]] = None) -> None:
        super().__init__()
        self.persons = list(persons)
        self._defaults = defaults or {}
        self._inputs: Dict[str, Input] = {}

    def compose(self) -> ComposeResult:
        yield Static("", classes="tag")
        yield Static("AMELO", classes="name")
        for p in self.persons:
            inp = CommitInput(value=self._defaults.get(p, ""),
                              placeholder="XY/XX", classes="geno",
                              id=f"amel_{_safe_id(p)}")
            self._inputs[p] = inp
            yield inp
        yield Static("", classes="pad")

    def get_sexes(self) -> Dict[str, str]:
        """Return ``{person: 'male'|'female'}`` (defaults to 'male' if unparseable)."""
        out: Dict[str, str] = {}
        for p, inp in self._inputs.items():
            v = inp.value.strip().upper()
            if v == "XX":
                out[p] = "female"
            elif v == "XY":
                out[p] = "male"
            else:
                out[p] = "male"
        return out
