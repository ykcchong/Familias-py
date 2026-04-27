"""Default 15-locus panel used by the duo / trio modes.

Names follow the user-facing convention; :func:`canonical_locus_name`
maps them to the keys actually used in the bundled
:data:`familias.NorwegianFrequencies` dictionary (which uses ``TH01`` and
``VWA``).
"""
from __future__ import annotations

DEFAULT_LOCI = [
    "D8S1179",
    "D21S11",
    "D7S820",
    "CSF1PO",
    "D3S1358",
    "THO1",
    "D13S317",
    "D16S539",
    "D2S1338",
    "D19S433",
    "vWA",
    "TPOX",
    "D18S51",
    "D5S818",
    "FGA",
]


# Dye-channel groupings (used to colour the leftmost column of the TUI).
LOCUS_COLOR_GROUPS = [
    ("blue",  "Blue",  ["D8S1179", "D21S11", "D7S820", "CSF1PO"]),
    ("green", "Green", ["D3S1358", "THO1", "D13S317", "D16S539", "D2S1338"]),
    ("white", "Black", ["D19S433", "vWA", "TPOX", "D18S51"]),
    ("red",   "Red",   ["D5S818", "FGA"]),
]


def color_for_locus(locus: str) -> str:
    for color, _label, members in LOCUS_COLOR_GROUPS:
        if locus in members:
            return color
    return "white"


def is_first_in_group(locus: str) -> tuple[bool, str, str]:
    """Return ``(is_first, color, label)`` for the leftmost colour-tag column."""
    for color, label, members in LOCUS_COLOR_GROUPS:
        if members and members[0] == locus:
            return True, color, label
    return False, color_for_locus(locus), ""


# user-typed name -> canonical name (case-sensitive lookup keys in databases)
_ALIASES = {
    "THO1": "TH01",
    "TH0": "TH01",
    "vWA": "VWA",
    "VWA": "VWA",
    "PENTA E": "PENTA_E",
    "PENTA D": "PENTA_D",
    "Penta E": "PENTA_E",
    "Penta D": "PENTA_D",
}


def canonical_locus_name(name: str) -> str:
    """Return the canonical key for a locus name (case-/spelling-tolerant)."""
    n = name.strip()
    if n in _ALIASES:
        return _ALIASES[n]
    upper = n.upper().replace(" ", "_")
    return _ALIASES.get(upper, upper)
