"""Port of ``pater.cpp`` and the relevant parts of ``aldata.cpp`` /
``family.cpp`` that orchestrate a single computation pass.

The :class:`Pater` is a stateful container that holds a :class:`Family`,
a list of :class:`AlleleSystem` objects and (during a computation) an
:class:`Odds` object. Methods mirror the C++ public interface of
``pater``: ``add_person``, ``add_parent``, ``remove_person``,
``remove_possible_parent``, ``add_person_to_cutset``, ``end_cutset``,
``execute``, ``get_results``, ``remove_cutsets``.
"""
from __future__ import annotations
from typing import List, Optional
from .family import Family, Person
from .allele_system import AlleleSystem
from .peeling import Odds, SystemData


class Pater:
    def __init__(self) -> None:
        self.family = Family()
        self.systems: List[AlleleSystem] = []
        self._systems_by_name: dict[str, AlleleSystem] = {}
        self.recalculate: bool = True
        self.cuts_added: bool = False
        self.previously_used_kinship: float = -1.0
        self._odds: Optional[Odds] = None
        self._results: List[float] = []

    # ---- family ops --------------------------------------------------
    def add_person(self, name: str, male: bool) -> Person:
        self.recalculate = True
        return self.family.add_person(name, male)

    def remove_person(self, name: str) -> None:
        p = self.family.get(name)
        if p is None:
            return
        # Also drop any DNA observations referring to this person.
        for sys in self.systems:
            sys.remove_data(p)
        self.family.remove_person(p)
        self.recalculate = True

    def add_parent(self, parent_name: str, child_name: str) -> None:
        p = self.family.get(parent_name)
        c = self.family.get(child_name)
        if p is None or c is None:
            raise KeyError(f"Unknown person {parent_name!r} or {child_name!r}")
        self.family.add_relation(p, c)
        self.recalculate = True

    def remove_possible_parent(self, parent_name: str, child_name: str) -> None:
        p = self.family.get(parent_name)
        c = self.family.get(child_name)
        if p is None or c is None:
            return
        self.family.remove_relation(p, c)
        self.recalculate = True

    # ---- system ops --------------------------------------------------
    def add_system(self, sys: AlleleSystem) -> None:
        self._systems_by_name[sys.name] = sys
        self.systems.append(sys)
        self.recalculate = True

    def get_number_of_systems(self) -> int:
        return len(self.systems)

    def set_kinship(self, kinship: float) -> None:
        for sys in self.systems:
            sys.set_kinship(kinship)
            sys.recalc_data = True
        if kinship != self.previously_used_kinship:
            self.recalculate = True
        self.previously_used_kinship = kinship

    def add_data(self, system_name: str, person_name: str,
                 allele1_name: str, allele2_name: str) -> None:
        sys = self._systems_by_name[system_name]
        p = self.family.get(person_name)
        if p is None:
            raise KeyError(f"Unknown person {person_name!r}")
        a1 = sys.allele_names.index(allele1_name)
        a2 = sys.allele_names.index(allele2_name)
        sys.add_data(p, a1, a2)
        self.recalculate = True

    # ---- cutsets -----------------------------------------------------
    def _ensure_odds(self) -> None:
        if self._odds is None:
            # build the Odds object on the fly the first time we need it
            first = self.family.persons[0] if self.family.persons else None
            self._odds = Odds(first, separate_components=True)

    def add_person_to_cutset(self, name: str) -> None:
        self._ensure_odds()
        p = self.family.get(name)
        if p is None:
            raise KeyError(name)
        # The Pers alias was set when Odds was constructed.
        pe = p.alias
        if pe is None:
            raise RuntimeError(f"No Pers alias for {name!r}; build Odds first.")
        self.cuts_added = True
        self._odds.add_person_to_cutset(name, pe)

    def end_cutset(self) -> None:
        self._odds.end_cutset()

    # ---- run ---------------------------------------------------------
    def execute(self) -> None:
        # Build a fresh Odds object every time (mirrors C++ behaviour where
        # the odds object is rebuilt for each computation pass).
        if self._odds is None:
            first = self.family.persons[0] if self.family.persons else None
            self._odds = Odds(first, separate_components=True)
        results: List[float] = []
        for sys in self.systems:
            if sys.recalc_data:
                sys.compute_dataprob()
            sd = SystemData(
                name=sys.name,
                n_dataalleles=sys.n_dataalleles,
                dataprobability=sys.dataprobability,
                dataprobmatrix_female=sys.dataprobmatrix_female,
                dataprobmatrix_male=sys.dataprobmatrix_male,
                kinship=sys.kinship,
                has_silent_allele=sys.has_silent_allele,
                silent_allele=int(sys.index[sys.silent_allele]) if sys.has_silent_allele else -1,
            )
            # Inject person data into the Pers aliases.
            injected: List[object] = []
            for person, (a1, a2) in sys.data.items():
                if person.alias is None:
                    continue
                ra1 = int(sys.index[a1])
                ra2 = int(sys.index[a2])
                person.alias.add_data(sd, ra1, ra2)
                injected.append(person.alias)
            r = self._odds.execute(sd)
            for alias in injected:
                alias.remove_data()
            sys.result = r
            results.append(r)
        self._results = results
        self.recalculate = False

    def get_results(self) -> List[float]:
        return list(self._results)

    def remove_cutsets(self) -> None:
        # Drop the Odds object so the next execute() rebuilds it cleanly.
        # Also clear Pers aliases on persons.
        for p in self.family.persons:
            p.alias = None
            p.collapsed_alias = None
        self._odds = None

    def total_family_size(self) -> int:
        return len(self.family)

    def in_family(self, name: str) -> bool:
        return self.family.in_family(name)

    def is_male(self, name: str) -> bool:
        p = self.family.get(name)
        if p is None:
            raise KeyError(name)
        return p.male
