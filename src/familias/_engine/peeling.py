"""Port of ``odds.cpp`` and ``cutset.cpp`` — the peeling/cutset engine.

The original C++ uses an intrusive doubly-linked Linked_list class that
plays both the role of a list of items and an item that can be a member
of another list (``branch`` is both ``Link`` and ``Linked_list``). The port
keeps the same conceptual structure but uses ordinary Python lists; we
manipulate ``link.belongs_to`` (the owning list) and the ordering of
``items`` directly.

The classes are :class:`Pers`, :class:`Branch`, :class:`Cutset`. The
top-level :class:`Odds` builds a ``primcut`` (root cutset) holding one
``branch`` per connected component. ``add_person_to_cutset`` /
``end_cutset`` lets users (or the orchestrator) declare cutsets among the
persons of a branch.
"""
from __future__ import annotations
from typing import List, Optional
import numpy as np


# ---------------------------------------------------------------------------
# SystemData (port of class systemdata in odds.h)
# ---------------------------------------------------------------------------
class SystemData:
    """Per-system context shared by the peeling recursion."""

    __slots__ = (
        "name", "n_dataalleles", "dataprobability",
        "dataprobmatrix_female", "dataprobmatrix_male",
        "kinship", "has_silent_allele", "silent_allele",
        "_used", "_total_used",
    )

    def __init__(
        self,
        name: str,
        n_dataalleles: int,
        dataprobability: np.ndarray,
        dataprobmatrix_female: np.ndarray,
        dataprobmatrix_male: np.ndarray,
        kinship: float,
        has_silent_allele: bool,
        silent_allele: int,
    ):
        self.name = name
        self.n_dataalleles = n_dataalleles
        self.dataprobability = dataprobability
        self.dataprobmatrix_female = dataprobmatrix_female
        self.dataprobmatrix_male = dataprobmatrix_male
        self.kinship = float(kinship)
        self.has_silent_allele = bool(has_silent_allele)
        self.silent_allele = int(silent_allele)
        self._used = [0] * n_dataalleles
        self._total_used = 0

    def number_of_alleles(self) -> int:
        return self.n_dataalleles

    def prob(self, a: int) -> float:
        return self.dataprobability[a]

    def set_allele(self, a: int) -> float:
        used_a = self._used[a]
        total = self._total_used
        kin = self.kinship
        val = (used_a * kin + self.dataprobability[a] * (1.0 - kin)) / (1.0 + (total - 1) * kin)
        self._used[a] = used_a + 1
        self._total_used = total + 1
        return val

    def unset_allele(self, a: int) -> None:
        self._used[a] -= 1
        self._total_used -= 1

    def use_cutsets(self) -> bool:
        return self.kinship == 0.0

    def inherit_prob_female(self, paternal: int, maternal: int, result: int) -> float:
        m = self.dataprobmatrix_female
        return 0.5 * (m[paternal, result] + m[maternal, result])

    def inherit_prob_male(self, paternal: int, maternal: int, result: int) -> float:
        m = self.dataprobmatrix_male
        return 0.5 * (m[paternal, result] + m[maternal, result])


# ---------------------------------------------------------------------------
# Helpers: pcopy / pers / branch / cutset
# ---------------------------------------------------------------------------
class _Link:
    """Common base. Stores the list it currently belongs to."""

    def __init__(self) -> None:
        self.belongs_to: Optional["_LinkedList"] = None

    # The "execute" / "add_tables" etc. methods are overridden in subclasses.
    def container_branch(self) -> Optional["Branch"]:
        bt = self.belongs_to
        if isinstance(bt, Branch):
            return bt
        owner = getattr(bt, "owner", None)
        return owner if isinstance(owner, Branch) else None

    def container_cutset(self) -> Optional["Cutset"]:
        bt = self.belongs_to
        if isinstance(bt, Cutset):
            return bt
        owner = getattr(bt, "owner", None)
        return owner if isinstance(owner, Cutset) else None

    def has_data(self) -> bool:
        return False


class _LinkedList:
    """Ordered list of ``_Link`` objects.

    ``owner`` is an optional back-pointer to whatever object "owns" this
    list (used by :class:`Cutset` so that ``container_cutset()`` can
    follow ``belongs_to.owner``).
    """

    def __init__(self) -> None:
        self._items: List[_Link] = []
        self.owner = None

    def items(self) -> List[_Link]:
        return self._items

    def add(self, lk: _Link) -> None:
        # Mirrors ``Linked_list::add`` which prepends to the front.
        if lk.belongs_to is not None:
            lk.belongs_to.remove(lk)
        self._items.insert(0, lk)
        lk.belongs_to = self

    def add_at_end(self, lk: _Link) -> None:
        if lk.belongs_to is not None:
            lk.belongs_to.remove(lk)
        self._items.append(lk)
        lk.belongs_to = self

    def remove(self, lk: _Link) -> None:
        self._items.remove(lk)
        lk.belongs_to = None

    def get_first(self) -> Optional[_Link]:
        return self._items[0] if self._items else None

    def get_next(self, lk: _Link) -> Optional[_Link]:
        i = self._items.index(lk)
        return self._items[i + 1] if i + 1 < len(self._items) else None

    def n_elements(self) -> int:
        return len(self._items)

    def empty(self) -> bool:
        return not self._items


class Pcopy:
    """Reduced-pedigree copy of a :class:`Person`."""

    def __init__(self, person, is_collapsed: bool):
        self.alias = person
        self.male: bool = person.male
        self.mother: Optional[Pcopy] = None
        self.father: Optional[Pcopy] = None
        self.child: Optional[Pcopy] = None
        self.paternal_sibling: Optional[Pcopy] = None
        self.maternal_sibling: Optional[Pcopy] = None
        if is_collapsed:
            person.collapsed_alias = self
        else:
            person.alias = self

    def set_relatives(self) -> None:
        a = self.alias
        if a.mother is not None:
            self.mother = a.mother.alias
        if a.father is not None:
            self.father = a.father.alias
        if a.child is not None:
            self.child = a.child.alias
        if a.paternal_sibling is not None:
            self.paternal_sibling = a.paternal_sibling.alias
        if a.maternal_sibling is not None:
            self.maternal_sibling = a.maternal_sibling.alias

    def get_next_relative(self, p: Optional["Pcopy"]) -> Optional["Pcopy"]:
        if p is None:
            return self.mother or self.father or self.child
        if p is self.mother:
            return self.father or self.child
        if p is self.father:
            return self.child
        if self.male:
            return p.paternal_sibling
        return p.maternal_sibling


class Pers(_Link, Pcopy):
    """A person within the peeling data structure."""

    def __init__(self, person, is_collapsed: bool = False):
        _Link.__init__(self)
        Pcopy.__init__(self, person, is_collapsed)
        self._has_data: bool = False
        self.allele1: int = -1
        self.allele2: int = -1
        self.paternal_allele: int = -1
        self.maternal_allele: int = -1
        self.is_processed: bool = False

    def has_data(self) -> bool:
        return self._has_data

    def add_data(self, sd: SystemData, all1: int, all2: int) -> bool:
        if self._has_data:
            if not (
                (all1 == self.allele1 and all2 == self.allele2)
                or (all1 == self.allele2 and all2 == self.allele1)
            ):
                raise ValueError(
                    f"In system {sd.name!r} the odds computation "
                    "is incompatible with the given data."
                )
            return False
        self._has_data = True
        self.allele1 = all1
        self.allele2 = all2
        return True

    def remove_data(self) -> None:
        self._has_data = False

    def get_owner_branch(self) -> "Branch":
        cb = self.container_branch()
        if cb is not None:
            return cb
        return self.container_cutset().container_branch()

    def get_owner_cutset(self) -> "Cutset":
        cc = self.container_cutset()
        if cc is not None:
            return cc
        return self.container_branch().container_cutset()

    def collect_from(self, oldbranch: "Branch") -> None:
        thisbranch = self.get_owner_branch()
        p = None
        while True:
            p = self.get_next_relative(p)
            if p is None:
                break
            cc = p.container_cutset()
            if cc is thisbranch.container_cutset() and cc is not None:
                continue
            if p.get_owner_branch() is not oldbranch:
                continue
            lk = cc if cc is not None else p
            oldbranch.remove(lk)
            thisbranch.add(lk)
            lk.collect_from(oldbranch)

    # --- the heart of the peeling algorithm -----------------------------
    def execute(self, sd: SystemData) -> float:
        self.is_processed = True
        sum_ = 0.0
        N = sd.number_of_alleles()
        fath = self.father
        moth = self.mother
        for pa in range(N):
            for ma in range(N):
                if self._has_data and not (
                    (ma == self.allele1 and pa == self.allele2)
                    or (pa == self.allele1 and ma == self.allele2)
                ):
                    if not sd.has_silent_allele:
                        continue
                    if (pa != sd.silent_allele and ma != sd.silent_allele) or self.allele1 != self.allele2:
                        continue
                    if pa != self.allele1 and ma != self.allele1:
                        continue
                self.paternal_allele = pa
                self.maternal_allele = ma
                prob = 1.0
                set_pa = False
                set_ma = False
                if fath is not None:
                    if fath.is_processed:
                        prob *= sd.inherit_prob_male(fath.paternal_allele, fath.maternal_allele, pa)
                else:
                    prob *= sd.set_allele(pa)
                    set_pa = True
                if moth is not None:
                    if moth.is_processed:
                        prob *= sd.inherit_prob_female(moth.paternal_allele, moth.maternal_allele, ma)
                else:
                    prob *= sd.set_allele(ma)
                    set_ma = True
                # Loop over already-processed children
                p = self.child
                if self.male:
                    while p is not None:
                        if p.is_processed:
                            prob *= sd.inherit_prob_male(pa, ma, p.paternal_allele)
                        p = p.paternal_sibling
                else:
                    while p is not None:
                        if p.is_processed:
                            prob *= sd.inherit_prob_female(pa, ma, p.maternal_allele)
                        p = p.maternal_sibling

                if prob > 0:
                    lk = self.container_branch().get_next(self)
                    if lk is not None:
                        prob *= lk.execute(sd)
                    sum_ += prob
                if set_pa:
                    sd.unset_allele(pa)
                if set_ma:
                    sd.unset_allele(ma)
        self.is_processed = False
        return sum_

    def execute_cutset_part(self, sd: SystemData, index: int) -> float:
        self.is_processed = True
        sum_ = 0.0
        N = sd.number_of_alleles()
        fath = self.father
        moth = self.mother
        index = index * N * N
        for pa in range(N):
            for ma in range(N):
                if self._has_data and not (
                    (ma == self.allele1 and pa == self.allele2)
                    or (pa == self.allele1 and ma == self.allele2)
                ):
                    if not sd.has_silent_allele:
                        continue
                    if (pa != sd.silent_allele and ma != sd.silent_allele) or self.allele1 != self.allele2:
                        continue
                    if pa != self.allele1 and ma != self.allele1:
                        continue
                self.paternal_allele = pa
                self.maternal_allele = ma
                prob = 1.0
                set_pa = False
                set_ma = False
                if fath is not None:
                    if fath.is_processed:
                        prob *= sd.inherit_prob_male(fath.paternal_allele, fath.maternal_allele, pa)
                else:
                    prob *= sd.set_allele(pa)
                    set_pa = True
                if moth is not None:
                    if moth.is_processed:
                        prob *= sd.inherit_prob_female(moth.paternal_allele, moth.maternal_allele, ma)
                else:
                    prob *= sd.set_allele(ma)
                    set_ma = True
                p = self.child
                if self.male:
                    while p is not None:
                        if p.is_processed:
                            prob *= sd.inherit_prob_male(pa, ma, p.paternal_allele)
                        p = p.paternal_sibling
                else:
                    while p is not None:
                        if p.is_processed:
                            prob *= sd.inherit_prob_female(pa, ma, p.maternal_allele)
                        p = p.maternal_sibling
                if prob > 0:
                    cs = self.container_cutset()
                    pr = cs.get_next_pers(self)
                    new_index = index + N * pa + ma
                    if pr is not None:
                        prob *= pr.execute_cutset_part(sd, new_index)
                    else:
                        prob *= cs.execute_cutset(sd, new_index)
                    sum_ += prob
                if set_pa:
                    sd.unset_allele(pa)
                if set_ma:
                    sd.unset_allele(ma)
        self.is_processed = False
        return sum_

    def add_tables(self, n_alleles: int) -> bool:
        return False

    def remove_tables(self) -> None:
        pass

    def sort(self) -> None:
        pass


class Branch(_Link, _LinkedList):
    """A branch is itself a Link (lives in a cutset's branch list) and a list."""

    def __init__(self) -> None:
        _Link.__init__(self)
        _LinkedList.__init__(self)

    def has_data(self) -> bool:
        return False

    def sort(self) -> None:
        # Move data-bearing pers to the front (preserving relative order).
        items = self._items
        must_move = False
        i = 0
        while i < len(items):
            lk = items[i]
            lk.sort()
            if lk.has_data():
                if must_move:
                    items.pop(i)
                    items.insert(0, lk)
                    # i stays the same since we moved one out and inserted at 0
                    # but the new item at i is the next original element.
                    continue
            else:
                must_move = True
            i += 1

    def add_tables(self, n_alleles: int) -> bool:
        for lk in self._items:
            if lk.add_tables(n_alleles):
                return True
        return False

    def remove_tables(self) -> None:
        for lk in self._items:
            lk.remove_tables()

    def remove_data(self) -> None:
        for lk in self._items:
            lk.remove_data()


# A "branch_list" / "pers_list" in C++ is just a list of Branch / Pers; we use
# simple Python lists for these inside Cutset.

class Cutset(_Link):
    """A cutset is a Link, plus a list of pers (its members) and a list of branches."""

    def __init__(self) -> None:
        _Link.__init__(self)
        self._pers_list = _LinkedList()   # type: ignore[assignment]
        self._branch_list = _LinkedList() # type: ignore[assignment]
        self._pers_list.owner = self
        self._branch_list.owner = self
        self._tab: Optional[np.ndarray] = None

    # mimic the dual nature of the C++ multi-inheritance:
    def get_first_pers(self) -> Optional[Pers]:
        first = self._pers_list.get_first()
        return first  # type: ignore[return-value]

    def get_next_pers(self, p: Pers) -> Optional[Pers]:
        return self._pers_list.get_next(p)  # type: ignore[return-value]

    def pers_add(self, p: Pers) -> None:
        self._pers_list.add(p)

    def pers_remove(self, p: Pers) -> None:
        self._pers_list.remove(p)

    def empty_pers(self) -> bool:
        return self._pers_list.empty()

    def get_first_branch(self) -> Optional[Branch]:
        return self._branch_list.get_first()  # type: ignore[return-value]

    def get_next_branch(self, b: Branch) -> Optional[Branch]:
        return self._branch_list.get_next(b)  # type: ignore[return-value]

    def branch_add(self, b: Branch) -> None:
        self._branch_list.add(b)

    def branch_remove(self, b: Branch) -> None:
        self._branch_list.remove(b)

    def empty_branch(self) -> bool:
        return self._branch_list.empty()

    def n_branches(self) -> int:
        return self._branch_list.n_elements()

    def n_pers(self) -> int:
        return self._pers_list.n_elements()

    def has_data(self) -> bool:
        return False

    def sort(self) -> None:
        for b in self._branch_list.items():
            b.sort()  # type: ignore[union-attr]

    def separate_branches(self) -> None:
        """Split this cutset's single branch into separate connected components.

        We assume the cutset currently has exactly one branch.
        """
        oldbranch = self.get_first_branch()
        assert oldbranch is not None
        while not oldbranch.empty():
            br = Branch()
            self.branch_add(br)
            lk = oldbranch.get_first()
            oldbranch.remove(lk)
            br.add(lk)
            lk.collect_from(oldbranch)  # type: ignore[union-attr]
        self.branch_remove(oldbranch)

    def find_relative_in_branch(self, br: Branch):
        p = self.get_first_pers()
        while p is not None:
            q = None
            while True:
                q = p.get_next_relative(q)
                if q is None:
                    break
                if q.get_owner_branch() is br:
                    cc = q.container_cutset()
                    return cc if cc is not None else q
            p = self.get_next_pers(p)
        return None

    def collect_from(self, oldbranch: Branch) -> None:
        p = self.get_first_pers()
        while p is not None:
            p.collect_from(oldbranch)
            p = self.get_next_pers(p)

    def add_tables(self, n_alleles: int) -> bool:
        from .special import MAX_TABLE_SIZE
        tablesize = 1
        for _ in range(2 * self.n_pers()):
            if tablesize > MAX_TABLE_SIZE // n_alleles:
                return True
            tablesize *= n_alleles
        self._tab = np.full(tablesize, -1.0, dtype=float)
        for b in self._branch_list.items():
            if b.add_tables(n_alleles):  # type: ignore[union-attr]
                return True
        return False

    def remove_tables(self) -> None:
        self._tab = None
        for b in self._branch_list.items():
            b.remove_tables()  # type: ignore[union-attr]

    def remove_data(self) -> None:
        p = self.get_first_pers()
        while p is not None:
            p.remove_data()
            p = self.get_next_pers(p)
        for b in self._branch_list.items():
            b.remove_data()  # type: ignore[union-attr]

    def execute(self, sd: SystemData) -> float:
        # Called only for non-primcut cutsets — primcut is handled by Odds.
        return self.get_first_pers().execute_cutset_part(sd, 0)

    def execute_cutset(self, sd: SystemData, index: int) -> float:
        if self._tab[index] < 0:
            result = 1.0
            for b in self._branch_list.items():
                result *= b.get_first().execute(sd)  # type: ignore[union-attr]
            self._tab[index] = result
        nxt = self.container_branch().get_next(self)
        if nxt is not None:
            return self._tab[index] * nxt.execute(sd)
        return float(self._tab[index])


# ---------------------------------------------------------------------------
# Odds (port of class odds in odds.h / odds.cpp)
# ---------------------------------------------------------------------------
class Odds:
    def __init__(self, family_first_person, separate_components: bool):
        self.primcut = Cutset()
        self.currcut: Optional[Cutset] = None
        self.currbranch: Optional[Branch] = None
        self.cutset_must_end: bool = False
        br = Branch()
        self.primcut.branch_add(br)
        p = family_first_person
        while p is not None:
            br.add_at_end(Pers(p, is_collapsed=False))
            p = p.next
        # set relatives
        for lk in br.items():
            lk.set_relatives()  # type: ignore[union-attr]
        self.collapsed_pers: Optional[Pers] = None
        if separate_components:
            self.primcut.separate_branches()

    # --- adding cutsets ------------------------------------------------
    def add_person_to_cutset(self, name: str, pe: Pers) -> bool:
        if self.cutset_must_end:
            raise ValueError(
                "Cutset persons cannot belong to two different cutsets "
                "unless the last contains only one person."
            )
        if pe.container_cutset() is not None:
            if pe.container_cutset() is self.currcut:
                return False
            if pe is self.collapsed_pers:
                if self.currcut is not None:
                    raise ValueError(
                        "Cutset persons cannot belong to two different cutsets."
                    )
                self.cutset_must_end = True
                return False
            raise ValueError(f"The person {name!r} cannot belong to two cutsets.")
        if self.currcut is None:
            self.currcut = Cutset()
            self.currbranch = pe.container_branch()
        else:
            if pe.container_branch() is not self.currbranch:
                raise ValueError(
                    f"Person {name!r} is not in the same component as the rest of the cutset."
                )
        self.currbranch.remove(pe)
        self.currcut.pers_add(pe)
        return False

    def end_cutset(self) -> None:
        if self.cutset_must_end:
            self.cutset_must_end = False
            return
        if self.currbranch.empty():
            # Unnecessary cutset: put pers back into branch.
            cc = self.currcut
            pr = cc.get_first_pers()
            while pr is not None:
                cc.pers_remove(pr)
                self.currbranch.add(pr)
                pr = cc.get_first_pers()
            self.currcut = None
            self.currbranch = None
            return
        oldcutset = self.currbranch.container_cutset()
        oldcutset.branch_remove(self.currbranch)
        newbranch = Branch()
        oldcutset.branch_add(newbranch)
        newbranch.add(self.currcut)
        if oldcutset is self.primcut:
            lk = self.currbranch.get_first()
        else:
            lk = oldcutset.find_relative_in_branch(self.currbranch)
        if lk is None:
            self.currcut.branch_add(self.currbranch)
            self.currcut.separate_branches()
        else:
            while lk is not None:
                self.currbranch.remove(lk)
                newbranch.add(lk)
                lk.collect_from(self.currbranch)
                if oldcutset is self.primcut:
                    break
                lk = oldcutset.find_relative_in_branch(self.currbranch)
            if self.currbranch.empty():
                pr = self.currcut.get_first_pers()
                while pr is not None:
                    self.currcut.pers_remove(pr)
                    newbranch.add(pr)
                    pr = self.currcut.get_first_pers()
                newbranch.remove(self.currcut)
            else:
                self.currcut.branch_add(self.currbranch)
                self.currcut.separate_branches()
        self.currcut = None
        self.currbranch = None

    # --- data --------------------------------------------------------
    def add_data(self, sd: SystemData, p: Pers, allele1: int, allele2: int) -> bool:
        return p.add_data(sd, allele1, allele2)

    def remove_data(self) -> None:
        self.primcut.remove_data()

    # --- run ---------------------------------------------------------
    def execute(self, sd: SystemData) -> float:
        if self.primcut.add_tables(sd.number_of_alleles()):
            self.primcut.remove_tables()
            raise MemoryError("cutset memoization table too large")
        self.primcut.sort()
        result = 1.0
        for br in self.primcut._branch_list.items():
            first = br.get_first()  # type: ignore[union-attr]
            result *= first.execute(sd)
        self.primcut.remove_tables()
        return result
