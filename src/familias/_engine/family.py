"""Port of the C++ ``person`` and ``family`` classes.

Each ``Person`` keeps direct ``mother``/``father`` references; the engine
also maintains ``children``/``paternal_sibling``/``maternal_sibling`` pointers
that mirror the C++ doubly-threaded sibling lists used by the peeling code.

The :class:`Family` is just an ordered linked list of :class:`Person` objects
(insertion-ordered like the original).
"""
from __future__ import annotations
from typing import Optional, List


class Person:
    __slots__ = (
        "name", "male",
        "mother", "father", "child",
        "paternal_sibling", "maternal_sibling",
        "next",
        "alias", "collapsed_alias",
    )

    def __init__(self, name: str, male: bool):
        self.name: str = name
        self.male: bool = bool(male)
        self.mother: Optional[Person] = None
        self.father: Optional[Person] = None
        self.child: Optional[Person] = None
        self.paternal_sibling: Optional[Person] = None
        self.maternal_sibling: Optional[Person] = None
        self.next: Optional[Person] = None
        self.alias = None              # set by pcopy / Pers
        self.collapsed_alias = None

    # --- ancestor check (person.cpp::has_ancestor) -----------------------
    def has_ancestor(self, p: "Person") -> bool:
        if self is p:
            return True
        stack = [self.mother, self.father]
        while stack:
            curr = stack.pop()
            if curr is None:
                continue
            if curr is p:
                return True
            stack.append(curr.mother)
            stack.append(curr.father)
        return False

    # --- relationship maintenance ---------------------------------------
    def add_parent(self, p: "Person") -> None:
        """Attach ``p`` as a parent (caller has verified parent does not exist)."""
        if p.male:
            self.father = p
            self.paternal_sibling = p.child
            p.child = self
        else:
            self.mother = p
            self.maternal_sibling = p.child
            p.child = self

    def remove_mother(self) -> None:
        m = self.mother
        if m is None:
            return
        # remove self from m.child / maternal_sibling chain
        if m.child is self:
            m.child = self.maternal_sibling
        else:
            q = m.child
            while q is not None and q.maternal_sibling is not self:
                q = q.maternal_sibling
            if q is not None:
                q.maternal_sibling = self.maternal_sibling
        self.mother = None
        self.maternal_sibling = None

    def remove_father(self) -> None:
        f = self.father
        if f is None:
            return
        if f.child is self:
            f.child = self.paternal_sibling
        else:
            q = f.child
            while q is not None and q.paternal_sibling is not self:
                q = q.paternal_sibling
            if q is not None:
                q.paternal_sibling = self.paternal_sibling
        self.father = None
        self.paternal_sibling = None


class Family:
    """Ordered collection of persons (preserves insertion order)."""

    def __init__(self) -> None:
        self.persons: List[Person] = []
        self._by_name: dict[str, Person] = {}

    def __iter__(self):
        return iter(self.persons)

    def __len__(self):
        return len(self.persons)

    def get(self, name: str) -> Optional[Person]:
        return self._by_name.get(name)

    def in_family(self, name: str) -> bool:
        return name in self._by_name

    def add_person(self, name: str, male: bool) -> Person:
        if name in self._by_name:
            raise ValueError(f"Person {name!r} already in family")
        p = Person(name, male)
        # Append at end and link via the C++-style 'next' pointer too.
        if self.persons:
            self.persons[-1].next = p
        self.persons.append(p)
        self._by_name[name] = p
        return p

    def remove_person(self, p: Person) -> None:
        # Detach from sibling chains of any parent
        if p.mother is not None:
            p.remove_mother()
        if p.father is not None:
            p.remove_father()
        # Fix children — orphan them
        c = p.child
        while c is not None:
            nxt = (c.paternal_sibling if p.male else c.maternal_sibling)
            if p.male:
                c.father = None
                c.paternal_sibling = None
            else:
                c.mother = None
                c.maternal_sibling = None
            c = nxt
        # Unlink from list
        idx = self.persons.index(p)
        self.persons.pop(idx)
        del self._by_name[p.name]
        # Fix only the single break point in the ``next`` chain.
        if idx > 0:
            self.persons[idx - 1].next = self.persons[idx] if idx < len(self.persons) else None
        p.next = None

    def add_relation(self, parent: Person, child: Person) -> None:
        """Attach ``parent`` as a parent of ``child``; raises on cycle."""
        if child.has_ancestor(parent) and parent is not child:
            # parent is descendant of child => would create cycle
            pass  # (this check is the wrong direction; see next line)
        if parent.has_ancestor(child):
            raise ValueError("Adding this relation would create a cycle")
        if parent.male:
            if child.father is not None and child.father is not parent:
                raise ValueError("Child already has a different father")
            if child.father is parent:
                return
        else:
            if child.mother is not None and child.mother is not parent:
                raise ValueError("Child already has a different mother")
            if child.mother is parent:
                return
        child.add_parent(parent)

    def remove_relation(self, parent: Person, child: Person) -> None:
        if parent.male and child.father is parent:
            child.remove_father()
        elif (not parent.male) and child.mother is parent:
            child.remove_mother()
