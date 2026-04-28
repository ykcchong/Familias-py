"""Port of ``PedigreeList.cpp`` — the public R interface only uses a small
subset of methods (no automated pedigree generation), so we only port:
``add_person``, ``add_pedigree``, ``add_relation``,
``remove_equivalent_pedigrees``, ``compute_prior`` and
``compute_posterior``.
"""
from __future__ import annotations
from typing import Dict, List, Tuple
from .pedigree_core import Pedigree
from .pater import Pater
from .special import mypow, get_name_prefix


class PedigreeList:
    def __init__(self) -> None:
        self.n_named_persons: int = 0
        self.male: List[int] = []
        # Sparse dict: (parent_index, child_index) -> 1 for fixed relations.
        # Replaces the dense n×n flat list to avoid O(n²) reallocation on
        # add_person / remove_person.
        self._fp: Dict[Tuple[int, int], int] = {}
        self.pedigrees: List[Pedigree] = []

    def n_pedigrees(self) -> int:
        return len(self.pedigrees)

    def get_pedigree(self, i: int) -> Pedigree:
        return self.pedigrees[i]

    def add_person(self, male: int) -> None:
        # With a sparse dict there is no matrix to resize — O(1).
        self.male.append(male)
        self.n_named_persons += 1
        for p in self.pedigrees:
            p.add_person(male)

    def add_pedigree(self, n_extra_females: int, n_extra_males: int) -> int:
        """Append a pedigree containing only the fixed relations. Returns its index."""
        n_named = self.n_named_persons
        n_total = n_named + n_extra_females + n_extra_males
        full = [0] * (n_total * n_total)
        for (p, c) in self._fp:
            full[p + c * n_total] = 1
        self.pedigrees.append(
            Pedigree(n_named, n_extra_females, n_extra_males, list(self.male), full)
        )
        return len(self.pedigrees) - 1

    def add_relation_to_pedigree(self, ped_index: int,
                                 parent_index: int, child_index: int) -> None:
        self.pedigrees[ped_index].add_relation(parent_index, child_index)

    # ---- fixed relations (parity with C++ PedigreeList) -------------
    def add_fixed_relation(self, parent_index: int, child_index: int) -> List[int]:
        """Add ``parent->child`` to all pedigrees as a fixed relation.

        Returns a list of length ``n_pedigrees`` where 1 means the pedigree
        was removed because the new relation conflicts.
        """
        n = self.n_named_persons
        if not (0 <= parent_index < n and 0 <= child_index < n):
            raise IndexError("parent/child index out of range")
        self._fp[(parent_index, child_index)] = 1
        kept: List[Pedigree] = []
        removed: List[int] = []
        for p in self.pedigrees:
            try:
                p.add_relation(parent_index, child_index)
                kept.append(p)
                removed.append(0)
            except ValueError:
                removed.append(1)
        self.pedigrees = kept
        return removed

    def remove_fixed_relation(self, parent_index: int, child_index: int) -> None:
        n = self.n_named_persons
        if not (0 <= parent_index < n and 0 <= child_index < n):
            raise IndexError("parent/child index out of range")
        self._fp.pop((parent_index, child_index), None)
        for p in self.pedigrees:
            p.remove_relation(parent_index, child_index)

    def is_fixed_parent(self, parent_index: int, child_index: int) -> bool:
        return (parent_index, child_index) in self._fp

    # ---- removal (parity with C++ PedigreeList) ----------------------
    def remove_person(self, index: int) -> None:
        n_old = self.n_named_persons
        if not (0 <= index < n_old):
            raise IndexError("person index out of range")
        # Rebuild sparse dict: drop entries touching `index`, renumber the rest.
        new_fp: Dict[Tuple[int, int], int] = {}
        for (p, c) in self._fp:
            if p == index or c == index:
                continue
            np_ = p - (1 if p > index else 0)
            nc = c - (1 if c > index else 0)
            new_fp[(np_, nc)] = 1
        self._fp = new_fp
        del self.male[index]
        self.n_named_persons -= 1
        for p in self.pedigrees:
            p.remove_person(index)

    def remove_pedigree(self, index: int) -> None:
        del self.pedigrees[index]

    def add_pedigree_object(self, p: Pedigree) -> int:
        """Append an already-built Pedigree (C++ ``addPedigree(Pedigree*)``)."""
        self.pedigrees.append(p)
        return len(self.pedigrees) - 1

    def remove_equivalent_pedigrees(self) -> List[int]:
        """Mark equivalent pedigrees as redundant and remove them.

        Returns a list of indicators of length ``n_pedigrees`` (1 = removed).
        Uses a hash set for O(1) duplicate lookup instead of the O(n²)
        pairwise comparison used in the original C++ code.
        """
        kept: List[Pedigree] = []
        removed: List[int] = []
        seen: set = set()
        for p in self.pedigrees:
            p.prune_and_remove()
            p.change_to_standard_form()
            key = (p.n_named_persons, tuple(p.mother), tuple(p.father), tuple(p.male))
            if key in seen:
                removed.append(1)
            else:
                seen.add(key)
                kept.append(p)
                removed.append(0)
        self.pedigrees = kept
        return removed

    def compute_prior(self,
                      generations_parameter: float,
                      max_generations: int,
                      inbreeding_parameter: float,
                      promiscuity_parameter: float,
                      is_child: List[int]) -> List[float]:
        """Return a list of normalised prior probabilities (one per pedigree).

        Returns an empty list if the total weight is zero.
        """
        probs: List[float] = []
        s = 0.0
        for p in self.pedigrees:
            n_gen = p.compute_generations(is_child)
            if max_generations != -1 and n_gen > max_generations:
                w = 0.0
            else:
                w = 1.0
                if generations_parameter != 1:
                    w *= mypow(generations_parameter, n_gen)
                if inbreeding_parameter != 1:
                    w *= mypow(inbreeding_parameter, p.compute_inbreeding())
                if promiscuity_parameter != 1:
                    w *= mypow(promiscuity_parameter, p.compute_promiscuity())
            probs.append(w)
            s += w
        if s == 0.0:
            return []
        return [w / s for w in probs]

    def compute_posterior(self, pat: Pater, make_cutsets: bool,
                          names: List[str]) -> List[List[float]]:
        """Return a 2-D list ``likelihoods[ped][system]``."""
        prefix = get_name_prefix(names)
        # Build a flat fixed-parent list once (used by each Pedigree.compute_probability).
        n = self.n_named_persons
        flat_fp = [0] * (n * n)
        for (p, c) in self._fp:
            flat_fp[p + c * n] = 1
        out: List[List[float]] = []
        for ped in self.pedigrees:
            res = ped.compute_probability(pat, flat_fp, names, prefix, make_cutsets)
            out.append(res)
        return out

    # =================================================================
    # Pedigree enumeration (port of PedigreeList::generatePedigrees and
    # the recursive helpers ``generateParentsForPerson`` /
    # ``generateFatherForPerson``).
    # =================================================================
    def generate_pedigrees(self, n_extra_females: int, n_extra_males: int,
                           possible_parent: List[int]) -> None:
        """Append all valid pedigrees that can be built with the given
        number of extra persons, subject to ``possible_parent`` constraints.

        ``possible_parent`` is a ``n_total*n_total`` flat matrix where
        ``possible_parent[i+j*n_total]`` is the maximum number of generations
        ``i`` may be ahead of ``j`` (0 means ``i`` cannot be an ancestor of
        ``j``; ``n_total`` means unrestricted).
        Pedigrees equivalent to one already in the list are skipped.
        """
        n_named = self.n_named_persons
        n_total = n_named + n_extra_females + n_extra_males
        if len(possible_parent) != n_total * n_total:
            raise ValueError(
                f"possible_parent must have length {n_total*n_total}"
            )
        parent = [0] * (n_total * n_total)
        # Seed with fixed relations (only between named persons).
        for (p, c) in self._fp:
            parent[p + c * n_total] = 1
        self._gen_parents_for_person(
            0, parent, n_total,
            n_extra_females, 0,
            n_extra_males, 0,
            possible_parent,
        )

    # ---- internal helpers --------------------------------------------
    def _is_ancestor(self, child: int, parent: int,
                    is_parent: List[int], n_total: int) -> bool:
        if child == parent:
            return True
        for k in range(n_total):
            if is_parent[k + parent * n_total] and self._is_ancestor(
                    child, k, is_parent, n_total):
                return True
        return False

    def _check_parents(self, j: int, i: int, n_gen: int, n_total: int,
                       parent: List[int], poss: List[int]) -> bool:
        for k in range(n_total):
            if parent[k + j * n_total]:
                if poss[k + i * n_total] < n_gen:
                    return False
                if not self._check_parents(k, i, n_gen + 1, n_total, parent, poss):
                    return False
        return True

    def _check_children(self, j: int, i: int, n_gen: int, n_total: int,
                        parent: List[int], poss: List[int]) -> bool:
        for k in range(n_total):
            if parent[i + k * n_total]:
                if poss[j + k * n_total] < n_gen:
                    return False
                if not self._check_children(j, k, n_gen + 1, n_total, parent, poss):
                    return False
        return True

    def _is_possible_parent(self, j: int, i: int, n_total: int,
                            parent: List[int], poss: List[int]) -> bool:
        if poss[j + i * n_total] < 1:
            return False
        if not self._check_parents(j, i, 2, n_total, parent, poss):
            return False
        if not self._check_children(j, i, 2, n_total, parent, poss):
            return False
        return not self._is_ancestor(i, j, parent, n_total)

    def _emit(self, parent: List[int], n_total: int,
              n_extra_females: int, n_extra_males: int) -> None:
        p = Pedigree(self.n_named_persons, n_extra_females, n_extra_males,
                     list(self.male), parent)
        p.prune_and_remove()
        p.change_to_standard_form()
        for q in self.pedigrees:
            if p.is_equal_to(q):
                return
        self.pedigrees.append(p)

    def _gen_parents_for_person(
        self, i: int, parent: List[int], n_total: int,
        n_extra_females: int, n_extra_females_used: int,
        n_extra_males: int, n_extra_males_used: int,
        poss: List[int],
    ) -> None:
        if i == n_total:
            self._emit(parent, n_total, n_extra_females, n_extra_males)
            return
        n_named = self.n_named_persons

        # Try with no (new) mother:
        self._gen_father_for_person(
            i, parent, n_total,
            n_extra_females, n_extra_females_used,
            n_extra_males, n_extra_males_used,
            poss,
        )

        # If a fixed mother is already set for a named person, stop.
        if i < n_named:
            for j in range(n_named):
                if parent[j + i * n_total] and not self.male[j]:
                    return

        # Try with a previously mentioned mother:
        for j in range(n_named + n_extra_females_used):
            if (j >= n_named or not self.male[j]) and \
               self._is_possible_parent(j, i, n_total, parent, poss):
                parent[j + i * n_total] = 1
                self._gen_father_for_person(
                    i, parent, n_total,
                    n_extra_females, n_extra_females_used,
                    n_extra_males, n_extra_males_used,
                    poss,
                )
                parent[j + i * n_total] = 0

        # Try with a new "extra" mother:
        if n_extra_females_used < n_extra_females:
            j = n_named + n_extra_females_used
            if self._is_possible_parent(j, i, n_total, parent, poss):
                parent[j + i * n_total] = 1
                self._gen_father_for_person(
                    i, parent, n_total,
                    n_extra_females, n_extra_females_used + 1,
                    n_extra_males, n_extra_males_used,
                    poss,
                )
                parent[j + i * n_total] = 0

    def _gen_father_for_person(
        self, i: int, parent: List[int], n_total: int,
        n_extra_females: int, n_extra_females_used: int,
        n_extra_males: int, n_extra_males_used: int,
        poss: List[int],
    ) -> None:
        n_named = self.n_named_persons
        # Try with no (new) father:
        self._gen_parents_for_person(
            i + 1, parent, n_total,
            n_extra_females, n_extra_females_used,
            n_extra_males, n_extra_males_used,
            poss,
        )
        if i < n_named:
            for j in range(n_named):
                if parent[j + i * n_total] and self.male[j]:
                    return

        # Try a previously mentioned father:
        upper = n_named + n_extra_females + n_extra_males_used
        for j in range(upper):
            if (j < n_named and self.male[j]) or j >= n_named + n_extra_females:
                if self._is_possible_parent(j, i, n_total, parent, poss):
                    parent[j + i * n_total] = 1
                    self._gen_parents_for_person(
                        i + 1, parent, n_total,
                        n_extra_females, n_extra_females_used,
                        n_extra_males, n_extra_males_used,
                        poss,
                    )
                    parent[j + i * n_total] = 0

        # Try a new "extra" father:
        if n_extra_males_used < n_extra_males:
            j = n_named + n_extra_females + n_extra_males_used
            if self._is_possible_parent(j, i, n_total, parent, poss):
                parent[j + i * n_total] = 1
                self._gen_parents_for_person(
                    i + 1, parent, n_total,
                    n_extra_females, n_extra_females_used,
                    n_extra_males, n_extra_males_used + 1,
                    poss,
                )
                parent[j + i * n_total] = 0
