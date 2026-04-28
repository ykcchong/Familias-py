"""Port of ``Pedigree.cpp`` — a single pedigree (mother/father arrays).

A :class:`Pedigree` stores the parents of each "named" person plus zero or
more "extra" (unnamed) persons. Methods compute pruning, cutsets,
inbreeding, promiscuity, generations and orchestrate the per-pedigree
likelihood computation by driving a :class:`Pater`.
"""
from __future__ import annotations
from typing import List
from .pater import Pater
from .special import get_name_prefix


class Pedigree:
    def __init__(
        self,
        n_named_persons: int,
        n_extra_females: int,
        n_extra_males: int,
        male: List[int],     # length n_named_persons
        parent: List[int],   # length n_total*n_total, parent[j+i*n_total]==1 if j is parent of i
    ):
        self.n_named_persons = n_named_persons
        self.n_total = n_named_persons + n_extra_males + n_extra_females
        self.male: List[int] = list(male) + [0] * n_extra_females + [1] * n_extra_males
        self.father: List[int] = [-1] * self.n_total
        self.mother: List[int] = [-1] * self.n_total
        n = self.n_total
        for i in range(n):
            for j in range(n):
                if parent[j + i * n]:
                    if self.male[j]:
                        self.father[i] = j
                    else:
                        self.mother[i] = j

    # ---- equality / standard form ----------------------------------
    def is_equal_to(self, other: "Pedigree") -> bool:
        if self.n_total != other.n_total or self.n_named_persons != other.n_named_persons:
            return False
        return (self.mother == other.mother and self.father == other.father
                and self.male == other.male)

    def unsafe_equals(self, other: "Pedigree") -> bool:
        return self.mother == other.mother and self.father == other.father

    # ---- dunders + trivial getters (parity with C++ Pedigree) -----
    def __len__(self) -> int:
        return self.n_total

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Pedigree) and self.is_equal_to(other)

    def __repr__(self) -> str:
        return (f"Pedigree(n_named={self.n_named_persons}, "
                f"n_extra_f={self.get_number_of_extra_females()}, "
                f"n_extra_m={self.get_number_of_extra_males()})")

    def get_pedigree_size(self) -> int:
        """C++ ``getPedigreeSize`` — total persons (named + extras)."""
        return self.n_total

    def get_parents(self) -> tuple:
        """C++ ``getParents`` — return ``(mother, father)`` index lists."""
        return (list(self.mother), list(self.father))

    def get_pedigree_matrix(self) -> List[List[int]]:
        """C++ ``getPedigree`` — adjacency matrix where ``m[j][i]==1`` iff
        ``j`` is a parent of ``i``."""
        n = self.n_total
        mat = [[0] * n for _ in range(n)]
        for i in range(n):
            if self.mother[i] >= 0:
                mat[self.mother[i]][i] = 1
            if self.father[i] >= 0:
                mat[self.father[i]][i] = 1
        return mat

    def has_children(self, persons: List[int]) -> bool:
        """C++ ``hasChildren`` — does any person in ``persons`` have a child?"""
        s = set(persons)
        for i in range(self.n_total):
            if self.mother[i] in s or self.father[i] in s:
                return True
        return False

    def get_number_of_extra_females(self) -> int:
        return sum(1 for i in range(self.n_named_persons, self.n_total) if not self.male[i])

    def get_number_of_extra_males(self) -> int:
        return sum(1 for i in range(self.n_named_persons, self.n_total) if self.male[i])

    def on_standard_form(self) -> bool:
        c_f = self.n_named_persons
        c_m = self.n_named_persons + self.get_number_of_extra_females()
        for i in range(self.n_total):
            if self.mother[i] >= c_f and i < c_f:
                if self.mother[i] > c_f:
                    return False
                c_f += 1
            if self.father[i] >= c_m and i < c_m:
                if self.father[i] > c_m:
                    return False
                c_m += 1
        return True

    def change_to_standard_form(self) -> None:
        n = self.n_total
        # permutation[i] gives original index for new position i; invPerm is inverse.
        # The C++ code uses index -1 as a sentinel, achieved via ``permutation++``.
        # We replicate that with a helper that maps -1 to -1 explicitly.
        permutation = list(range(n))
        invPerm = list(range(n))

        def perm(idx: int) -> int:
            return -1 if idx < 0 else permutation[idx]

        def invp(idx: int) -> int:
            return -1 if idx < 0 else invPerm[idx]

        c_f = self.n_named_persons
        c_m = self.n_named_persons + self.get_number_of_extra_females()
        made_change = False
        for i in range(n):
            mp = invp(self.mother[perm(i)])
            if mp >= c_f and i < c_f:
                if mp > c_f:
                    made_change = True
                    a = self.mother[perm(i)]
                    b = invPerm[a]
                    c = c_f
                    d = permutation[c]
                    invPerm[a] = c
                    invPerm[d] = b
                    permutation[c] = a
                    permutation[b] = d
                c_f += 1
            fp = invp(self.father[perm(i)])
            if fp >= c_m and i < c_m:
                if fp > c_m:
                    made_change = True
                    a = self.father[perm(i)]
                    b = invPerm[a]
                    c = c_m
                    d = permutation[c]
                    invPerm[a] = c
                    invPerm[d] = b
                    permutation[c] = a
                    permutation[b] = d
                c_m += 1
        if made_change:
            new_mother = [invp(self.mother[permutation[i]]) for i in range(n)]
            new_father = [invp(self.father[permutation[i]]) for i in range(n)]
            new_male = [self.male[permutation[i]] for i in range(n)]
            self.mother = new_mother
            self.father = new_father
            self.male = new_male

    # ---- inbreeding / promiscuity / generations -------------------
    def is_ancestor(self, ancestor: int, descendant: int) -> bool:
        if ancestor == descendant:
            return True
        stack = [descendant]
        while stack:
            curr = stack.pop()
            f = self.father[curr]
            m = self.mother[curr]
            if f == ancestor or m == ancestor:
                return True
            if f != -1:
                stack.append(f)
            if m != -1:
                stack.append(m)
        return False

    def has_common_ancestor(self, p1: int, p2: int) -> bool:
        ancestors = set()
        stack = [p1]
        while stack:
            curr = stack.pop()
            if curr in ancestors:
                continue
            ancestors.add(curr)
            f = self.father[curr]
            m = self.mother[curr]
            if f != -1:
                stack.append(f)
            if m != -1:
                stack.append(m)
        stack = [p2]
        while stack:
            curr = stack.pop()
            if curr in ancestors:
                return True
            f = self.father[curr]
            m = self.mother[curr]
            if f != -1:
                stack.append(f)
            if m != -1:
                stack.append(m)
        return False

    def compute_inbreeding(self) -> int:
        n = self.n_total
        # Precompute the full ancestor set for every person once, then check
        # pairs via set intersection — O(n²) worst case but avoids redundant
        # traversals from the original O(n) calls to has_common_ancestor().
        ancestor_sets: List[set] = []
        for i in range(n):
            visited: set = set()
            stack = [i]
            while stack:
                curr = stack.pop()
                if curr in visited:
                    continue
                visited.add(curr)
                f = self.father[curr]
                m = self.mother[curr]
                if f != -1:
                    stack.append(f)
                if m != -1:
                    stack.append(m)
            ancestor_sets.append(visited)
        n_inb = 0
        for i in range(n):
            if self.father[i] != -1 and self.mother[i] != -1:
                if ancestor_sets[self.father[i]] & ancestor_sets[self.mother[i]]:
                    n_inb += 1
        return n_inb

    def compute_promiscuity(self) -> int:
        # Group children by parent, then count half-sibling pairs within each
        # group — O(s²) per group (s = group size) rather than O(n²) overall.
        children_by_mother: dict = {}
        children_by_father: dict = {}
        for i in range(self.n_total):
            m = self.mother[i]
            if m >= 0:
                if m not in children_by_mother:
                    children_by_mother[m] = []
                children_by_mother[m].append(i)
            f = self.father[i]
            if f >= 0:
                if f not in children_by_father:
                    children_by_father[f] = []
                children_by_father[f].append(i)
        n_pairs = 0
        # Same mother, different fathers (or both father-less)
        for children in children_by_mother.values():
            for a in range(len(children)):
                for b in range(a + 1, len(children)):
                    ci, cj = children[a], children[b]
                    if self.father[ci] != self.father[cj] or (
                        self.father[ci] == -1 and self.father[cj] == -1
                    ):
                        n_pairs += 1
        # Same father, different mothers (or both mother-less)
        for children in children_by_father.values():
            for a in range(len(children)):
                for b in range(a + 1, len(children)):
                    ci, cj = children[a], children[b]
                    if self.mother[ci] != self.mother[cj] or (
                        self.mother[ci] == -1 and self.mother[cj] == -1
                    ):
                        n_pairs += 1
        return n_pairs

    # ---- helpers ----------------------------------------------------
    def _children_list(self) -> List[List[int]]:
        """Return a list where ``children[i]`` contains all indices of
        children of person ``i``."""
        n = self.n_total
        children: List[List[int]] = [[] for _ in range(n)]
        for j in range(n):
            m = self.mother[j]
            if m >= 0:
                children[m].append(j)
            f = self.father[j]
            if f >= 0:
                children[f].append(j)
        return children

    def _get_max_generations(self, i: int, memo: dict | None = None) -> int:
        if memo is None:
            memo = {}
        if i in memo:
            return memo[i]
        result = 0
        if self.father[i] != -1:
            tmp = self._get_max_generations(self.father[i], memo)
            if tmp > 0 or self.father[i] < self.n_named_persons:
                result = tmp + 1
        if self.mother[i] != -1:
            tmp = self._get_max_generations(self.mother[i], memo)
            if tmp > 0 or self.mother[i] < self.n_named_persons:
                if tmp + 1 > result:
                    result = tmp + 1
        memo[i] = result
        return result

    def compute_generations(self, is_child: List[int]) -> int:
        max_gen = 0
        memo: dict = {}
        for i in range(self.n_named_persons):
            if not is_child[i]:
                g = self._get_max_generations(i, memo)
                if g > max_gen:
                    max_gen = g
        return max_gen

    # ---- pruning / cutsets ---------------------------------------
    def get_pruning(self) -> List[int]:
        n = self.n_total
        result = [0] * n
        children = self._children_list()
        finished = False
        while not finished:
            finished = True
            # Step 1: remove extras with no children
            for i in range(self.n_named_persons, n):
                if result[i] == 0:
                    has_child = any(result[j] == 0 for j in children[i])
                    if not has_child:
                        result[i] = 1
                        finished = False
            # Step 2: remove extras with no parents and at most one child
            for i in range(self.n_named_persons, n):
                if result[i] == 0:
                    if self.father[i] >= 0 and result[self.father[i]] == 0:
                        continue
                    if self.mother[i] >= 0 and result[self.mother[i]] == 0:
                        continue
                    own_children = sum(1 for j in children[i] if result[j] == 0)
                    if own_children > 1:
                        continue
                    result[i] = 1
                    finished = False
        return result

    def prune_and_remove(self) -> None:
        pr = self.get_pruning()
        i = self.n_named_persons
        while i < self.n_total:
            if pr[i]:
                for j in range(self.n_total):
                    if self.mother[j] == i:
                        self.mother[j] = -1
                    if self.father[j] == i:
                        self.father[j] = -1
                    if self.mother[j] > i:
                        self.mother[j] -= 1
                    if self.father[j] > i:
                        self.father[j] -= 1
                # Shift arrays
                del self.mother[i]
                del self.father[i]
                del self.male[i]
                del pr[i]
                self.n_total -= 1
            else:
                i += 1

    def _mark(self, i: int, pruned: List[int], marks: List[int],
              children: List[List[int]]) -> None:
        marks[i] = 1
        for j in children[i]:
            if pruned[j] != 1 and marks[j] == 0:
                self._mark(j, pruned, marks, children)
        if self.father[i] != -1 and pruned[self.father[i]] != 1 and marks[self.father[i]] == 0:
            self._mark(self.father[i], pruned, marks, children)
        if self.mother[i] != -1 and pruned[self.mother[i]] != 1 and marks[self.mother[i]] == 0:
            self._mark(self.mother[i], pruned, marks, children)

    def get_cutsets(self) -> List[int]:
        pruned = self.get_pruning()
        n = self.n_total
        children = self._children_list()
        for i in range(n):
            if pruned[i] != 1:
                marks = [0] * n
                marks[i] = 1
                # find a starting neighbour that is not pruned
                start = n
                for j in children[i]:
                    if pruned[j] != 1:
                        start = j
                        break
                if start == n:
                    if self.father[i] != -1 and pruned[self.father[i]] != 1:
                        start = self.father[i]
                    elif self.mother[i] != -1 and pruned[self.mother[i]] != 1:
                        start = self.mother[i]
                if start == n:
                    continue  # i is unconnected
                self._mark(start, pruned, marks, children)
                # If anything is unmarked, ``i`` is a cutset
                for j in range(n):
                    if marks[j] == 0:
                        pruned[i] = 2
                        break
        return pruned

    # ---- relations / persons (used during pedigree generation) ---
    def add_person(self, male: int) -> None:
        n = self.n_total
        for i in range(n):
            if self.mother[i] >= self.n_named_persons:
                self.mother[i] += 1
            if self.father[i] >= self.n_named_persons:
                self.father[i] += 1
        self.male.insert(self.n_named_persons, male)
        self.mother.insert(self.n_named_persons, -1)
        self.father.insert(self.n_named_persons, -1)
        self.n_total += 1
        self.n_named_persons += 1

    def remove_person(self, index: int) -> None:
        for i in range(self.n_total):
            if self.father[i] > index:
                self.father[i] -= 1
            elif self.father[i] == index:
                self.father[i] = -1
            if self.mother[i] > index:
                self.mother[i] -= 1
            elif self.mother[i] == index:
                self.mother[i] = -1
        del self.male[index]
        del self.father[index]
        del self.mother[index]
        if index < self.n_named_persons:
            self.n_named_persons -= 1
        self.n_total -= 1

    def add_relation(self, parent_index: int, child_index: int) -> None:
        if self.is_ancestor(child_index, parent_index):
            raise ValueError("Adding this relation would create a cycle")
        if self.male[parent_index]:
            if self.father[child_index] >= 0 and self.father[child_index] != parent_index:
                raise ValueError("Child already has a different father")
            self.father[child_index] = parent_index
        else:
            if self.mother[child_index] >= 0 and self.mother[child_index] != parent_index:
                raise ValueError("Child already has a different mother")
            self.mother[child_index] = parent_index

    def remove_relation(self, parent_index: int, child_index: int) -> None:
        if self.male[parent_index]:
            self.father[child_index] = -1
        else:
            self.mother[child_index] = -1

    # ---- the per-pedigree probability orchestration -------------
    def compute_probability(
        self,
        pat: Pater,
        fixed_parent: List[int],
        names: List[str],
        name_prefix: str,
        make_cutsets: bool,
    ) -> List[float]:
        """Replays the C++ ``Pedigree::computeProbability`` orchestration.

        Returns one likelihood per allele system.
        """
        n_named = self.n_named_persons
        n = self.n_total
        # 1. extend 'names' with fresh names for extras, add as persons in pat
        all_names: List[str] = list(names)
        for i in range(n_named, n):
            extra = f"{name_prefix}{i}"
            all_names.append(extra)
            pat.add_person(extra, bool(self.male[i]))
        # 2. add parent relations (those not part of fixedParent)
        added_relations: List[tuple] = []
        for i in range(n):
            if self.mother[i] >= 0 and not (
                self.mother[i] < n_named and i < n_named
                and fixed_parent[self.mother[i] + i * n_named]
            ):
                pat.add_parent(all_names[self.mother[i]], all_names[i])
                added_relations.append((all_names[self.mother[i]], all_names[i]))
            if self.father[i] >= 0 and not (
                self.father[i] < n_named and i < n_named
                and fixed_parent[self.father[i] + i * n_named]
            ):
                pat.add_parent(all_names[self.father[i]], all_names[i])
                added_relations.append((all_names[self.father[i]], all_names[i]))
        # 3. cutsets (only when there is no kinship correction)
        if make_cutsets:
            prune = self.get_cutsets()
            # Ensure the Odds object is built before declaring cutsets.
            pat._ensure_odds()
            for i in range(n):
                if prune[i] == 2:
                    pat.add_person_to_cutset(all_names[i])
                    pat.end_cutset()
        # 4. execute, collect results
        pat.execute()
        result = pat.get_results()
        # 5. tear down cutsets, parent relations, extras (mirroring C++ for state preservation)
        pat.remove_cutsets()
        for parent_name, child_name in added_relations:
            pat.remove_possible_parent(parent_name, child_name)
        for i in range(n_named, n):
            pat.remove_person(all_names[i])
        return result
