"""Port of ``R/FamiliasPedigree.R``.

A :class:`FamiliasPedigree` describes a family graph. Each person has an
``id``, a father id (``dadid``) and a mother id (``momid``) — both may be
``"NA"``/``None`` to indicate "unknown" — and a ``sex`` ("male" or
"female"). The constructor checks that ids are unique, that referenced
parents exist (and have the correct sex), and that the relation graph is
acyclic.
"""
from __future__ import annotations
from typing import List, Optional, Sequence, Union


_NA_LIKE = {None, "NA", "na", "Na", "<NA>"}


def _normalise_parent(value) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, str) and value in _NA_LIKE:
        return None
    return str(value)


class FamiliasPedigree:
    """A pedigree definition (analog of the R ``FamiliasPedigree`` object)."""

    def __init__(
        self,
        id: Sequence[Union[str, int]],
        dadid: Sequence,
        momid: Sequence,
        sex: Sequence[str],
    ):
        n = len(id)
        if not (len(dadid) == n and len(momid) == n and len(sex) == n):
            raise ValueError("id, dadid, momid and sex must have equal length.")
        ids = [str(x) for x in id]
        if len(set(ids)) != n:
            raise ValueError("All ids must be unique.")
        sexes = []
        for s in sex:
            s = str(s).lower()
            if s not in ("male", "female"):
                raise ValueError(f"sex must be 'male' or 'female', got {s!r}")
            sexes.append(s)
        dad = [_normalise_parent(x) for x in dadid]
        mom = [_normalise_parent(x) for x in momid]
        # parents must be valid ids and have the correct sex
        index = {pid: i for i, pid in enumerate(ids)}
        findex = []  # 0 means none, otherwise 1-based index (mirrors R)
        mindex = []
        for i, (d, m) in enumerate(zip(dad, mom)):
            if d is None:
                findex.append(0)
            else:
                if d not in index:
                    raise ValueError(f"Father {d!r} of {ids[i]!r} not in id list.")
                if sexes[index[d]] != "male":
                    raise ValueError(f"Father {d!r} is not male.")
                findex.append(index[d] + 1)
            if m is None:
                mindex.append(0)
            else:
                if m not in index:
                    raise ValueError(f"Mother {m!r} of {ids[i]!r} not in id list.")
                if sexes[index[m]] != "female":
                    raise ValueError(f"Mother {m!r} is not female.")
                mindex.append(index[m] + 1)
        # Acyclicity check
        for i in range(n):
            if _is_ancestor(i, i, findex, mindex, set()):
                raise ValueError("Pedigree is cyclic.")
        self.id: List[str] = ids
        self.dadid: List[Optional[str]] = dad
        self.momid: List[Optional[str]] = mom
        self.sex: List[str] = sexes
        # 1-based parent indices (0 means none) — matches R semantics
        self.findex: List[int] = findex
        self.mindex: List[int] = mindex

    def __len__(self) -> int:
        return len(self.id)

    def __repr__(self) -> str:
        return f"FamiliasPedigree(n={len(self)})"


def _is_ancestor(target: int, person: int, findex, mindex, seen) -> bool:
    """Return True if ``target`` (0-based) is a strict ancestor of ``person``.

    Used to detect cycles: ``person`` is itself a member of ``seen`` initially,
    and we recurse via parents.
    """
    f = findex[person] - 1
    m = mindex[person] - 1
    if f >= 0:
        if f == target and person != target:
            return True
        if f == target and person == target:
            return True
        if f not in seen:
            seen.add(f)
            if _is_ancestor(target, f, findex, mindex, seen):
                return True
    if m >= 0:
        if m == target:
            return True
        if m not in seen:
            seen.add(m)
            if _is_ancestor(target, m, findex, mindex, seen):
                return True
    return False
