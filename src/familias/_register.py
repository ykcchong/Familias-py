"""Common helper used by :mod:`familias.prior` and :mod:`familias.posterior`.

Encapsulates the registration of persons + pedigrees into a shared
:class:`FamInterface` instance — the equivalent of the R code that walks
a list of ``FamiliasPedigree`` objects and issues ``AddPerson`` /
``AddPedigree`` / ``AddRelation`` calls.
"""
from __future__ import annotations
from typing import List, Sequence, Tuple
from ._engine.interface import FamInterface
from .pedigree import FamiliasPedigree


def normalise_pedigrees(pedigrees) -> List[FamiliasPedigree]:
    if isinstance(pedigrees, FamiliasPedigree):
        return [pedigrees]
    pedigrees = list(pedigrees)
    for p in pedigrees:
        if not isinstance(p, FamiliasPedigree):
            raise TypeError("All pedigrees must be FamiliasPedigree instances.")
    if not pedigrees:
        raise ValueError("At least one pedigree required.")
    return pedigrees


def _common_persons(pedigrees: List[FamiliasPedigree]) -> List[str]:
    common = list(pedigrees[0].id)
    for p in pedigrees[1:]:
        ids = set(p.id)
        common = [x for x in common if x in ids]
    return common


def register(
    pedigrees: List[FamiliasPedigree],
    fi: FamInterface,
    persons: Sequence[str],
    require_consistent_sex: bool = True,
) -> None:
    """Register ``persons`` (in order) and all pedigrees in ``fi``.

    ``persons`` is the (ordered) list of common persons that will become
    indices 0..len(persons)-1 inside the engine; remaining pedigree members
    become extra females / males.
    """
    # Determine each common person's sex (use first pedigree they appear in).
    person_sex: List[str] = []
    for name in persons:
        for p in pedigrees:
            if name in p.id:
                idx = p.id.index(name)
                person_sex.append(p.sex[idx])
                break
        else:  # pragma: no cover - guarded by caller
            raise KeyError(name)

    # Cross-check the sex across pedigrees
    if require_consistent_sex:
        for j, name in enumerate(persons):
            for p in pedigrees:
                if name in p.id:
                    idx = p.id.index(name)
                    if p.sex[idx] != person_sex[j]:
                        raise ValueError(
                            f"Person {name!r} has inconsistent sex across pedigrees."
                        )

    for s in person_sex:
        fi.add_person(male=(s != "female"))

    # Add each pedigree
    npersons_common = len(persons)
    for ped in pedigrees:
        n_persons_ped = len(ped)
        # Compute extra-female / extra-male counts and the new-order mapping
        new_order = [0] * n_persons_ped
        n_ex_females = 0
        n_ex_males = 0
        for j in range(n_persons_ped):
            pid = ped.id[j]
            if pid in persons:
                new_order[j] = persons.index(pid)  # 0-based
            elif ped.sex[j] == "female":
                new_order[j] = npersons_common + n_ex_females
                n_ex_females += 1
            else:
                # placeholder; we'll patch below once we know n_ex_females
                new_order[j] = -(n_ex_males + 1)
                n_ex_males += 1
        # Patch males now that n_ex_females is known
        for j in range(n_persons_ped):
            if new_order[j] < 0:
                m_idx = -new_order[j] - 1
                new_order[j] = npersons_common + n_ex_females + m_idx

        index = fi.add_pedigree(n_ex_females, n_ex_males)
        for j in range(n_persons_ped):
            f = ped.findex[j]  # 1-based, 0 = none
            if f > 0:
                fi.add_relation(new_order[f - 1], new_order[j], index)
            m = ped.mindex[j]
            if m > 0:
                fi.add_relation(new_order[m - 1], new_order[j], index)
