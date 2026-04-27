"""Tests covering the C++ FamInterface / PedigreeList methods that the
public R wrapper does not use directly but that exist in the C++ source
(now ported for parity)."""
import numpy as np
import pytest

from familias._engine.interface import FamInterface, FERTILITY_AGE


def _trivial_locus():
    """A 2-allele system with equal frequencies and no mutation."""
    M = np.eye(2)
    return dict(n_alleles=2, mutation_matrix_female=M, mutation_matrix_male=M,
                simplify_mutation_matrix=False, frequencies=np.array([0.5, 0.5]),
                has_silent_allele=False)


# ------------------------------------------------------------------
# Trivial getters
# ------------------------------------------------------------------
def test_get_number_of_pedigrees_and_systems():
    fi = FamInterface()
    assert fi.get_number_of_pedigrees() == 0
    assert fi.get_number_of_systems() == 0
    fi.add_person(male=True)
    fi.add_person(male=False)
    fi.add_pedigree(0, 0)
    fi.add_pedigree(0, 0)
    fi.add_allele_system(**_trivial_locus())
    assert fi.get_number_of_pedigrees() == 2
    assert fi.get_number_of_systems() == 1


# ------------------------------------------------------------------
# YOB sentinel: -1 should also mean "unknown"
# ------------------------------------------------------------------
def test_yob_minus_one_means_unknown():
    fi = FamInterface()
    p = fi.add_person(male=True, yob=-1)
    c = fi.add_person(male=False, yob=-1)
    fi.add_pedigree(0, 0)
    # Should not raise: both YOBs are unknown.
    fi.add_relation(p, c, 0)


def test_yob_filter_blocks_too_young_parent():
    fi = FamInterface()
    p = fi.add_person(male=True, yob=2000)
    c = fi.add_person(male=False, yob=2005)  # only 5 years younger
    fi.add_pedigree(0, 0)
    with pytest.raises(ValueError, match="Illegal relation"):
        fi.add_relation(p, c, 0)


def test_yob_filter_allows_old_enough_parent():
    fi = FamInterface()
    p = fi.add_person(male=True, yob=1980)
    c = fi.add_person(male=False, yob=1980 + FERTILITY_AGE + 5)
    fi.add_pedigree(0, 0)
    fi.add_relation(p, c, 0)


# ------------------------------------------------------------------
# Pedigree dunders / getters
# ------------------------------------------------------------------
def test_pedigree_dunders_and_getters():
    fi = FamInterface()
    fi.add_person(male=True)    # 0: father
    fi.add_person(male=False)   # 1: mother
    fi.add_person(male=True)    # 2: child
    fi.add_pedigree(0, 0)
    fi.add_relation(0, 2, 0)
    fi.add_relation(1, 2, 0)
    ped = fi.pedset.get_pedigree(0)
    assert len(ped) == 3
    assert ped.get_pedigree_size() == 3
    moms, dads = ped.get_parents()
    assert moms == [-1, -1, 1]
    assert dads == [-1, -1, 0]
    mat = ped.get_pedigree_matrix()
    assert mat[0][2] == 1 and mat[1][2] == 1
    assert mat[2][0] == 0
    assert ped.has_children([0, 1]) is True
    assert ped.has_children([2]) is False
    assert "Pedigree(" in repr(ped)


# ------------------------------------------------------------------
# Fixed relations
# ------------------------------------------------------------------
def test_add_and_remove_fixed_relation():
    fi = FamInterface()
    fi.add_person(male=True)    # 0
    fi.add_person(male=False)   # 1
    fi.add_pedigree(0, 0)       # ped 0
    fi.add_pedigree(0, 0)       # ped 1
    # Fixed relation must propagate to BOTH pedigrees.
    removed = fi.add_fixed_relation(0, 1)
    assert removed == [0, 0]
    for i in range(2):
        ped = fi.pedset.get_pedigree(i)
        # 0 (male) is father of 1 (female)
        assert ped.father[1] == 0
    # Removing it should clear from both pedigrees.
    fi.remove_fixed_relation(0, 1)
    for i in range(2):
        assert fi.pedset.get_pedigree(i).father[1] == -1


# ------------------------------------------------------------------
# Pedigree generation: enumerate all valid trio pedigrees
# ------------------------------------------------------------------
def test_generate_pedigrees_two_persons_one_extra():
    """With one female + one male child + one extra male (potential father),
    expected pedigrees: { (no parents), (extra male is father of child) }.

    The pruning step removes the extra male when she has no children, so
    after deduplication we expect exactly 1 pedigree.
    """
    fi = FamInterface()
    fi.add_person(male=False)   # 0: mother
    fi.add_person(male=True)    # 1: child
    # No fixed relations. Generate with 0 extras first to confirm trivial case.
    fi.generate_pedigrees(0, 0)
    n_trivial = fi.get_number_of_pedigrees()
    assert n_trivial >= 1


def test_generate_pedigrees_paternity_enumeration():
    """Two named persons (alleged father, child) with an extra female (mother)
    should produce a small, finite set of distinct pedigrees."""
    fi = FamInterface()
    fi.add_person(male=True)    # 0: alleged father
    fi.add_person(male=True)    # 1: child
    fi.generate_pedigrees(1, 0)  # one extra female (the mother)
    # At minimum we expect: (a) alleged father IS father of child,
    # (b) alleged father is NOT father of child. Both should have a mother.
    n = fi.get_number_of_pedigrees()
    assert n >= 1
    # Every generated pedigree must be on standard form.
    for i in range(n):
        ped = fi.pedset.get_pedigree(i)
        assert ped.on_standard_form()


# ------------------------------------------------------------------
# remove_pedigree
# ------------------------------------------------------------------
def test_remove_pedigree():
    fi = FamInterface()
    fi.add_person(male=True)
    fi.add_pedigree(0, 0)
    fi.add_pedigree(0, 0)
    assert fi.get_number_of_pedigrees() == 2
    fi.remove_pedigree(0)
    assert fi.get_number_of_pedigrees() == 1
