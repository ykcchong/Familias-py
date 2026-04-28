"""Port of ``FamInterface.cpp`` — the top-level stateful engine.

The original kept a ``static FamInterface* fint`` global; here we use
plain instance state. Methods mirror the C wrappers in ``main.cpp``:

    NewFamilias / __init__
    AddPerson  -> add_person
    AddPedigree -> add_pedigree
    AddRelation -> add_relation
    AddAlleleSystem -> add_allele_system
    AddDNAObservation -> add_dna_observation
    GetProbabilities -> get_probabilities
"""
from __future__ import annotations
from typing import Dict, List, Tuple, Optional
import numpy as np

from .pater import Pater
from .pedigree_list import PedigreeList
from .allele_system import AlleleSystem


# Constant from FamInterface: ages below this are considered children
FERTILITY_AGE: int = 12


class FamInterface:
    def __init__(self) -> None:
        self.pat = Pater()
        self.pedset = PedigreeList()
        self.n_persons: int = 0
        self.male: List[bool] = []
        self.is_child: List[int] = []
        self.yob: List[int] = []  # year of birth (0 means unknown)
        self.internal_person_name: List[str] = []
        # Allele systems are owned by ``self.pat``, but we also keep a
        # parallel list so we can refer to them by index.
        self.systems: List[AlleleSystem] = []
        self._counter: int = 0
        # ``index[ped]`` maps R-side pedigree index -> internal Pedigree index
        self.pedigrees: List[int] = []

    # --- helpers -----------------------------------------------------
    def _new_internal_name(self) -> str:
        self._counter += 1
        return f"_p{self._counter}"

    # --- person ------------------------------------------------------
    def add_person(self, male: bool, yob: int = 0, is_child: Optional[bool] = None) -> int:
        """Register a new person and propagate to all pedigrees + Pater.

        Returns the new person's index.
        """
        idx = self.n_persons
        self.male.append(bool(male))
        self.yob.append(int(yob))
        if is_child is None:
            is_child = 0
        self.is_child.append(int(is_child))
        name = self._new_internal_name()
        self.internal_person_name.append(name)
        self.pat.add_person(name, bool(male))
        self.pedset.add_person(int(bool(male)))
        self.n_persons += 1
        return idx

    # --- pedigree ----------------------------------------------------
    def add_pedigree(self, n_extra_females: int = 0, n_extra_males: int = 0) -> int:
        return self.pedset.add_pedigree(n_extra_females, n_extra_males)

    def add_relation(self, parent_index: int, child_index: int, ped_index: int) -> None:
        # YOB / fertility / is-child checks apply only to "common" persons
        # registered via add_person. Extra (pedigree-local) persons have
        # no YOB/is_child, so skip those checks.
        # Both 0 and -1 are accepted as "unknown" YOB (R passes -1).
        n_common = self.n_persons
        if parent_index < n_common and child_index < n_common:
            yp, yc = self.yob[parent_index], self.yob[child_index]
            if yp not in (0, -1) and yc not in (0, -1):
                if yp >= yc - FERTILITY_AGE + 1:
                    raise ValueError(
                        "Illegal relation: parent's YOB rules them out as a parent."
                    )
        if parent_index < n_common and self.is_child[parent_index]:
            raise ValueError("A person flagged as a child cannot be a parent.")
        self.pedset.add_relation_to_pedigree(ped_index, parent_index, child_index)

    # --- trivial getters (parity with C++ FamInterface) -------------
    def get_number_of_pedigrees(self) -> int:
        return self.pedset.n_pedigrees()

    def get_number_of_systems(self) -> int:
        return len(self.systems)

    # --- fixed relations & removal (parity with C++ FamInterface) ----
    def add_fixed_relation(self, parent_index: int, child_index: int) -> List[int]:
        n = self.n_persons
        if not (0 <= parent_index < n and 0 <= child_index < n):
            raise IndexError("person index out of range")
        yp, yc = self.yob[parent_index], self.yob[child_index]
        if yp not in (0, -1) and yc not in (0, -1) and yp >= yc:
            raise ValueError("Illegal fixed relation: parent must be older than child.")
        if self.is_child[parent_index]:
            raise ValueError("A person flagged as a child cannot be a parent.")
        self.pat.add_parent(self.internal_person_name[parent_index],
                            self.internal_person_name[child_index])
        return self.pedset.add_fixed_relation(parent_index, child_index)

    def remove_fixed_relation(self, parent_index: int, child_index: int) -> None:
        self.pat.remove_possible_parent(self.internal_person_name[parent_index],
                                        self.internal_person_name[child_index])
        self.pedset.remove_fixed_relation(parent_index, child_index)

    def remove_person(self, index: int) -> None:
        if not (0 <= index < self.n_persons):
            raise IndexError("person index out of range")
        self.pat.remove_person(self.internal_person_name[index])
        self.pedset.remove_person(index)
        del self.male[index]
        del self.is_child[index]
        del self.yob[index]
        del self.internal_person_name[index]
        self.n_persons -= 1

    def remove_pedigree(self, index: int) -> None:
        self.pedset.remove_pedigree(index)

    def generate_pedigrees(self, n_extra_females: int, n_extra_males: int) -> None:
        """Enumerate all pedigrees consistent with current fixed relations,
        adding extras and respecting YOB / is-child constraints (port of
        ``FamInterface::GeneratePedigrees``).
        """
        if n_extra_females < 0 or n_extra_males < 0:
            raise ValueError("extras counts must be non-negative")
        n = self.n_persons
        n_total = n + n_extra_females + n_extra_males
        # Build possibleParent[i+j*n_total]: max generations i may precede j.
        poss = [0] * (n_total * n_total)
        for i in range(n_total):
            if i >= n:
                # Extras have no YOB/is_child constraint -> anything goes.
                for j in range(n_total):
                    poss[i + j * n_total] = n_total
            elif self.is_child[i]:
                # A child cannot be anyone's ancestor.
                for j in range(n_total):
                    poss[i + j * n_total] = 0
            elif self.yob[i] in (0, -1):
                for j in range(n_total):
                    poss[i + j * n_total] = n_total
            else:
                for j in range(n_total):
                    if j >= n or self.yob[j] in (0, -1):
                        poss[i + j * n_total] = n_total
                    else:
                        n_g = (self.yob[j] - self.yob[i]) // FERTILITY_AGE
                        poss[i + j * n_total] = n_g if n_g > 0 else 0
        self.pedset.generate_pedigrees(n_extra_females, n_extra_males, poss)

    # --- allele systems ---------------------------------------------
    def add_allele_system(
        self,
        n_alleles: int,
        mutation_matrix_female: np.ndarray,   # (n_alleles, n_alleles) row-stochastic
        mutation_matrix_male: np.ndarray,
        simplify_mutation_matrix: bool,
        frequencies: np.ndarray,
        has_silent_allele: bool = False,
    ) -> int:
        if n_alleles < 1:
            raise ValueError("nAlleles must be >= 1")
        if (frequencies <= 0).any():
            raise ValueError("All allele frequencies must be > 0")
        # Internal names
        sys_name = self._new_internal_name()
        allele_names = [self._new_internal_name() for _ in range(n_alleles)]
        sys = AlleleSystem(
            name=sys_name,
            allele_names=allele_names,
            probability=np.asarray(frequencies, dtype=float),
            mutation_matrix_female=np.asarray(mutation_matrix_female, dtype=float),
            mutation_matrix_male=np.asarray(mutation_matrix_male, dtype=float),
            simplify_mutation_matrix=bool(simplify_mutation_matrix),
            has_silent_allele=bool(has_silent_allele),
            silent_allele=(n_alleles - 1) if has_silent_allele else -1,
        )
        self.pat.add_system(sys)
        self.systems.append(sys)
        return len(self.systems) - 1

    def add_dna_observation(self, person_index: int, system_index: int,
                            allele1_index: int, allele2_index: int) -> None:
        sys = self.systems[system_index]
        person_name = self.internal_person_name[person_index]
        person = self.pat.family.get(person_name)
        if person is None:
            raise KeyError(person_name)
        sys.add_data(person, allele1_index, allele2_index)

    # --- the main entry point --------------------------------------
    def get_probabilities(
        self,
        generations_parameter: float,
        max_generations: int,
        inbreeding_parameter: float,
        promiscuity_parameter: float,
        using_dna_observations: bool,
        kinship: float,
    ) -> Tuple[List[int], List[float], np.ndarray]:
        """Return ``(redundant, probabilities, likelihoods)``.

        ``redundant`` is a list of length ``n_pedigrees`` (1 if removed),
        ``probabilities`` are the normalised priors over kept pedigrees,
        ``likelihoods`` is a 2-D array of shape ``(n_kept_peds, n_systems)``.
        """
        redundant = self.pedset.remove_equivalent_pedigrees()
        probabilities = self.pedset.compute_prior(
            generations_parameter, max_generations,
            inbreeding_parameter, promiscuity_parameter, self.is_child,
        )
        if not probabilities:
            raise RuntimeError("All pedigrees have zero prior probability.")
        if not using_dna_observations:
            n_sys = len(self.systems)
            likelihoods = np.ones((self.pedset.n_pedigrees(), n_sys))
            return redundant, probabilities, likelihoods
        # Posterior: peeling for each (pedigree, system)
        self.pat.set_kinship(kinship)
        names = list(self.internal_person_name)
        like_2d = self.pedset.compute_posterior(
            self.pat,
            make_cutsets=(kinship == 0.0),
            names=names,
        )
        likelihoods = np.array(like_2d, dtype=float)
        return redundant, probabilities, likelihoods
