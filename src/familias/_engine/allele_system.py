"""Port of the C++ ``allelesystem`` (alsys.cpp) and ``alleledata`` (aldata.cpp).

Each :class:`AlleleSystem` keeps a list of allele frequencies plus 2-D
mutation matrices (``M[i, j]`` = probability of mutating from allele ``i``
to allele ``j``; rows sum to 1). DNA observations are stored per person
as a tuple of two integer indices into the full allele list.

When called for likelihood computation, :meth:`compute_dataprob` produces a
*reduced* allele system that contains only the alleles observed in the
data (plus, when ``simplify_mutation_matrix`` is true and there are at
least two unobserved non-silent alleles, a single "collapsed" allele that
absorbs the rest), together with the matching reduced mutation matrices.
"""
from __future__ import annotations
from typing import Dict, List, Optional, Tuple
import numpy as np


class AlleleSystem:
    def __init__(
        self,
        name: str,
        allele_names: List[str],
        probability: np.ndarray,
        mutation_matrix_female: np.ndarray,
        mutation_matrix_male: np.ndarray,
        simplify_mutation_matrix: bool,
        has_silent_allele: bool = False,
        silent_allele: int = -1,
    ):
        n = len(allele_names)
        assert probability.shape == (n,)
        assert mutation_matrix_female.shape == (n, n)
        assert mutation_matrix_male.shape == (n, n)
        self.name = name
        self.allele_names: List[str] = list(allele_names)
        self.n_alleles: int = n
        self.probability = np.asarray(probability, dtype=float)
        self.mutation_matrix_female = np.asarray(mutation_matrix_female, dtype=float)
        self.mutation_matrix_male = np.asarray(mutation_matrix_male, dtype=float)
        self.simplify_mutation_matrix = bool(simplify_mutation_matrix)
        self.has_silent_allele = bool(has_silent_allele)
        self.silent_allele = int(silent_allele)
        self.kinship: float = 0.0

        # data: dict[person -> (allele1_idx, allele2_idx)] (full-system indices)
        self.data: Dict[object, Tuple[int, int]] = {}

        # Reduced system, computed lazily via compute_dataprob():
        self.recalc_data: bool = True
        self.n_dataalleles: int = 0
        self.index: np.ndarray = np.zeros(0, dtype=int)
        self.dataprobability: np.ndarray = np.zeros(0)
        self.dataprobmatrix_female: np.ndarray = np.zeros((0, 0))
        self.dataprobmatrix_male: np.ndarray = np.zeros((0, 0))

        # Result of last execute()
        self.result: float = 1.0

    # ---- mutation / kinship setters -----------------------------------
    def set_kinship(self, kinship: float) -> None:
        self.kinship = float(kinship)

    # ---- data manipulation --------------------------------------------
    def add_data(self, person, allele1_idx: int, allele2_idx: int) -> None:
        self.data[person] = (int(allele1_idx), int(allele2_idx))
        self.recalc_data = True

    def remove_data(self, person) -> None:
        if person in self.data:
            del self.data[person]
            self.recalc_data = True

    # ---- core: compute reduced (data) system --------------------------
    def compute_dataprob(self) -> None:
        """Build the reduced allele system used by the peeling algorithm.

        Faithful port of :func:`allelesystem::compute_dataprob` in
        ``src/alsys.cpp`` (the post-2014 version using flat mutation matrices).
        """
        n_alleles = self.n_alleles
        exists_in_data = np.zeros(n_alleles, dtype=bool)

        n_extra = 0
        if self.simplify_mutation_matrix:
            if self.data:
                vals = np.fromiter(
                    (a for pair in self.data.values() for a in pair),
                    dtype=int,
                    count=len(self.data) * 2,
                )
                exists_in_data[vals] = True
            if self.has_silent_allele:
                exists_in_data[self.silent_allele] = True
            n_extra = int((~exists_in_data).sum())

        if n_extra > 1:
            obs = np.where(exists_in_data)[0]
            unobs = np.where(~exists_in_data)[0]
            n_data = len(obs) + 1
            index = np.zeros(n_alleles, dtype=int)
            index[obs] = np.arange(1, len(obs) + 1)
            dataprob = np.zeros(n_data)
            dataprob[0] = 1.0 - self.probability[obs].sum()
            dataprob[1:] = self.probability[obs]

            mf = self.mutation_matrix_female
            mm = self.mutation_matrix_male
            Mf = np.zeros((n_data, n_data))
            Mm = np.zeros((n_data, n_data))

            Mf[1:, 1:] = mf[np.ix_(obs, obs)]
            Mf[1:, 0] = mf[np.ix_(obs, unobs)].sum(axis=1)
            Mf[0, 1:] = self.probability[unobs] @ mf[np.ix_(unobs, obs)]
            Mf[0, 0] = self.probability[unobs] @ mf[np.ix_(unobs, unobs)].sum(axis=1)

            Mm[1:, 1:] = mm[np.ix_(obs, obs)]
            Mm[1:, 0] = mm[np.ix_(obs, unobs)].sum(axis=1)
            Mm[0, 1:] = self.probability[unobs] @ mm[np.ix_(unobs, obs)]
            Mm[0, 0] = self.probability[unobs] @ mm[np.ix_(unobs, unobs)].sum(axis=1)

            s = Mf[0].sum()
            if s > 0:
                Mf[0] /= s
            s = Mm[0].sum()
            if s > 0:
                Mm[0] /= s
        else:
            n_data = n_alleles
            index = np.arange(n_alleles)
            dataprob = self.probability.copy()
            Mf = self.mutation_matrix_female.copy()
            Mm = self.mutation_matrix_male.copy()

        self.index = index
        self.n_dataalleles = n_data
        self.dataprobability = dataprob
        self.dataprobmatrix_female = Mf
        self.dataprobmatrix_male = Mm
        self.recalc_data = False
