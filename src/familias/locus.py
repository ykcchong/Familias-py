"""Port of ``R/FamiliasLocus.R`` — locus / mutation-model construction.

The result mimics the R ``FamiliasLocus`` list with these keys::

    locusname, alleles, femaleMutationType, femaleMutationMatrix,
    maleMutationType, maleMutationMatrix, simpleMutationMatrices,
    Stabilization

where ``alleles`` is a dict ``{name: frequency}``, both mutation matrices
are square :class:`numpy.ndarray` row-stochastic with index order matching
the allele dict.
"""
from __future__ import annotations
from typing import Dict, Optional, Sequence, Union
import warnings
import numpy as np
from scipy.optimize import brentq, minimize


# ---------------------------------------------------------------------------
# Stabilisation helpers (port of stabilize / pKCompute / theOpting* / SgivenMopt*)
# ---------------------------------------------------------------------------

def _fratio(A: np.ndarray, B: np.ndarray) -> float:
    return float(max(np.max(A / B), np.max(B / A)))


def _pK_compute(paj: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Return p satisfying p*(I-D(p)) = K * paj * (I-D(d))."""
    paj = np.asarray(paj, dtype=float).flatten()
    n = paj.size
    w = paj * (1 - d)
    if 2 * w.max() > w.sum():
        raise ValueError("Task impossible with this input")
    ord_ = np.argsort(-w)
    Perm = np.zeros((n, n))
    for i in range(n):
        Perm[i, ord_[i]] = 1
    w_sorted = Perm @ w
    eps = 1e-5
    p = np.zeros(n)
    Kmax = 0.25 / w_sorted[0]

    def ffun(x):
        return 0.5 * np.sum(1 - np.sqrt(1 - 4 * x * w_sorted)) - 1

    def gfun(x):
        return ffun(x) + np.sqrt(1 - 4 * x * w_sorted[0])

    if ffun(Kmax) > 0:
        K = brentq(ffun, 0.0, Kmax)
        p = 0.5 * (1 - np.sqrt(1 - 4 * K * w_sorted))
    else:
        while gfun(eps) < 0:
            eps /= 2
        K = brentq(gfun, eps, Kmax)
        p[0] = 0.5 * (1 + np.sqrt(1 - 4 * K * w_sorted[0]))
        p[1:] = 0.5 * (1 - np.sqrt(1 - 4 * K * w_sorted[1:]))
    return Perm.T @ p


def _SgivenMopt(M, S0, p, *, rm_lower=None):
    n = M.shape[0]
    steg = 1e-2 if rm_lower is None else 1e-3
    S = S0.copy()
    best = _fratio(M, S)
    changed = False
    if np.max(M / S0) >= np.max(S0 / M):
        idx = np.argwhere(np.isclose(M / S0, best))
    else:
        idx = np.argwhere(np.isclose(S0 / M, best))
        steg = -steg
    iis = idx[:, 0]
    jjs = idx[:, 1]
    for i in iis:
        for j in jjs:
            if rm_lower is None and j == i:
                continue
            for i1 in range(n):
                if i1 == i:
                    continue
                if rm_lower is None and i1 == j:
                    continue
                for j1 in range(n):
                    if j1 == j:
                        continue
                    if rm_lower is None and (j1 == i or j1 == i1):
                        continue
                    while True:
                        S_temp = S.copy()
                        S_temp[i, j] = (1 + steg) * S[i, j]
                        S_temp[i, j1] = S[i, j] + S[i, j1] - S_temp[i, j]
                        S_temp[i1, j] = (S[i, j] * p[i] + S[i1, j] * p[i1]
                                         - S_temp[i, j] * p[i]) / p[i1]
                        S_temp[i1, j1] = S[i1, j] + S[i1, j1] - S_temp[i1, j]
                        f1 = _fratio(M[np.ix_([i, i1], [j, j1])],
                                     S_temp[np.ix_([i, i1], [j, j1])])
                        f2 = _fratio(M[np.ix_([i, i1], [j, j1])],
                                     S[np.ix_([i, i1], [j, j1])])
                        ok = S_temp.min() > 0 and (f1 - f2) < 0
                        if rm_lower is not None:
                            ok = ok and np.diag(S_temp).min() > rm_lower
                        if ok:
                            changed = True
                            S = S_temp
                            best = _fratio(M, S)
                        else:
                            break
    return S, changed


def _theOpting(M, S, p, *, rm_lower=None):
    go = True
    while go:
        S, go = _SgivenMopt(M, S, p, rm_lower=rm_lower)
    return S


def _stabilize(M: np.ndarray, pe: np.ndarray, method: str, t: float = 1.0) -> dict:
    n = M.shape[0]
    R = 1 - float(np.sum(np.diag(M) * pe))
    tol = 1e-10
    if np.allclose(M, np.eye(n)):
        return {"stabilized": M, "fratio": 1.0, "mindiag": np.diag(M).min(), "error": ""}
    if (M == 0).any():
        return {"stabilized": M, "fratio": 1.0, "mindiag": np.diag(M).min(),
                "error": "Cannot stabilize non-identity matrices with zero entries."}
    if method == "DP":
        if 2 * np.max(pe * (1 - np.diag(M))) > np.sum(pe * (1 - np.diag(M))):
            return {"stabilized": M, "fratio": 1.0, "mindiag": np.diag(M).min(),
                    "error": "DP stabilization doesn't exist."}
        pnew = _pK_compute(pe, np.diag(M))
        with np.errstate(divide="ignore", invalid="ignore"):
            P0 = ((np.eye(n) - np.diag(np.diag(M)))
                  @ np.diag(np.where(pnew == 1, 0, 1.0 / (1 - pnew)))
                  @ (np.outer(np.ones(n), pnew) - np.diag(pnew))) + np.diag(np.diag(M))
        P = _theOpting(M, P0, pe)
        if np.max(np.abs(pe @ P - pe)) > tol:
            return {"error": "The proposed stabilization doesn't have the desired stationary distribution."}
        if np.max(np.abs(P @ np.ones(n) - np.ones(n))) > tol:
            return {"error": "The proposed stabilization isn't a proper mutation matrix."}
        if P.min() < 0:
            return {"error": "The proposed stabilization has negative elements."}
        return {"stabilized": P, "fratio": _fratio(P, M), "mindiag": np.diag(P).min(), "error": ""}
    if method == "RM":
        if t < R + tol or t > 1:
            return {"error": "MaxStabilizedMutrate parameter out of bounds."}
        # Build constraint matrix C and rhs b (analogous to R formulation)
        m = n * n
        C = np.zeros((2 * n, m))
        for i in range(n):
            C[0:n, n * i:n * (i + 1)] = np.eye(n)
            C[n + i, n * i:n * (i + 1)] = pe
        # Last constraint: sum_i pe[i] * P[i, n-1] for last column == something  -- mimic R
        C[2 * n - 1, :] = 0
        for i in range(n):
            C[2 * n - 1, n * (n - 1) + i] = 0  # cleared by R as well
        b = np.concatenate([np.ones(n), pe[:-1], [1 - R]])
        xM = M.flatten(order="F")
        lb = ((1 - t) * np.eye(n)).flatten(order="F")
        ub = np.ones(m)

        def make_obj(target):
            def obj(x):
                x = np.clip(x, 1e-12, None)
                return float(max(np.sum(x / target), np.sum(target / x)))
            return obj

        cons = [{"type": "eq", "fun": lambda x, C=C, b=b: C @ x - b}]
        bounds = list(zip(lb, ub))
        res1 = minimize(make_obj(xM), xM, bounds=bounds, constraints=cons, method="SLSQP",
                        options={"maxiter": 200, "disp": False})
        x0 = res1.x

        def obj2(x):
            x = np.clip(x, 1e-12, None)
            return float(np.max(np.abs(np.log(x) - np.log(xM))))

        res2 = minimize(obj2, x0, bounds=bounds, constraints=cons, method="SLSQP",
                        options={"maxiter": 200, "disp": False})
        P0 = res2.x.reshape((n, n), order="F")
        P = _theOpting(M, P0, pe, rm_lower=1 - t)
        return {"stabilized": P, "fratio": _fratio(P, M),
                "mindiag": np.diag(P).min(), "error": ""}
    if method == "PM":
        X = M.T - np.eye(n)
        X[-1, :] = 1
        rhs = np.zeros(n)
        rhs[-1] = 1
        v = np.linalg.solve(X, rhs)
        d = R * v / ((1 - np.sum(np.diag(M) * v)) * pe)
        if (d * (1 - np.diag(M)) > t).any():
            return {"stabilized": M, "fratio": 1.0, "mindiag": np.diag(M).min(),
                    "error": "PM stabilization doesn't exist."}
        P = np.diag(d) @ (M - np.eye(n)) + np.eye(n)
        return {"stabilized": P, "fratio": _fratio(P, M),
                "mindiag": np.diag(P).min(), "error": ""}
    raise ValueError("Stabilization method must be 'DP', 'RM' or 'PM'.")


# ---------------------------------------------------------------------------
# Mutation matrix builders
# ---------------------------------------------------------------------------

def _mutation_equal(n_all: int, n_alleles: int, rate: float) -> np.ndarray:
    M = np.eye(n_alleles)
    if rate == 0 or n_all <= 1:
        return M
    off = rate / (n_all - 1)
    for i in range(n_all):
        for j in range(n_all):
            M[i, j] = (1 - rate) if i == j else off
    return M


def _mutation_proportional(freq: np.ndarray, n_alleles: int, rate: float) -> np.ndarray:
    M = np.eye(n_alleles)
    n_all = freq.size
    if rate == 0:
        return M
    sumfreq = freq.sum()
    frq = freq / sumfreq
    alpha = rate / sumfreq / float(np.sum(frq * (1 - frq)))
    for i in range(n_all):
        for j in range(n_all):
            M[i, j] = (1 - alpha + alpha * frq[j]) if i == j else alpha * frq[j]
    return M


def _mutation_stepwise(allele_names: Sequence[str], freq: np.ndarray, n_alleles: int,
                       rate: float, range_: float, rate2: float) -> np.ndarray:
    M = np.eye(n_alleles)
    n_all = freq.size
    if rate == 0 and rate2 == 0:
        return M
    try:
        numfreq = np.array([float(s) for s in allele_names[:n_all]], dtype=float)
    except ValueError as e:
        raise ValueError(
            "The 'Stepwise' mutation model requires all non-silent alleles to have numerical names."
        ) from e
    if not np.allclose(np.round(numfreq, 1), numfreq):
        raise ValueError("Microvariants must be named as a decimal number with one decimal.")
    microgroup = np.round((numfreq - np.round(numfreq)) * 10).astype(int)
    for i in range(n_all):
        microcompats = (microgroup == microgroup[i])
        for j in range(n_all):
            if i == j:
                if microcompats.all():
                    M[i, j] = 1 - rate
                elif microcompats.sum() == 1:
                    M[i, j] = 1 - rate2
                else:
                    M[i, j] = 1 - rate - rate2
            elif microcompats[j]:
                M[i, j] = range_ ** abs(numfreq[i] - numfreq[j])
            else:
                M[i, j] = rate2 / (n_all - microcompats.sum())
        # Re-normalise the microcompat off-diagonal entries to sum to ``rate``
        mc = microcompats.copy()
        mc[i] = False
        if mc.any():
            s = M[i, mc].sum()
            if s > 0:
                M[i, mc] = M[i, mc] / s * rate
    return M


def _build_mutation_matrix(model: str, allele_names, freq: np.ndarray,
                           n_alleles: int, rate: float, range_: float,
                           rate2: float, custom: Optional[np.ndarray]):
    model_l = model.lower()
    if model_l == "custom":
        if custom is None:
            raise ValueError("Custom mutation model requires a mutation matrix.")
        M = np.asarray(custom, dtype=float)
        if M.shape != (n_alleles, n_alleles):
            raise ValueError("Mutation matrix dimension must match number of alleles.")
        if (M < 0).any():
            raise ValueError("Mutation matrix entries must be >= 0.")
        if not np.allclose(M.sum(axis=1), 1.0, atol=1e-6):
            raise ValueError("Rows of mutation matrix must sum to 1.")
        return M, "A 'Custom' specified mutation matrix"
    if model_l == "equal":
        M = _mutation_equal(freq.size, n_alleles, rate)
        if rate == 0:
            return M, "No mutations"
        return M, f"An 'Equal' mutation model with mutation rate {rate}"
    if model_l == "proportional":
        M = _mutation_proportional(freq, n_alleles, rate)
        if rate == 0:
            return M, "No mutations"
        return M, f"A 'Proportional' mutation model with expected mutation rate {rate}"
    if model_l == "stepwise":
        M = _mutation_stepwise(list(allele_names), freq, n_alleles, rate, range_, rate2)
        if rate == 0 and rate2 == 0:
            return M, "No mutations"
        # microgroup all-zero => simple message
        try:
            numfreq = np.array([float(s) for s in allele_names[:freq.size]])
            microgroup = np.round((numfreq - np.round(numfreq)) * 10).astype(int)
        except Exception:
            microgroup = np.array([1])
        if (microgroup == 0).all():
            return M, (f"A 'Stepwise' mutation model with mutation rate {rate} "
                       f"and mutation range {range_}")
        return M, (f"A 'Stepwise' mutation model with mutation rate {rate}, "
                   f"range {range_}, and fractional mutation rate {rate2}")
    raise ValueError(f"Unknown mutation model: {model!r}")


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

class _FamiliasLocusObj(dict):
    """A locus definition (also subscriptable like the R list).

    Attribute access (``loc.alleles``) and item access (``loc['alleles']``)
    both work.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    def __repr__(self) -> str:
        return f"FamiliasLocus({self.get('locusname')!r}, n_alleles={len(self.get('alleles', {}))})"


def _validate_mutation(M: np.ndarray) -> np.ndarray:
    if (M < 0).any() or (M > 1).any():
        raise ValueError("Mutation matrix entries must be in [0, 1].")
    if not np.allclose(M.sum(axis=1), 1.0, atol=1e-3):
        raise ValueError(f"Rows of mutation matrix do not sum to 1: {M.sum(axis=1)}")
    return M


def familias_locus(
    frequencies,
    allelenames: Optional[Sequence[str]] = None,
    name: Optional[str] = None,
    *,
    MutationModel: str = "Stepwise",
    MutationRate: float = 0.0,
    MutationRange: float = 0.5,
    MutationRate2: float = 0.0,
    MutationMatrix: Optional[np.ndarray] = None,
    Stabilization: str = "None",
    MaxStabilizedMutrate: float = 1.0,
    femaleMutationModel: Optional[str] = None,
    femaleMutationRate: Optional[float] = None,
    femaleMutationRange: Optional[float] = None,
    femaleMutationRate2: Optional[float] = None,
    femaleMutationMatrix: Optional[np.ndarray] = None,
    maleMutationModel: Optional[str] = None,
    maleMutationRate: Optional[float] = None,
    maleMutationRange: Optional[float] = None,
    maleMutationRate2: Optional[float] = None,
    maleMutationMatrix: Optional[np.ndarray] = None,
) -> _FamiliasLocusObj:
    """Construct a :class:`FamiliasLocus`.

    ``frequencies`` may be a sequence of floats, a dict ``{name: freq}`` or
    an existing :class:`FamiliasLocus` (in which case mutation / stabilisation
    parameters can be edited).
    """
    # Editing-mode: pass an existing FamiliasLocus
    if isinstance(frequencies, _FamiliasLocusObj):
        x = frequencies
        if name is not None or allelenames is not None:
            raise ValueError("Only mutation parameters can be edited.")
        alleles = dict(x["alleles"])
        frequencies = list(alleles.values())
        allelenames = list(alleles.keys())
        name = x["locusname"]
    elif isinstance(frequencies, dict):
        if allelenames is None:
            allelenames = list(frequencies.keys())
        frequencies = [frequencies[k] for k in allelenames]

    freq_arr = np.asarray(frequencies, dtype=float)
    if (freq_arr <= 0).any():
        raise ValueError("frequencies must all be positive.")
    if round(float(freq_arr.sum()), 6) != 1:
        raise ValueError("Frequencies must sum to 1.")
    if allelenames is None:
        allelenames = [str(i + 1) for i in range(freq_arr.size)]
    elif len(allelenames) != freq_arr.size:
        raise ValueError("len(allelenames) must equal len(frequencies).")
    if len(set(allelenames)) != len(allelenames):
        raise ValueError("Duplicate allele names.")
    if any(s in ("silent", "Silent") for s in allelenames[:-1]):
        raise ValueError("Only the last allele may be silent.")

    n_alleles = freq_arr.size
    has_silent = allelenames[-1] in ("silent", "Silent")
    n_all = n_alleles - int(has_silent)
    freq = freq_arr[:n_all]

    # Per-sex parameter inheritance
    fmm = femaleMutationModel if femaleMutationModel is not None else MutationModel
    mmm = maleMutationModel if maleMutationModel is not None else MutationModel
    fmr = femaleMutationRate if femaleMutationRate is not None else MutationRate
    mmr = maleMutationRate if maleMutationRate is not None else MutationRate
    fmrng = femaleMutationRange if femaleMutationRange is not None else MutationRange
    mmrng = maleMutationRange if maleMutationRange is not None else MutationRange
    fmr2 = femaleMutationRate2 if femaleMutationRate2 is not None else MutationRate2
    mmr2 = maleMutationRate2 if maleMutationRate2 is not None else MutationRate2
    fmmat = femaleMutationMatrix if femaleMutationMatrix is not None else MutationMatrix
    mmmat = maleMutationMatrix if maleMutationMatrix is not None else MutationMatrix

    for v, label in [(fmr, "fmr"), (mmr, "mmr"), (fmr2, "fmr2"), (mmr2, "mmr2")]:
        if v < 0 or v > 1:
            raise ValueError(f"Mutation rate {label}={v} must be in [0, 1].")
    for v, label in [(fmrng, "fmrng"), (mmrng, "mmrng")]:
        if v <= 0:
            raise ValueError(f"Mutation range {label}={v} must be > 0.")

    fmM, fmType = _build_mutation_matrix(fmm, allelenames, freq, n_alleles,
                                         fmr, fmrng, fmr2, fmmat)
    mmM, mmType = _build_mutation_matrix(mmm, allelenames, freq, n_alleles,
                                         mmr, mmrng, mmr2, mmmat)

    # Stabilisation
    Stab = Stabilization.upper()
    if Stab not in ("NONE", "DP", "RM", "PM"):
        raise ValueError("Stabilization must be one of None, DP, RM, PM.")
    if Stab != "NONE":
        zero_silent = (
            has_silent
            and np.all(mmM[:n_all, n_alleles - 1] == 0)
            and np.all(mmM[n_alleles - 1, :n_all] == 0)
            and np.all(fmM[:n_all, n_alleles - 1] == 0)
            and np.all(fmM[n_alleles - 1, :n_all] == 0)
        )
        if zero_silent:
            sub_pe = freq_arr[:n_all] / freq_arr[:n_all].sum()
            res = _stabilize(fmM[:n_all, :n_all], sub_pe, Stab, MaxStabilizedMutrate)
            if res["error"]:
                warnings.warn(f"Female matrix not stabilised: {res['error']}")
            else:
                fmM[:n_all, :n_all] = res["stabilized"]
            res = _stabilize(mmM[:n_all, :n_all], sub_pe, Stab, MaxStabilizedMutrate)
            if res["error"]:
                warnings.warn(f"Male matrix not stabilised: {res['error']}")
            else:
                mmM[:n_all, :n_all] = res["stabilized"]
        else:
            res = _stabilize(fmM, freq_arr, Stab, MaxStabilizedMutrate)
            if res["error"]:
                warnings.warn(f"Female matrix not stabilised: {res['error']}")
            else:
                fmM = res["stabilized"]
            res = _stabilize(mmM, freq_arr, Stab, MaxStabilizedMutrate)
            if res["error"]:
                warnings.warn(f"Male matrix not stabilised: {res['error']}")
            else:
                mmM = res["stabilized"]

    if n_alleles == 1:
        fmM[0, 0] = 1.0
        mmM[0, 0] = 1.0

    simple = True
    for j in range(n_alleles):
        col = np.delete(fmM[:, j], j)
        if not np.allclose(np.round(col - col[0], 6), 0):
            simple = False
            break
    if simple:
        for j in range(n_alleles):
            col = np.delete(mmM[:, j], j)
            if not np.allclose(np.round(col - col[0], 6), 0):
                simple = False
                break

    alleles = {n: float(f) for n, f in zip(allelenames, freq_arr)}
    if name is None:
        name = "locus"
    res = _FamiliasLocusObj(
        locusname=name,
        alleles=alleles,
        femaleMutationType=fmType,
        femaleMutationMatrix=_validate_mutation(fmM),
        maleMutationType=mmType,
        maleMutationMatrix=_validate_mutation(mmM),
        simpleMutationMatrices=simple,
        Stabilization=Stab,
        hasSilentAllele=has_silent,
    )
    return res


# Public name (matches R)
FamiliasLocus = familias_locus
