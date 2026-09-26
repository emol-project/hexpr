"""hexpr -- Python ctypes binding for the hexpr C library.

Quick start::

    import hexpr

    hexpr.init()   # idempotent; safe to call multiple times

    drt = hexpr.Drt(kind='soci', nve=4, spin_x2=0, ncore=0, nval=4, next_=0)
    print(drt.norb, drt.nnodes, drt.ncsfs)

    csfs = drt.csfs()          # shape (ncsfs, norb), uint8 (0=E,1=U,2=D,3=F)
    expr = hexpr.Expression(drt)
    ij, pqrs, coef = expr.terms_full()   # 0-based indices, float64 coef
"""
from __future__ import annotations

import ctypes
import sys
from typing import Optional, Tuple

import numpy as np

from ._ctypes import (
    lib,
    DrtHandle, ExprHandle,
    HEXPR_DRT_SOCI, HEXPR_DRT_FOCI, HEXPR_DRT_CSFSET,
    HEXPR_EXPR_FULL, HEXPR_EXPR_ONE,
)
from ._errors import (
    HexprError, InvalidArg, OutOfMemory, RangeError, NotBuilt, InternalError,
    check_status,
)
from ._errors import IOError as HexprIOError  # avoid shadowing builtin in this scope

__all__ = [
    "init", "shutdown", "version",
    "HexprError", "InvalidArg", "OutOfMemory", "IOError", "RangeError",
    "NotBuilt", "InternalError",
    "Drt", "Expression",
    "step_parse",
    "csf_pair_eval", "csf_block_eval", "csf_subspace_eval",
]

# Re-export IOError under the spec-required name in this module's namespace.
IOError = HexprIOError  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# C stdout FILE* for the *_show* wrappers
# ---------------------------------------------------------------------------

_c_stdout_ptr: Optional[ctypes.c_void_p] = None

# libc exports the standard stdout FILE* as ``stdout`` on Linux but as
# ``__stdoutp`` on macOS.
_STDOUT_SYM = "__stdoutp" if sys.platform == "darwin" else "stdout"


def _c_stdout() -> ctypes.c_void_p:
    """Return the C library's ``stdout`` FILE* as a ctypes pointer.

    Resolved lazily from the global symbol table (``CDLL(None)``) and cached.
    """
    global _c_stdout_ptr
    if _c_stdout_ptr is None:
        _c_stdout_ptr = ctypes.c_void_p.in_dll(ctypes.CDLL(None), _STDOUT_SYM)
    return _c_stdout_ptr


# ---------------------------------------------------------------------------
# Lifecycle
# ---------------------------------------------------------------------------

def init() -> None:
    """Initialize the hexpr library (idempotent)."""
    check_status(lib.hexpr_init(), lib)


def shutdown() -> None:
    """Release library resources (optional; mainly for leak-audit tools)."""
    lib.hexpr_shutdown()


def version() -> str:
    """Return the library version string, e.g. '1.0.0'."""
    raw = lib.hexpr_version_string()
    return raw.decode("utf-8") if raw else ""


# ---------------------------------------------------------------------------
# DRT
# ---------------------------------------------------------------------------

def _kind_to_int(kind) -> int:
    if kind in ("soci", HEXPR_DRT_SOCI):
        return HEXPR_DRT_SOCI
    if kind in ("foci", HEXPR_DRT_FOCI):
        return HEXPR_DRT_FOCI
    if kind in ("csfset", HEXPR_DRT_CSFSET):
        return HEXPR_DRT_CSFSET
    raise ValueError(f"Unknown DRT kind {kind!r}; expected 'soci', 'foci', or 'csfset'")


def _kind_to_str(kind_int: int) -> str:
    if kind_int == HEXPR_DRT_SOCI:   return "soci"
    if kind_int == HEXPR_DRT_FOCI:   return "foci"
    if kind_int == HEXPR_DRT_CSFSET: return "csfset"
    return f"unknown({kind_int})"


def _which_to_int(which: str) -> int:
    w = which.lower()
    if w == "full":
        return HEXPR_EXPR_FULL
    if w == "one":
        return HEXPR_EXPR_ONE
    raise ValueError(f"which must be 'full' or 'one', got {which!r}")


def _as_int32_1d(arr, name: str) -> np.ndarray:
    a = np.ascontiguousarray(arr, dtype=np.int32)
    if a.ndim != 1:
        raise ValueError(f"{name} must be 1-D, got {a.ndim}-D")
    return a


def _ptr32(arr: np.ndarray):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))


def _ptri(arr: np.ndarray):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


def _ptru64(arr: np.ndarray):
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))


class Drt:
    """Distinct Row Table built from hexpr_drt_create.

    Parameters
    ----------
    kind : 'soci' | 'foci' | int
        CI type (HEXPR_DRT_SOCI=0 or HEXPR_DRT_FOCI=1).
    nve : int
        Number of valence electrons.
    spin_x2 : int
        2 * S (0 for singlet, 2 for triplet, ...).
    ncore : int
        Number of core (doubly-occupied) MOs.
    nval : int
        Number of valence (active) MOs.
    next_ : int
        Number of external (virtual) MOs.
        Named ``next_`` to avoid collision with the Python builtin ``next``.
    """

    def __init__(self, kind, nve: int, spin_x2: int,
                 ncore: int, nval: int, next_: int) -> None:
        kind_int = _kind_to_int(kind)
        h = lib.hexpr_drt_create(kind_int, nve, spin_x2, ncore, nval, next_)
        if not h:
            raw = lib.hexpr_last_error()
            msg = raw.decode("utf-8", errors="replace") if raw else "hexpr_drt_create returned NULL"
            raise HexprError(msg)
        self._handle: DrtHandle = h

    # --- scalar properties ---

    @property
    def norb(self) -> int:
        return lib.hexpr_drt_norb(self._handle)

    @property
    def nnodes(self) -> int:
        return lib.hexpr_drt_nnodes(self._handle)

    @property
    def ncsfs(self) -> int:
        return lib.hexpr_drt_ncsfs(self._handle)

    @property
    def walk_ncsfs(self) -> int:
        """Canonical GUGA-walk CSF count.

        For a CSFSET sub-DRT this is the internal canonical count, which may
        be a superset (M >= N) of the input ``ncsfs`` (or a dedup, M <= N);
        for SOCI/FOCI it equals ``ncsfs``. See :meth:`walk_csfs`.
        """
        return lib.hexpr_drt_walk_ncsfs(self._handle)

    @property
    def kind(self) -> str:
        return _kind_to_str(lib.hexpr_drt_kind(self._handle))

    @property
    def ncore(self) -> int:
        return lib.hexpr_drt_ncore(self._handle)

    @property
    def nval(self) -> int:
        return lib.hexpr_drt_nval(self._handle)

    @property
    def next_(self) -> int:
        return lib.hexpr_drt_next(self._handle)

    @property
    def nve_state(self) -> int:
        return lib.hexpr_drt_nve_state(self._handle)

    @property
    def spin_x2(self) -> int:
        return lib.hexpr_drt_spin_x2(self._handle)

    # --- bulk data ---

    def csfs(self) -> np.ndarray:
        """Return CSF step codes as shape (ncsfs, norb) uint8 array.

        Values: 0=E, 1=U, 2=D, 3=F  (matching HEXPR_STEP_* constants).

        Ordering: for a CSFSET sub-DRT from ``Drt.from_csfset`` (or the
        high-level ``csf_*_eval``), ``csfs()`` / ``ncsfs`` / ``csf_index``
        / ``csf_steps`` speak the **input set in input order**: this returns
        exactly the CSFs you supplied, in the order you supplied them, and
        ``ncsfs`` is that input count. The ``bra_pos``/``ket_pos`` from
        ``block_eval``/``subspace_eval`` index that same input array,
        consistent with ``csfs()``. For SOCI/FOCI there is no input set, so
        these are the parametric CSFs.

        The internal DRT may expand to a canonical **superset** (more rows
        than the input, when input CSFs share nodes) or dedup. That
        canonical view (GUGA order, used internally for evaluation) is
        available via :meth:`walk_csfs` / :attr:`walk_ncsfs`.
        """
        ncsfs = self.ncsfs
        norb  = self.norb
        n = ncsfs * norb
        if n == 0:
            return np.empty((ncsfs, norb), dtype=np.uint8)
        ptr = lib.hexpr_drt_csfs(self._handle)
        addr = ctypes.cast(ptr, ctypes.c_void_p).value
        raw = ctypes.string_at(addr, n)
        return np.frombuffer(raw, dtype=np.uint8).reshape(ncsfs, norb).copy()

    def walk_csfs(self) -> np.ndarray:
        """Return the canonical GUGA-walk CSF list, shape (walk_ncsfs, norb).

        For a CSFSET sub-DRT this is the internal canonical view, which may
        be a **superset** of the input ``csfs()`` (``walk_ncsfs >= ncsfs``,
        in canonical GUGA order e<u<d<f) or a dedup; it is the form used
        internally for evaluation. For SOCI/FOCI it equals ``csfs()``.

        Values: 0=E, 1=U, 2=D, 3=F  (matching HEXPR_STEP_* constants).
        """
        ncsfs = self.walk_ncsfs
        norb  = self.norb
        n = ncsfs * norb
        if n == 0:
            return np.empty((ncsfs, norb), dtype=np.uint8)
        ptr = lib.hexpr_drt_walk_csfs(self._handle)
        addr = ctypes.cast(ptr, ctypes.c_void_p).value
        raw = ctypes.string_at(addr, n)
        return np.frombuffer(raw, dtype=np.uint8).reshape(ncsfs, norb).copy()

    def node_info(self, node_index: int) -> Tuple[int, int, int]:
        """Return (orb, nve, is) for the given node index."""
        out_orb = ctypes.c_int()
        out_nve = ctypes.c_int()
        out_is  = ctypes.c_int()
        status = lib.hexpr_drt_node_info(
            self._handle, node_index,
            ctypes.byref(out_orb), ctypes.byref(out_nve), ctypes.byref(out_is),
        )
        check_status(status, lib)
        return (out_orb.value, out_nve.value, out_is.value)

    @classmethod
    def from_csfset(cls, csfs, norb: int) -> "Drt":
        """Build a sub-DRT from an explicit CSF set.

        Parameters
        ----------
        csfs : array-like
            CSF step codes, shape (n_csfs, norb) or flat (n_csfs * norb,).
            Each element is 0..3 (HEXPR_STEP_E/U/D/F).
        norb : int
            Number of orbitals (length of each CSF).

        Returns
        -------
        Drt
            A new Drt instance with kind 'csfset'.

        Notes
        -----
        All CSFs must share a common top node. Invalid step values or
        topology errors cause HexprError.

        Ordering: the resulting sub-DRT's ``csfs()`` / ``ncsfs`` /
        ``csf_index`` / ``csf_steps`` speak the **input set in input
        order** -- ``csfs()`` returns exactly these rows, in this order,
        and ``ncsfs`` equals ``n_csfs``. The ``bra_pos``/``ket_pos`` from
        ``block_eval``/``subspace_eval`` index this same input array. The
        internal DRT may expand to a canonical **superset** (more rows,
        when input CSFs share nodes) or dedup; that canonical view is
        available via :meth:`walk_csfs` / :attr:`walk_ncsfs`.
        """
        arr = np.ascontiguousarray(csfs, dtype=np.uint8)
        if arr.ndim == 2:
            if arr.shape[1] != norb:
                raise ValueError(
                    f"csfs.shape[1]={arr.shape[1]} != norb={norb}")
            n_csfs = arr.shape[0]
            arr = arr.ravel()
        elif arr.ndim == 1:
            if arr.size % norb != 0:
                raise ValueError(
                    f"csfs.size={arr.size} is not a multiple of norb={norb}")
            n_csfs = arr.size // norb
        else:
            raise ValueError(f"csfs must be 1-D or 2-D, got {arr.ndim}-D")

        h = lib.hexpr_drt_from_csfset(
            int(norb),
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
            int(n_csfs))
        if not h:
            raw = lib.hexpr_last_error()
            msg = raw.decode("utf-8", errors="replace") if raw \
                else "hexpr_drt_from_csfset returned NULL"
            raise HexprError(msg)

        inst = cls.__new__(cls)
        inst._handle = h
        return inst

    def csf_index(self, steps) -> int:
        """Return the CSF row index (0-based) for the given step-code sequence."""
        arr = np.ascontiguousarray(steps, dtype=np.uint8)
        if arr.ndim != 1:
            raise ValueError(f"steps must be 1-D, got {arr.ndim}-D")
        if len(arr) != self.norb:
            raise ValueError(f"steps length {len(arr)} != norb {self.norb}")
        out = ctypes.c_int()
        check_status(lib.hexpr_csf_index(
            self._handle,
            arr.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
            ctypes.byref(out),
        ), lib)
        return out.value

    def csf_steps(self, index: int) -> np.ndarray:
        """Return the step codes for CSF row ``index`` as a uint8 array.

        The result has length ``norb`` with values 0..3 (E/U/D/F) and
        round-trips back into :meth:`csf_index`.
        """
        norb = self.norb
        out = np.empty(norb, dtype=np.uint8)
        check_status(lib.hexpr_csf_steps(
            self._handle, int(index),
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        ), lib)
        return out

    # --- CSF -> determinant expansion ---

    def csf_ndets(self, index: int, ms_x2: int) -> int:
        """Number of determinants in the expansion of CSF ``index`` at ``ms_x2``.

        Exact, not an upper bound: it is C(n_open, n_alpha_open), and
        follows from the step vector and ``ms_x2`` alone.

        ``ms_x2`` is 2*Ms. Any |Ms| <= S with the same parity as
        :attr:`spin_x2` is accepted, negative values included; Ms is not
        assumed to equal S. ``index`` follows the same convention as
        :meth:`csf_steps` -- the input position for a CSFSET sub-DRT, the
        canonical row for SOCI/FOCI.
        """
        out = ctypes.c_int()
        check_status(lib.hexpr_csf_ndets(
            self._handle, int(index), int(ms_x2), ctypes.byref(out),
        ), lib)
        return out.value

    def csf_to_dets(self, index: int, ms_x2: int
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Expand CSF ``index`` into Slater determinants at ``ms_x2``.

        Returns ``(alpha, beta, coef)``. ``alpha`` and ``beta`` are uint64
        arrays of shape ``(ndets, nwords)`` with
        ``nwords = (norb + 63) // 64``; ``coef`` is float64 of length
        ``ndets``. Bit ``p`` of word ``w`` is spatial orbital ``64*w + p``,
        so orbital ``p`` is occupied by an alpha electron in determinant
        ``k`` iff ``(alpha[k, p >> 6] >> (p & 63)) & 1``.

        The masks stay 2-D even when ``nwords == 1``; the axis is not
        collapsed, because that would change the caller's indexing at
        norb > 64.

        ``index`` follows the same convention as :meth:`csf_steps` -- the
        input position, in ``range(ncsfs)``. So do :meth:`pair_eval`,
        :meth:`block_eval` and :meth:`subspace_eval`: there is one index
        space across the whole API.

        Conventions
        -----------
        These are fixed by the C API and derived in c/src/csf/csf2det.c.

        * Spin coupling: genealogical, in the Eq. (5) order -- the new
          one-orbital state couples on the LEFT of the accumulated partial
          CSF. A CSF expressed in another coupling scheme will not
          reproduce these coefficients.
        * Clebsch-Gordan phase: Condon-Shortley.
        * Operator order: all alpha in ascending spatial orbital, then all
          beta in ascending spatial orbital. An interleaved-by-orbital
          convention (used by some qubit mappings) differs from this by a
          determinant-dependent permutation sign; csf2det.c gives the rule.
        * Determinant order: lexicographically ascending in the set of
          open-shell positions (walk order, 0-based) carrying alpha.
        * Normalization: ``sum(coef ** 2) == 1`` for every CSF at every Ms.
        * Zero coefficients are emitted, not dropped, so the count is exact
          and no coefficient threshold enters this API.

        The OVERALL sign of a CSF is a convention, not a measured quantity:
        it is whatever the coupling order and the operator order produce.
        Relative signs -- within one CSF, and between CSFs of one DRT --
        are meaningful and are verified against PySCF by
        validation/csf2det_vs_pyscf.py. The overall sign is not something
        to assert against.
        """
        ndets  = self.csf_ndets(index, ms_x2)
        nwords = (self.norb + 63) // 64
        alpha = np.zeros((ndets, nwords), dtype=np.uint64)
        beta  = np.zeros((ndets, nwords), dtype=np.uint64)
        coef  = np.zeros(ndets, dtype=np.float64)
        out_count = ctypes.c_int()
        check_status(lib.hexpr_csf_to_dets(
            self._handle, int(index), int(ms_x2),
            ctypes.byref(out_count),
            _ptru64(alpha), _ptru64(beta),
            coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ), lib)
        n = out_count.value
        if n != ndets:
            raise HexprError(
                f"hexpr_csf_to_dets wrote {n} determinants but "
                f"hexpr_csf_ndets promised {ndets}")
        return alpha, beta, coef

    def show(self, file=None) -> None:
        """Dump a human-readable DRT summary to stdout.

        ``file`` must be ``None`` (stdout) in v1.0.0; arbitrary writable
        files are not yet supported.
        """
        if file is not None:
            raise NotImplementedError(
                "Drt.show supports only stdout (file=None) in v1.0.0")
        check_status(lib.hexpr_drt_show(self._handle, _c_stdout()), lib)

    # --- per-pair evaluation ---

    def pair_count(self, bra: int, ket: int, which: str = "full") -> int:
        """Number of expression terms for the (bra, ket) CSF pair."""
        out = ctypes.c_int()
        check_status(lib.hexpr_pair_count(
            self._handle, int(bra), int(ket), _which_to_int(which),
            ctypes.byref(out),
        ), lib)
        return out.value

    def pair_max_terms(self, which: str = "full") -> int:
        """Upper bound on the term count for any single pair in this DRT."""
        out = ctypes.c_int()
        check_status(lib.hexpr_pair_max_terms(
            self._handle, _which_to_int(which), ctypes.byref(out),
        ), lib)
        return out.value

    def pair_eval(self, bra: int, ket: int,
                  which: str = "full") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate the (bra, ket) CSF pair.

        ``bra`` and ``ket`` are input-order CSF indices in ``range(ncsfs)``,
        the same space as :meth:`csf_steps` and :meth:`csf_index`.

        Returns (ij, pqrs, coef) as three numpy arrays (int32, int32, float64),
        trimmed to the actual term count. ``ij`` is a packed canonical pair
        address, not a pair of indices into this DRT -- see
        docs/PQRS_INDEX_SPEC.md. Do not decode it and feed the result back.
        """
        which_int = _which_to_int(which)
        max_t = self.pair_max_terms(which)
        ij   = np.empty(max_t, dtype=np.int32)
        pqrs = np.empty(max_t, dtype=np.int32)
        coef = np.empty(max_t, dtype=np.float64)
        out_count = ctypes.c_int()
        check_status(lib.hexpr_pair_eval(
            self._handle, int(bra), int(ket), which_int,
            ctypes.byref(out_count),
            _ptr32(ij), _ptr32(pqrs),
            coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ), lib)
        n = out_count.value
        return ij[:n], pqrs[:n], coef[:n]

    # --- block evaluation ---

    def block_count(self, bra_indices, ket_indices, which: str = "full") -> int:
        """Number of expression terms across the rectangular block."""
        bra = _as_int32_1d(bra_indices, "bra_indices")
        ket = _as_int32_1d(ket_indices, "ket_indices")
        out = ctypes.c_int()
        check_status(lib.hexpr_block_count(
            self._handle,
            _ptri(bra), len(bra),
            _ptri(ket), len(ket),
            _which_to_int(which),
            ctypes.byref(out),
        ), lib)
        return out.value

    def block_max_terms(self, n_bra: int, n_ket: int, which: str = "full") -> int:
        """Upper bound on the term count for a block of this size."""
        out = ctypes.c_int()
        check_status(lib.hexpr_block_max_terms(
            self._handle, int(n_bra), int(n_ket), _which_to_int(which),
            ctypes.byref(out),
        ), lib)
        return out.value

    def block_eval(self, bra_indices, ket_indices,
                   which: str = "full") -> Tuple[np.ndarray, ...]:
        """Evaluate the rectangular block.

        ``bra_indices`` and ``ket_indices`` hold input-order CSF indices in
        ``range(ncsfs)``.

        Returns (ij, pqrs, coef, bra_pos, ket_pos) as five numpy arrays.
        bra_pos[i] and ket_pos[i] are 0-based positions into bra_indices /
        ket_indices, not raw CSF indices; they, not ``ij``, are how a term is
        attributed to a pair.
        """
        which_int = _which_to_int(which)
        bra = _as_int32_1d(bra_indices, "bra_indices")
        ket = _as_int32_1d(ket_indices, "ket_indices")
        max_t = self.block_max_terms(len(bra), len(ket), which)
        bra_pos = np.empty(max_t, dtype=np.int32)
        ket_pos = np.empty(max_t, dtype=np.int32)
        ij      = np.empty(max_t, dtype=np.int32)
        pqrs    = np.empty(max_t, dtype=np.int32)
        coef    = np.empty(max_t, dtype=np.float64)
        out_count = ctypes.c_int()
        check_status(lib.hexpr_block_eval(
            self._handle,
            _ptri(bra), len(bra),
            _ptri(ket), len(ket),
            which_int,
            ctypes.byref(out_count),
            _ptr32(bra_pos), _ptr32(ket_pos),
            _ptr32(ij), _ptr32(pqrs),
            coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ), lib)
        n = out_count.value
        return ij[:n], pqrs[:n], coef[:n], bra_pos[:n], ket_pos[:n]

    # --- subspace evaluation ---

    def subspace_count(self, indices, which: str = "full") -> int:
        """Number of expression terms for the square subspace (bra == ket == indices)."""
        idx = _as_int32_1d(indices, "indices")
        out = ctypes.c_int()
        check_status(lib.hexpr_subspace_count(
            self._handle, _ptri(idx), len(idx),
            _which_to_int(which), ctypes.byref(out),
        ), lib)
        return out.value

    def subspace_max_terms(self, n: int, which: str = "full") -> int:
        """Upper bound on the term count for a subspace of size n."""
        out = ctypes.c_int()
        check_status(lib.hexpr_subspace_max_terms(
            self._handle, int(n), _which_to_int(which), ctypes.byref(out),
        ), lib)
        return out.value

    def subspace_eval(self, indices,
                      which: str = "full") -> Tuple[np.ndarray, ...]:
        """Evaluate the square subspace.

        ``indices`` holds input-order CSF indices in ``range(ncsfs)``.

        Returns (ij, pqrs, coef, bra_pos, ket_pos) as five numpy arrays.
        bra_pos[i] and ket_pos[i] are 0-based positions into ``indices``.
        """
        which_int = _which_to_int(which)
        idx = _as_int32_1d(indices, "indices")
        max_t = self.subspace_max_terms(len(idx), which)
        bra_pos = np.empty(max_t, dtype=np.int32)
        ket_pos = np.empty(max_t, dtype=np.int32)
        ij      = np.empty(max_t, dtype=np.int32)
        pqrs    = np.empty(max_t, dtype=np.int32)
        coef    = np.empty(max_t, dtype=np.float64)
        out_count = ctypes.c_int()
        check_status(lib.hexpr_subspace_eval(
            self._handle, _ptri(idx), len(idx),
            which_int,
            ctypes.byref(out_count),
            _ptr32(bra_pos), _ptr32(ket_pos),
            _ptr32(ij), _ptr32(pqrs),
            coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ), lib)
        n = out_count.value
        return ij[:n], pqrs[:n], coef[:n], bra_pos[:n], ket_pos[:n]

    def __del__(self) -> None:
        h = getattr(self, "_handle", None)
        if h is not None:
            lib.hexpr_drt_destroy(h)
            self._handle = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Expression
# ---------------------------------------------------------------------------

class Expression:
    """Brooks energy expression for a given DRT.

    Build via ``hexpr_expr_build``; term arrays are read on demand
    and returned as numpy arrays (0-based indices).
    """

    def __init__(self, drt: Drt) -> None:
        h = lib.hexpr_expr_build(drt._handle)
        if not h:
            raw = lib.hexpr_last_error()
            msg = raw.decode("utf-8", errors="replace") if raw else "hexpr_expr_build returned NULL"
            raise HexprError(msg)
        self._handle: ExprHandle = h

    def nterms_full(self) -> int:
        return int(lib.hexpr_expr_nterms_full(self._handle))

    def nterms_one(self) -> int:
        return int(lib.hexpr_expr_nterms_one(self._handle))

    def terms_full(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (ij, pqrs, coef) for the full (1e+2e) expression.

        Both ij and pqrs are 0-based int32 arrays.  coef is float64.

        ij is a 0-based *packed triangle index* of the (bra, ket) CSF pair, not a
        row index: decode with I = int((sqrt(8*ij + 1) - 1) // 2), J = ij - I*(I+1)//2
        (0-based I >= J), then add 1 in-language for 1-based array frameworks. This
        0-based convention is what the PySCF / Psi4 / Fermi.jl example seams decode.
        The text dump format adds 1 to ij (to match the Ruby corpus); this binary
        API does not, and neither the dump nor this API adjusts pqrs.
        """
        return self._read(HEXPR_EXPR_FULL, self.nterms_full())

    def terms_one(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (ij, pqrs, coef) for the 1-electron subset only."""
        return self._read(HEXPR_EXPR_ONE, self.nterms_one())

    def show_full(self, file=None) -> None:
        """Dump the full (1e+2e) expression to stdout.

        ``file`` must be ``None`` (stdout) in v1.0.0.
        """
        if file is not None:
            raise NotImplementedError(
                "Expression.show_full supports only stdout (file=None) in v1.0.0")
        check_status(lib.hexpr_expr_show_full(self._handle, _c_stdout()), lib)

    def show_one(self, file=None) -> None:
        """Dump the 1-electron expression subset to stdout.

        ``file`` must be ``None`` (stdout) in v1.0.0.
        """
        if file is not None:
            raise NotImplementedError(
                "Expression.show_one supports only stdout (file=None) in v1.0.0")
        check_status(lib.hexpr_expr_show_one(self._handle, _c_stdout()), lib)

    def _read(self, which: int, nterms: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ij   = np.empty(nterms, dtype=np.int32)
        pqrs = np.empty(nterms, dtype=np.int32)
        coef = np.empty(nterms, dtype=np.float64)
        if nterms > 0:
            status = lib.hexpr_expr_read(
                self._handle, which,
                ij.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                pqrs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                coef.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )
            check_status(status, lib)
        return ij, pqrs, coef

    def __del__(self) -> None:
        h = getattr(self, "_handle", None)
        if h is not None:
            lib.hexpr_expr_destroy(h)
            self._handle = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Step-string parsing
# ---------------------------------------------------------------------------

def step_parse(s: str, norb: int) -> np.ndarray:
    """Parse an ``"eudf"``-style step string of length ``norb`` into codes.

    Returns a uint8 array of length ``norb`` with values 0..3
    (E/U/D/F). For example, ``step_parse("fude", 4)`` -> ``[3, 1, 2, 0]``.
    """
    out = np.empty(norb, dtype=np.uint8)
    check_status(lib.hexpr_step_parse(
        s.encode("utf-8"), int(norb),
        out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
    ), lib)
    return out


# ---------------------------------------------------------------------------
# High-level CSF step-code functions
# ---------------------------------------------------------------------------

# The three functions below build a sub-DRT from the supplied rows and then
# address those rows by their INPUT positions, 0..n-1. Earlier versions kept a
# {row -> canonical index} lookup here because the eval API took canonical
# (walk) indices, which differ from input order for a CSFSET sub-DRT. The
# library does that translation now, so the lookup is gone.


def csf_pair_eval(bra_steps, ket_steps,
                  which: str = "full") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a CSF pair given step-code sequences.

    Internally builds a 2-CSF sub-Drt and discards it.

    Parameters
    ----------
    bra_steps, ket_steps : array-like
        1-D uint8 arrays of length norb (values 0..3).
    which : 'full' | 'one'

    Returns
    -------
    (ij, pqrs, coef) as three numpy arrays (int32, int32, float64).
    """
    bra = np.ascontiguousarray(bra_steps, dtype=np.uint8)
    ket = np.ascontiguousarray(ket_steps, dtype=np.uint8)
    if bra.ndim != 1:
        raise ValueError(f"bra_steps must be 1-D, got {bra.ndim}-D")
    if ket.ndim != 1:
        raise ValueError(f"ket_steps must be 1-D, got {ket.ndim}-D")
    norb = len(bra)
    if len(ket) != norb:
        raise ValueError(
            f"bra_steps length {norb} != ket_steps length {len(ket)}")
    stacked = np.stack([bra, ket], axis=0)
    drt = Drt.from_csfset(stacked, norb)
    # Row 0 is the bra, row 1 is the ket, whatever the canonical order is.
    return drt.pair_eval(0, 1, which)


def csf_block_eval(bra_csfs, ket_csfs,
                   which: str = "full") -> Tuple[np.ndarray, ...]:
    """Evaluate a rectangular block given CSF step-code matrices.

    Internally builds a sub-Drt from the union of bra_csfs and ket_csfs.

    Parameters
    ----------
    bra_csfs : array-like
        2-D uint8 array, shape (n_bra, norb).
    ket_csfs : array-like
        2-D uint8 array, shape (n_ket, norb).
    which : 'full' | 'one'

    Returns
    -------
    (ij, pqrs, coef, bra_pos, ket_pos) as five numpy arrays.
    bra_pos[i] and ket_pos[i] are 0-based positions into bra_csfs / ket_csfs.
    """
    bra = np.ascontiguousarray(bra_csfs, dtype=np.uint8)
    ket = np.ascontiguousarray(ket_csfs, dtype=np.uint8)
    if bra.ndim != 2:
        raise ValueError(f"bra_csfs must be 2-D, got {bra.ndim}-D")
    if ket.ndim != 2:
        raise ValueError(f"ket_csfs must be 2-D, got {ket.ndim}-D")
    n_bra, norb = bra.shape
    n_ket, norb2 = ket.shape
    if norb != norb2:
        raise ValueError(
            f"bra_csfs and ket_csfs have different norb: {norb} vs {norb2}")
    all_csfs = np.concatenate([bra, ket], axis=0)
    drt = Drt.from_csfset(all_csfs, norb)
    # The sub-DRT's input rows are bra_csfs followed by ket_csfs, so the two
    # index sets are the two halves of range(n_bra + n_ket).
    bra_indices = np.arange(n_bra, dtype=np.int32)
    ket_indices = np.arange(n_bra, n_bra + n_ket, dtype=np.int32)
    return drt.block_eval(bra_indices, ket_indices, which)


def csf_subspace_eval(csfs,
                      which: str = "full") -> Tuple[np.ndarray, ...]:
    """Evaluate a square subspace given a CSF step-code matrix.

    Parameters
    ----------
    csfs : array-like
        2-D uint8 array, shape (n, norb).
    which : 'full' | 'one'

    Returns
    -------
    (ij, pqrs, coef, bra_pos, ket_pos) as five numpy arrays.
    """
    arr = np.ascontiguousarray(csfs, dtype=np.uint8)
    if arr.ndim != 2:
        raise ValueError(f"csfs must be 2-D, got {arr.ndim}-D")
    n, norb = arr.shape
    drt = Drt.from_csfset(arr, norb)
    indices = np.arange(n, dtype=np.int32)
    return drt.subspace_eval(indices, which)
