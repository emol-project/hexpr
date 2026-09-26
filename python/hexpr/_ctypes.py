"""ctypes function signatures for libhexpr.so.

Locates the shared library using this search order:
  1. HEXPR_LIBRARY environment variable (full path to .so)
  2. ../../c/libhexpr.so relative to this file (in-tree build)
  3. LD_LIBRARY_PATH / system linker path via ctypes.util.find_library
"""
from __future__ import annotations

import ctypes
import ctypes.util
import os
import sys
from ctypes import c_char_p, c_int, c_int32, c_int64, c_double, c_void_p, POINTER


def _lib_names(stem: str) -> list[str]:
    """In-tree shared-lib filenames to try, most-specific first."""
    if sys.platform == "darwin":
        return [stem + ".dylib", stem + ".so"]
    return [stem + ".so"]


# --- opaque handle types ---

class DrtHandle(c_void_p):
    """Opaque handle for hexpr_drt_t *."""


class ExprHandle(c_void_p):
    """Opaque handle for hexpr_expr_t *."""


# --- status constants ---
HEXPR_OK              = 0
HEXPR_ERR_INVALID_ARG = 1
HEXPR_ERR_OOM         = 2
HEXPR_ERR_IO          = 3
HEXPR_ERR_RANGE       = 4
HEXPR_ERR_NOT_BUILT   = 5
HEXPR_ERR_INTERNAL    = 99

# --- DRT kind ---
HEXPR_DRT_SOCI   = 0
HEXPR_DRT_FOCI   = 1
HEXPR_DRT_CSFSET = 2

# --- step codes ---
HEXPR_STEP_E = 0
HEXPR_STEP_U = 1
HEXPR_STEP_D = 2
HEXPR_STEP_F = 3

# --- expression selector ---
HEXPR_EXPR_FULL = 0
HEXPR_EXPR_ONE  = 1


def _find_library() -> str:
    env = os.environ.get("HEXPR_LIBRARY")
    if env:
        return env

    here = os.path.dirname(os.path.abspath(__file__))
    for name in _lib_names("libhexpr"):
        cand = os.path.normpath(os.path.join(here, "..", "..", "c", name))
        if os.path.isfile(cand):
            return cand

    found = ctypes.util.find_library("hexpr")
    if found:
        return found

    soname = _lib_names("libhexpr")[0]
    envvar = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    raise OSError(
        f"Cannot locate {soname}. Set the HEXPR_LIBRARY environment variable "
        f"to the full path of {soname}, or add its directory to {envvar}."
    )


def _load() -> ctypes.CDLL:
    lib = ctypes.CDLL(_find_library())

    # --- version ---
    lib.hexpr_version_string.restype = c_char_p
    lib.hexpr_version_string.argtypes = []
    lib.hexpr_version_major.restype = c_int
    lib.hexpr_version_major.argtypes = []
    lib.hexpr_version_minor.restype = c_int
    lib.hexpr_version_minor.argtypes = []
    lib.hexpr_version_patch.restype = c_int
    lib.hexpr_version_patch.argtypes = []

    # --- last error ---
    lib.hexpr_last_error.restype = c_char_p
    lib.hexpr_last_error.argtypes = []

    # --- lifecycle ---
    lib.hexpr_init.restype = c_int
    lib.hexpr_init.argtypes = []
    lib.hexpr_shutdown.restype = None
    lib.hexpr_shutdown.argtypes = []

    # --- DRT create/destroy ---
    lib.hexpr_drt_create.restype = DrtHandle
    lib.hexpr_drt_create.argtypes = [c_int, c_int, c_int, c_int, c_int, c_int]
    lib.hexpr_drt_destroy.restype = None
    lib.hexpr_drt_destroy.argtypes = [DrtHandle]

    # --- DRT from CSF set ---
    lib.hexpr_drt_from_csfset.restype  = DrtHandle
    lib.hexpr_drt_from_csfset.argtypes = [c_int, POINTER(ctypes.c_uint8), c_int]

    # --- DRT queries ---
    lib.hexpr_drt_norb.restype   = c_int
    lib.hexpr_drt_norb.argtypes  = [DrtHandle]
    lib.hexpr_drt_nnodes.restype = c_int
    lib.hexpr_drt_nnodes.argtypes = [DrtHandle]
    lib.hexpr_drt_ncsfs.restype  = c_int
    lib.hexpr_drt_ncsfs.argtypes = [DrtHandle]
    lib.hexpr_drt_walk_ncsfs.restype  = c_int
    lib.hexpr_drt_walk_ncsfs.argtypes = [DrtHandle]
    lib.hexpr_drt_kind.restype   = c_int
    lib.hexpr_drt_kind.argtypes  = [DrtHandle]
    lib.hexpr_drt_ncore.restype  = c_int
    lib.hexpr_drt_ncore.argtypes = [DrtHandle]
    lib.hexpr_drt_nval.restype   = c_int
    lib.hexpr_drt_nval.argtypes  = [DrtHandle]
    lib.hexpr_drt_next.restype   = c_int
    lib.hexpr_drt_next.argtypes  = [DrtHandle]
    lib.hexpr_drt_nve_state.restype  = c_int
    lib.hexpr_drt_nve_state.argtypes = [DrtHandle]
    lib.hexpr_drt_spin_x2.restype    = c_int
    lib.hexpr_drt_spin_x2.argtypes   = [DrtHandle]

    # --- CSF buffer (owned by DRT) ---
    lib.hexpr_drt_csfs.restype  = POINTER(ctypes.c_uint8)
    lib.hexpr_drt_csfs.argtypes = [DrtHandle]
    lib.hexpr_drt_walk_csfs.restype  = POINTER(ctypes.c_uint8)
    lib.hexpr_drt_walk_csfs.argtypes = [DrtHandle]

    # --- node info ---
    lib.hexpr_drt_node_info.restype  = c_int
    lib.hexpr_drt_node_info.argtypes = [DrtHandle, c_int,
                                        POINTER(c_int), POINTER(c_int), POINTER(c_int)]

    # --- expression build/destroy ---
    lib.hexpr_expr_build.restype  = ExprHandle
    lib.hexpr_expr_build.argtypes = [DrtHandle]
    lib.hexpr_expr_destroy.restype  = None
    lib.hexpr_expr_destroy.argtypes = [ExprHandle]

    # --- term counts ---
    lib.hexpr_expr_nterms_full.restype  = c_int64
    lib.hexpr_expr_nterms_full.argtypes = [ExprHandle]
    lib.hexpr_expr_nterms_one.restype   = c_int64
    lib.hexpr_expr_nterms_one.argtypes  = [ExprHandle]

    # --- bulk read ---
    lib.hexpr_expr_read.restype  = c_int
    lib.hexpr_expr_read.argtypes = [ExprHandle, c_int,
                                    POINTER(c_int32), POINTER(c_int32), POINTER(c_double)]

    # --- CSF lookup ---
    lib.hexpr_csf_index.restype  = c_int
    lib.hexpr_csf_index.argtypes = [DrtHandle, POINTER(ctypes.c_uint8), POINTER(c_int)]

    # --- CSF steps / step-string parse ---
    lib.hexpr_csf_steps.restype  = c_int
    lib.hexpr_csf_steps.argtypes = [DrtHandle, c_int, POINTER(ctypes.c_uint8)]
    lib.hexpr_step_parse.restype  = c_int
    lib.hexpr_step_parse.argtypes = [c_char_p, c_int, POINTER(ctypes.c_uint8)]

    # --- CSF -> determinant expansion ---
    lib.hexpr_csf_ndets.restype  = c_int
    lib.hexpr_csf_ndets.argtypes = [DrtHandle, c_int, c_int, POINTER(c_int)]

    lib.hexpr_csf_to_dets.restype  = c_int
    lib.hexpr_csf_to_dets.argtypes = [DrtHandle, c_int, c_int,
                                      POINTER(c_int),
                                      POINTER(ctypes.c_uint64),
                                      POINTER(ctypes.c_uint64),
                                      POINTER(c_double)]

    # --- show (dump to a C FILE* stream) ---
    lib.hexpr_drt_show.restype  = c_int
    lib.hexpr_drt_show.argtypes = [DrtHandle, c_void_p]
    lib.hexpr_expr_show_full.restype  = c_int
    lib.hexpr_expr_show_full.argtypes = [ExprHandle, c_void_p]
    lib.hexpr_expr_show_one.restype  = c_int
    lib.hexpr_expr_show_one.argtypes = [ExprHandle, c_void_p]

    # --- per-pair evaluation ---
    lib.hexpr_pair_count.restype  = c_int
    lib.hexpr_pair_count.argtypes = [DrtHandle, c_int, c_int, c_int, POINTER(c_int)]

    lib.hexpr_pair_eval.restype  = c_int
    lib.hexpr_pair_eval.argtypes = [DrtHandle, c_int, c_int, c_int, POINTER(c_int),
                                     POINTER(c_int32), POINTER(c_int32), POINTER(c_double)]

    lib.hexpr_pair_max_terms.restype  = c_int
    lib.hexpr_pair_max_terms.argtypes = [DrtHandle, c_int, POINTER(c_int)]

    # --- block evaluation ---
    lib.hexpr_block_count.restype  = c_int
    lib.hexpr_block_count.argtypes = [DrtHandle,
                                       POINTER(c_int), c_int,
                                       POINTER(c_int), c_int,
                                       c_int, POINTER(c_int)]

    lib.hexpr_block_eval.restype  = c_int
    lib.hexpr_block_eval.argtypes = [DrtHandle,
                                      POINTER(c_int), c_int,
                                      POINTER(c_int), c_int,
                                      c_int, POINTER(c_int),
                                      POINTER(c_int32), POINTER(c_int32),
                                      POINTER(c_int32), POINTER(c_int32),
                                      POINTER(c_double)]

    lib.hexpr_block_max_terms.restype  = c_int
    lib.hexpr_block_max_terms.argtypes = [DrtHandle, c_int, c_int, c_int, POINTER(c_int)]

    # --- subspace evaluation ---
    lib.hexpr_subspace_count.restype  = c_int
    lib.hexpr_subspace_count.argtypes = [DrtHandle, POINTER(c_int), c_int, c_int, POINTER(c_int)]

    lib.hexpr_subspace_eval.restype  = c_int
    lib.hexpr_subspace_eval.argtypes = [DrtHandle,
                                         POINTER(c_int), c_int,
                                         c_int, POINTER(c_int),
                                         POINTER(c_int32), POINTER(c_int32),
                                         POINTER(c_int32), POINTER(c_int32),
                                         POINTER(c_double)]

    lib.hexpr_subspace_max_terms.restype  = c_int
    lib.hexpr_subspace_max_terms.argtypes = [DrtHandle, c_int, c_int, POINTER(c_int)]

    # Initialize the library (idempotent; sets up Wigner tables).
    lib.hexpr_init()

    return lib


lib = _load()
