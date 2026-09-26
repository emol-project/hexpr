"""Tests for the Expression class against the reference text files.

Index-base note
---------------
The text dump (hexpr_expr_show_full/one) prints ``ij + 1`` (1-based).
The binary API (hexpr_expr_read) returns ``ij`` 0-based.
So when parsing reference text we subtract 1 from the IJ column.

pqrs in the text dump is already 0-based (no +1 is applied by the C code).
This was confirmed by reading hexpr_api.c / expr.c directly:
  fprintf(f, " %7d  %10d  ...", t->ij + 1, t->pqrs, t->coef)
"""
from __future__ import annotations

import os

import numpy as np
import pytest
import hexpr

from conftest import CASES, REF_DIR


# ---------------------------------------------------------------------------
# Reference file presence
# ---------------------------------------------------------------------------

def _reference_readable(path) -> bool:
    """Whether a reference file is present, using the same rule as the C
    harness.

    The C suite (c/tests/test_full.c) decides a reference file is missing
    with ``access(ref, R_OK) != 0`` and prints
    ``SKIP (reference file not bundled)`` instead of FAIL. We mirror that
    exactly with ``os.access(path, os.R_OK)``, so both suites agree on what
    "absent" means (covers non-existence and unreadable files alike). The
    large ``*_med.expr.txt`` files are fetched separately via
    ``reference.tar.gz`` (see INSTALL.md section 5).
    """
    return os.access(path, os.R_OK)


# ---------------------------------------------------------------------------
# Reference file parser
# ---------------------------------------------------------------------------

def parse_expr_file(path) -> tuple:
    """Return (ij, pqrs, coef) arrays from a reference .expr.txt or .expr_one.txt.

    ij is returned 0-based (text is 1-based; we subtract 1).
    pqrs is returned as-is (text is already 0-based).
    """
    ij_list, pqrs_list, coef_list = [], [], []
    in_expr = False
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped == 'Energy Expression':
            in_expr = True
            continue
        if in_expr:
            if stripped.startswith('IJ'):
                continue   # header line
            parts = stripped.split()
            if len(parts) == 3:
                try:
                    ij_list.append(int(parts[0]) - 1)   # 1-based -> 0-based
                    pqrs_list.append(int(parts[1]))      # already 0-based
                    coef_list.append(float(parts[2]))
                except ValueError:
                    pass
    return (
        np.array(ij_list,  dtype=np.int32),
        np.array(pqrs_list, dtype=np.int32),
        np.array(coef_list, dtype=np.float64),
    )


# ---------------------------------------------------------------------------
# Parametrised tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case_name,kwargs", CASES)
def test_terms_full(case_name, kwargs):
    """terms_full() numerically matches the reference .expr.txt."""
    drt  = hexpr.Drt(**kwargs)
    expr = hexpr.Expression(drt)

    got_ij, got_pqrs, got_coef = expr.terms_full()

    ref_path = REF_DIR / f"{case_name}.expr.txt"
    if not _reference_readable(ref_path):
        pytest.skip(
            f"{ref_path.name} not bundled; fetch reference.tar.gz for full coverage")
    ref_ij, ref_pqrs, ref_coef = parse_expr_file(ref_path)

    assert len(got_ij) == len(ref_ij), \
        f"{case_name}: nterms_full {len(got_ij)} != ref {len(ref_ij)}"

    np.testing.assert_array_equal(got_ij,   ref_ij,
        err_msg=f"{case_name}: ij arrays differ")
    np.testing.assert_array_equal(got_pqrs, ref_pqrs,
        err_msg=f"{case_name}: pqrs arrays differ")
    np.testing.assert_allclose(got_coef, ref_coef, rtol=1e-12, atol=1e-14,
        err_msg=f"{case_name}: coef arrays differ")


@pytest.mark.parametrize("case_name,kwargs", CASES)
def test_terms_one(case_name, kwargs):
    """terms_one() numerically matches the reference .expr_one.txt."""
    drt  = hexpr.Drt(**kwargs)
    expr = hexpr.Expression(drt)

    got_ij, got_pqrs, got_coef = expr.terms_one()

    ref_path = REF_DIR / f"{case_name}.expr_one.txt"
    if not _reference_readable(ref_path):
        pytest.skip(
            f"{ref_path.name} not bundled; fetch reference.tar.gz for full coverage")
    ref_ij, ref_pqrs, ref_coef = parse_expr_file(ref_path)

    assert len(got_ij) == len(ref_ij), \
        f"{case_name}: nterms_one {len(got_ij)} != ref {len(ref_ij)}"

    np.testing.assert_array_equal(got_ij,   ref_ij,
        err_msg=f"{case_name}: ij arrays differ")
    np.testing.assert_array_equal(got_pqrs, ref_pqrs,
        err_msg=f"{case_name}: pqrs arrays differ")
    np.testing.assert_allclose(got_coef, ref_coef, rtol=1e-12, atol=1e-14,
        err_msg=f"{case_name}: coef arrays differ")


@pytest.mark.parametrize("case_name,kwargs", CASES)
def test_nterms_one_subset_of_full(case_name, kwargs):
    """nterms_one <= nterms_full and one-electron terms are a subset."""
    drt  = hexpr.Drt(**kwargs)
    expr = hexpr.Expression(drt)
    assert expr.nterms_one() <= expr.nterms_full()
