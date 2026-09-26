"""Tests for the Drt class against the reference text files."""
from __future__ import annotations

import re
import numpy as np
import pytest
import hexpr

from conftest import CASES, REF_DIR


# ---------------------------------------------------------------------------
# Reference file parsers
# ---------------------------------------------------------------------------

def parse_drt_header(path) -> tuple:
    """Return (norb, nnodes, ncsfs) from the .drt.txt header."""
    text = path.read_text()
    norb   = int(re.search(r'NORB\s*=\s*(\d+)', text).group(1))
    nnodes = int(re.search(r'Number of generated nodes\s*=\s*(\d+)', text).group(1))
    ncsfs  = int(re.search(r'Number of generated csfs\s*=\s*(\d+)', text).group(1))
    return norb, nnodes, ncsfs


_STEP_MAP = {'e': 0, 'u': 1, 'd': 2, 'f': 3}


def parse_csf_list(path) -> np.ndarray:
    """Return (ncsfs, norb) uint8 array of step codes from a .drt.txt file."""
    rows = []
    in_csf = False
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped == 'CSF List':
            in_csf = True
            continue
        if in_csf:
            if stripped.startswith('['):
                steps = re.findall(r'"([eudf])"', stripped)
                rows.append([_STEP_MAP[s] for s in steps])
            else:
                break   # blank line or new section
    return np.array(rows, dtype=np.uint8)


# ---------------------------------------------------------------------------
# Parametrised tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case_name,kwargs", CASES)
def test_drt_counts(case_name, kwargs):
    """norb, nnodes, ncsfs match the reference header."""
    drt = hexpr.Drt(**kwargs)
    ref_path = REF_DIR / f"{case_name}.drt.txt"
    norb, nnodes, ncsfs = parse_drt_header(ref_path)

    assert drt.norb   == norb,   f"{case_name}: norb mismatch"
    assert drt.nnodes == nnodes, f"{case_name}: nnodes mismatch"
    assert drt.ncsfs  == ncsfs,  f"{case_name}: ncsfs mismatch"


@pytest.mark.parametrize("case_name,kwargs", CASES)
def test_drt_csfs(case_name, kwargs):
    """csfs() matches the CSF list in the reference .drt.txt."""
    drt = hexpr.Drt(**kwargs)
    ref_path = REF_DIR / f"{case_name}.drt.txt"
    ref_csfs = parse_csf_list(ref_path)

    got = drt.csfs()
    assert got.shape == ref_csfs.shape, \
        f"{case_name}: CSF shape {got.shape} != ref {ref_csfs.shape}"
    np.testing.assert_array_equal(got, ref_csfs,
        err_msg=f"{case_name}: CSF step codes differ")


@pytest.mark.parametrize("case_name,kwargs", CASES)
def test_drt_properties(case_name, kwargs):
    """kind, ncore, nval, next_, spin_x2 round-trip correctly."""
    drt = hexpr.Drt(**kwargs)
    assert drt.kind    == kwargs["kind"]
    assert drt.ncore   == kwargs["ncore"]
    assert drt.nval    == kwargs["nval"]
    assert drt.next_   == kwargs["next_"]
    assert drt.spin_x2 == kwargs["spin_x2"]
