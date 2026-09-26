"""Tests for Drt.from_csfset."""
from pathlib import Path
import numpy as np
import pytest
import hexpr

REPO_ROOT      = Path(__file__).parent.parent.parent  # .../hexpr/
REF_CSFSET_DIR = REPO_ROOT / "reference_csfset"


def _parse_oracle(path):
    """Parse a reference_csfset/<case>.txt oracle file.

    Returns (norb, rows[uint8 (N,norb)], sections) where each section is
    {'which': str, 'nterms': int, 'terms': sorted list of (bra_pos, ket_pos,
    ij, pqrs, coef)} — bra_pos/ket_pos are input positions (input-order
    contract).
    """
    rows, sections, norb, cur = [], [], None, None
    for line in path.read_text().splitlines():
        t = line.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "norb":
            norb = int(t[1])
        elif t[0] == "ncsfs":
            pass
        elif t[0] == "csf":
            rows.append([int(x) for x in t[2:]])
        elif t[0] == "which":
            cur = {"which": t[1], "terms": []}
            sections.append(cur)
        elif t[0] == "nterms":
            cur["nterms"] = int(t[1])
        else:
            cur["terms"].append(
                (int(t[0]), int(t[1]), int(t[2]), int(t[3]), float(t[4])))
    for sec in sections:
        sec["terms"].sort()
    return norb, np.array(rows, dtype=np.uint8), sections


@pytest.mark.parametrize("case", ["O1", "O2", "O5", "OF", "OFs", "OS"])
def test_csf_subspace_eval_matches_oracle(case):
    """Fixed high-level csf_subspace_eval reproduces the frozen
    reference_csfset/ oracle bit-exactly (input-order bra_pos/ket_pos
    contract; covers superset O2, dedup O5, FOCI OF/OFs, SOCI OS)."""
    path = REF_CSFSET_DIR / f"{case}.txt"
    if not path.is_file():
        pytest.skip(f"oracle file {path} not bundled")
    _, rows, sections = _parse_oracle(path)
    for sec in sections:
        ij, pqrs, coef, bp, kp = hexpr.csf_subspace_eval(rows, which=sec["which"])
        got = sorted(zip(bp.tolist(), kp.tolist(), ij.tolist(),
                         pqrs.tolist(), coef.tolist()))
        exp = sec["terms"]
        assert len(got) == sec["nterms"], \
            f"{case}/{sec['which']}: nterms {len(got)} != {sec['nterms']}"
        for g, e in zip(got, exp):
            assert g[:4] == e[:4], f"{case}/{sec['which']}: key mismatch {g[:4]} vs {e[:4]}"
            assert g[4] == e[4], f"{case}/{sec['which']}: coef {g[4]!r} != {e[4]!r}"


def test_from_csfset_round_trip_full_4o4e():
    """Full CSF set sub-DRT matches the SOCI DRT in shape."""
    full = hexpr.Drt(kind='soci', nve=4, spin_x2=0,
                     ncore=0, nval=4, next_=0)
    assert full.norb == 4
    assert full.ncsfs == 20
    assert full.nnodes == 14

    csfs = full.csfs()
    sub = hexpr.Drt.from_csfset(csfs, norb=4)

    assert sub.norb == full.norb
    assert sub.ncsfs == full.ncsfs
    assert sub.nnodes == full.nnodes
    assert sub.kind == 'csfset'


def test_from_csfset_2_csf_subset_smaller():
    """2-CSF sub-DRT has fewer nodes than full DRT."""
    full = hexpr.Drt(kind='soci', nve=4, spin_x2=0,
                     ncore=0, nval=4, next_=0)
    csfs = full.csfs()
    sub = hexpr.Drt.from_csfset(csfs[:2], norb=4)

    assert sub.norb == 4
    assert sub.ncsfs >= 2
    assert sub.nnodes < full.nnodes
    assert sub.kind == 'csfset'


def test_from_csfset_flat_input():
    """Accepts flat 1-D array of shape (n_csfs * norb,)."""
    full = hexpr.Drt(kind='soci', nve=4, spin_x2=0,
                     ncore=0, nval=4, next_=0)
    csfs_flat = full.csfs().ravel()
    sub = hexpr.Drt.from_csfset(csfs_flat, norb=4)
    assert sub.ncsfs == 20


def test_from_csfset_input_order_and_walk_superset():
    """CSFSET: csfs()/ncsfs speak input order; walk_* exposes the superset."""
    # Non-canonical input that triggers a canonical superset: N=2, M=4.
    inp = np.array([[1, 0, 3, 2],   # uefd
                    [0, 1, 2, 3]],  # eudf
                   dtype=np.uint8)
    sub = hexpr.Drt.from_csfset(inp, norb=4)

    # csfs()/ncsfs == the input rows in input order
    assert sub.ncsfs == 2
    assert np.array_equal(sub.csfs(), inp)

    # walk_* == canonical superset M >= N
    assert sub.walk_ncsfs == 4
    assert sub.walk_ncsfs >= sub.ncsfs
    walk = sub.walk_csfs()
    assert walk.shape == (4, 4)
    # every input row appears somewhere in the canonical walk set
    walk_rows = {tuple(r) for r in walk}
    for row in inp:
        assert tuple(row) in walk_rows


def test_walk_equals_csfs_for_soci():
    """SOCI/FOCI: walk_csfs()/walk_ncsfs == csfs()/ncsfs (no input set)."""
    soci = hexpr.Drt(kind='soci', nve=4, spin_x2=0, ncore=0, nval=4, next_=0)
    assert soci.walk_ncsfs == soci.ncsfs
    assert np.array_equal(soci.walk_csfs(), soci.csfs())


def test_from_csfset_invalid_step_raises():
    """Step value outside 0..3 raises HexprError."""
    bad = np.array([[0, 1, 2, 4]], dtype=np.uint8)  # step 4 invalid
    with pytest.raises(hexpr.HexprError):
        hexpr.Drt.from_csfset(bad, norb=4)


def test_from_csfset_shape_mismatch_raises():
    """csfs.shape[1] != norb raises ValueError."""
    arr = np.zeros((2, 3), dtype=np.uint8)
    with pytest.raises(ValueError):
        hexpr.Drt.from_csfset(arr, norb=4)
