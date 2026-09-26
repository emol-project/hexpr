"""Tests for per-pair / block / subspace evaluation."""
from __future__ import annotations

import numpy as np
import pytest
import hexpr

from conftest import CASES, REF_DIR


# ---------------------------------------------------------------------------
# Helpers shared across tests
# ---------------------------------------------------------------------------

def triangle(a: int, b: int) -> int:
    hi, lo = max(a, b), min(a, b)
    return hi * (hi + 1) // 2 + lo


def parse_expr_file(path) -> tuple:
    """Return (ij, pqrs, coef) from a reference .expr.txt (ij is 0-based)."""
    ij_list, pqrs_list, coef_list = [], [], []
    in_expr = False
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped == "Energy Expression":
            in_expr = True
            continue
        if in_expr:
            if stripped.startswith("IJ"):
                continue
            parts = stripped.split()
            if len(parts) == 3:
                try:
                    ij_list.append(int(parts[0]) - 1)
                    pqrs_list.append(int(parts[1]))
                    coef_list.append(float(parts[2]))
                except ValueError:
                    pass
    return (
        np.array(ij_list,   dtype=np.int32),
        np.array(pqrs_list, dtype=np.int32),
        np.array(coef_list, dtype=np.float64),
    )


def aggregate_all_pairs(drt: hexpr.Drt, which: str = "full"):
    """Aggregate pair_eval over all (bra, ket), first-wins dedup on (ij, pqrs)."""
    ncsfs = drt.ncsfs
    seen: dict[tuple, float] = {}
    for bra in range(ncsfs):
        for ket in range(ncsfs):
            ij_arr, pqrs_arr, coef_arr = drt.pair_eval(bra, ket, which)
            for k in range(len(ij_arr)):
                key = (int(ij_arr[k]), int(pqrs_arr[k]))
                if key not in seen:
                    seen[key] = float(coef_arr[k])
    if not seen:
        return (np.empty(0, np.int32), np.empty(0, np.int32), np.empty(0, np.float64))
    keys = sorted(seen)
    ij_out   = np.array([k[0] for k in keys], dtype=np.int32)
    pqrs_out = np.array([k[1] for k in keys], dtype=np.int32)
    coef_out = np.array([seen[k] for k in keys], dtype=np.float64)
    return ij_out, pqrs_out, coef_out


# ---------------------------------------------------------------------------
# pair_count / pair_max_terms / pair_eval -- basic
# ---------------------------------------------------------------------------

class TestPairBasic:
    def setup_method(self):
        self.drt = hexpr.Drt(kind="soci", nve=4, spin_x2=0,
                              ncore=0, nval=4, next_=0)

    def test_pair_count_is_nonneg(self):
        c = self.drt.pair_count(0, 1)
        assert c >= 0

    def test_pair_max_terms_full_is_positive(self):
        m = self.drt.pair_max_terms("full")
        assert m > 0

    def test_pair_max_terms_one_leq_full(self):
        assert self.drt.pair_max_terms("one") <= self.drt.pair_max_terms("full")

    def test_pair_eval_returns_correct_dtypes(self):
        ij, pqrs, coef = self.drt.pair_eval(0, 1)
        assert ij.dtype   == np.int32
        assert pqrs.dtype == np.int32
        assert coef.dtype == np.float64

    def test_pair_eval_trimmed_length_matches_count(self):
        cnt = self.drt.pair_count(0, 1)
        ij, pqrs, coef = self.drt.pair_eval(0, 1)
        assert len(ij) == len(pqrs) == len(coef) == cnt

    def test_pair_eval_diagonal_has_terms(self):
        # diagonal (bra == ket) should produce at least some terms
        total = sum(len(self.drt.pair_eval(i, i)[0])
                    for i in range(self.drt.ncsfs))
        assert total > 0

    def test_pair_eval_one_subset_of_full(self):
        ij_f, pqrs_f, _ = self.drt.pair_eval(0, 1, "full")
        ij_o, pqrs_o, _ = self.drt.pair_eval(0, 1, "one")
        assert len(ij_o) <= len(ij_f)
        # every (ij, pqrs) in 'one' must appear in 'full'
        full_keys = {(int(ij_f[k]), int(pqrs_f[k])) for k in range(len(ij_f))}
        for k in range(len(ij_o)):
            assert (int(ij_o[k]), int(pqrs_o[k])) in full_keys

    def test_pair_eval_ij_equals_triangle(self):
        """ij values returned by pair_eval equal triangle(bra, ket)."""
        for bra in range(3):
            for ket in range(3):
                ij_arr, _, _ = self.drt.pair_eval(bra, ket)
                expected = triangle(bra, ket)
                assert all(v == expected for v in ij_arr.tolist()), \
                    f"pair ({bra},{ket}): expected ij={expected}, got {ij_arr.tolist()}"


# ---------------------------------------------------------------------------
# Aggregate test: sum of all pairs matches expr_build
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case_name,kwargs", [
    ("s_full_4o4e", dict(kind="soci", nve=4, spin_x2=0, ncore=0, nval=4, next_=0)),
    ("t_full_4o4e", dict(kind="soci", nve=4, spin_x2=2, ncore=0, nval=4, next_=0)),
    ("s_full_6o6e", dict(kind="soci", nve=6, spin_x2=0, ncore=0, nval=6, next_=0)),
])
def test_pair_eval_aggregate_matches_reference(case_name, kwargs):
    """Aggregate of pair_eval over all pairs matches the reference .expr.txt."""
    drt = hexpr.Drt(**kwargs)
    got_ij, got_pqrs, got_coef = aggregate_all_pairs(drt, "full")

    ref_path = REF_DIR / f"{case_name}.expr.txt"
    ref_ij, ref_pqrs, ref_coef = parse_expr_file(ref_path)
    # reference is already sorted by (ij, pqrs)
    ref_keys = sorted(zip(ref_ij.tolist(), ref_pqrs.tolist()))
    ref_ij_s   = np.array([k[0] for k in ref_keys], dtype=np.int32)
    ref_pqrs_s = np.array([k[1] for k in ref_keys], dtype=np.int32)
    ref_coef_s = np.array(
        [ref_coef[i] for i in np.lexsort((ref_pqrs, ref_ij))],
        dtype=np.float64)

    assert len(got_ij) == len(ref_ij_s), \
        f"{case_name}: term count {len(got_ij)} != ref {len(ref_ij_s)}"
    np.testing.assert_array_equal(got_ij,   ref_ij_s,
        err_msg=f"{case_name}: ij arrays differ")
    np.testing.assert_array_equal(got_pqrs, ref_pqrs_s,
        err_msg=f"{case_name}: pqrs arrays differ")
    np.testing.assert_allclose(got_coef, ref_coef_s, rtol=1e-10, atol=1e-12,
        err_msg=f"{case_name}: coef arrays differ")


# ---------------------------------------------------------------------------
# block_count / block_max_terms / block_eval
# ---------------------------------------------------------------------------

class TestBlockEval:
    def setup_method(self):
        self.drt = hexpr.Drt(kind="soci", nve=4, spin_x2=0,
                              ncore=0, nval=4, next_=0)

    def test_block_count_nonneg(self):
        c = self.drt.block_count([0, 1], [2, 3])
        assert c >= 0

    def test_block_max_terms_geq_count(self):
        bra = [0, 1]
        ket = [2, 3]
        mx = self.drt.block_max_terms(len(bra), len(ket))
        cnt = self.drt.block_count(bra, ket)
        assert mx >= cnt

    def test_block_eval_returns_five_arrays(self):
        result = self.drt.block_eval([0, 1], [2, 3])
        assert len(result) == 5

    def test_block_eval_dtypes(self):
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.block_eval([0, 1], [2, 3])
        assert ij.dtype      == np.int32
        assert pqrs.dtype    == np.int32
        assert coef.dtype    == np.float64
        assert bra_pos.dtype == np.int32
        assert ket_pos.dtype == np.int32

    def test_block_eval_length_matches_count(self):
        bra, ket = [0, 1, 2], [3, 4]
        cnt = self.drt.block_count(bra, ket)
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.block_eval(bra, ket)
        assert len(ij) == len(pqrs) == len(coef) == len(bra_pos) == len(ket_pos) == cnt

    def test_block_eval_bra_pos_ket_pos_are_within_array_positions(self):
        """bra_pos values must be in 0..n_bra-1, ket_pos in 0..n_ket-1."""
        bra_list = [5, 9, 12]
        ket_list = [2, 7]
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.block_eval(bra_list, ket_list)
        assert all(0 <= v < len(bra_list) for v in bra_pos.tolist()), \
            f"bra_pos out of range: {bra_pos.tolist()}"
        assert all(0 <= v < len(ket_list) for v in ket_pos.tolist()), \
            f"ket_pos out of range: {ket_pos.tolist()}"

    def test_block_eval_pos_values_not_raw_csf_indices(self):
        """bra_pos should be 0-based into bra_indices, not the raw CSF indices.

        If bra_indices={5,9,12}, bra_pos values must be in {0,1,2}, not in {5,9,12}.
        """
        bra_list = [5, 9, 12]
        ket_list = [2, 7]
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.block_eval(bra_list, ket_list)
        if len(bra_pos) > 0:
            assert set(bra_pos.tolist()).issubset({0, 1, 2}), \
                f"Expected bra_pos in {{0,1,2}}, got values {set(bra_pos.tolist())}"

    def test_block_eval_matches_individual_pairs(self):
        """block_eval result agrees with per-pair pair_eval calls."""
        bra_list = [0, 1]
        ket_list = [2, 3]
        ij_b, pqrs_b, coef_b, bra_pos, ket_pos = \
            self.drt.block_eval(bra_list, ket_list)

        # Reconstruct from individual pair_eval calls
        for ip, bra in enumerate(bra_list):
            for jp, ket in enumerate(ket_list):
                ij_p, pqrs_p, coef_p = self.drt.pair_eval(bra, ket)
                # Extract block terms for this (ip, jp)
                mask = (bra_pos == ip) & (ket_pos == jp)
                ij_got   = ij_b[mask]
                pqrs_got = pqrs_b[mask]
                coef_got = coef_b[mask]
                # Sort both by pqrs for comparison
                order_b = np.argsort(pqrs_got)
                order_p = np.argsort(pqrs_p)
                np.testing.assert_array_equal(
                    ij_got[order_b], ij_p[order_p],
                    err_msg=f"ij mismatch for bra={bra},ket={ket}")
                np.testing.assert_array_equal(
                    pqrs_got[order_b], pqrs_p[order_p],
                    err_msg=f"pqrs mismatch for bra={bra},ket={ket}")
                np.testing.assert_allclose(
                    coef_got[order_b], coef_p[order_p], rtol=1e-13,
                    err_msg=f"coef mismatch for bra={bra},ket={ket}")

    def test_block_eval_empty_bra(self):
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.block_eval([], [0, 1])
        assert len(ij) == 0

    def test_block_eval_empty_ket(self):
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.block_eval([0, 1], [])
        assert len(ij) == 0


# ---------------------------------------------------------------------------
# subspace_count / subspace_max_terms / subspace_eval
# ---------------------------------------------------------------------------

class TestSubspaceEval:
    def setup_method(self):
        self.drt = hexpr.Drt(kind="soci", nve=4, spin_x2=0,
                              ncore=0, nval=4, next_=0)

    def test_subspace_count_nonneg(self):
        c = self.drt.subspace_count([0, 1, 2])
        assert c >= 0

    def test_subspace_eval_matches_block_eval(self):
        """subspace_eval([a,b,c]) must equal block_eval([a,b,c],[a,b,c])."""
        idx = [0, 3, 7]
        ij_s, pqrs_s, coef_s, bra_s, ket_s = self.drt.subspace_eval(idx)
        ij_b, pqrs_b, coef_b, bra_b, ket_b = self.drt.block_eval(idx, idx)

        np.testing.assert_array_equal(ij_s,   ij_b)
        np.testing.assert_array_equal(pqrs_s, pqrs_b)
        np.testing.assert_allclose(coef_s, coef_b, rtol=1e-13)
        np.testing.assert_array_equal(bra_s, bra_b)
        np.testing.assert_array_equal(ket_s, ket_b)

    def test_subspace_eval_pos_within_range(self):
        idx = [0, 1, 2, 5]
        ij, pqrs, coef, bra_pos, ket_pos = self.drt.subspace_eval(idx)
        assert all(0 <= v < len(idx) for v in bra_pos.tolist())
        assert all(0 <= v < len(idx) for v in ket_pos.tolist())

    def test_subspace_max_terms_geq_count(self):
        idx = [0, 1, 2]
        mx  = self.drt.subspace_max_terms(len(idx))
        cnt = self.drt.subspace_count(idx)
        assert mx >= cnt


# ---------------------------------------------------------------------------
# High-level csf_pair_eval / csf_block_eval / csf_subspace_eval
# ---------------------------------------------------------------------------

class TestHighLevelFunctions:
    def setup_method(self):
        self.drt = hexpr.Drt(kind="soci", nve=4, spin_x2=0,
                              ncore=0, nval=4, next_=0)
        self.csfs = self.drt.csfs()

    def test_csf_pair_eval_matches_low_level(self):
        bra, ket = 0, 3
        bra_steps = self.csfs[bra]
        ket_steps = self.csfs[ket]

        ij_hi, pqrs_hi, coef_hi = hexpr.csf_pair_eval(bra_steps, ket_steps)
        ij_lo, pqrs_lo, coef_lo = self.drt.pair_eval(bra, ket)

        # Sort both by pqrs for comparison (order may differ between sub-DRT)
        order_hi = np.argsort(pqrs_hi)
        order_lo = np.argsort(pqrs_lo)
        assert len(ij_hi) == len(ij_lo), \
            f"term count mismatch: high={len(ij_hi)}, low={len(ij_lo)}"
        np.testing.assert_array_equal(pqrs_hi[order_hi], pqrs_lo[order_lo])
        np.testing.assert_allclose(coef_hi[order_hi], coef_lo[order_lo], rtol=1e-13)

    def test_csf_pair_eval_diagonal(self):
        """bra == ket case returns same result as low-level."""
        idx = 5
        steps = self.csfs[idx]
        ij_hi, pqrs_hi, coef_hi = hexpr.csf_pair_eval(steps, steps)
        ij_lo, pqrs_lo, coef_lo = self.drt.pair_eval(idx, idx)
        assert len(ij_hi) == len(ij_lo)
        order_hi = np.argsort(pqrs_hi)
        order_lo = np.argsort(pqrs_lo)
        np.testing.assert_array_equal(pqrs_hi[order_hi], pqrs_lo[order_lo])
        np.testing.assert_allclose(coef_hi[order_hi], coef_lo[order_lo], rtol=1e-13)

    def test_csf_pair_eval_reordering(self):
        """bra,ket whose 2-CSF sub-DRT reorders (canonical order != input order).

        The eval API takes input-order indices, so csf_pair_eval passes rows
        0 and 1 and the library maps them onto the reordered canonical walk.
        Feeding the CANONICAL index here would now evaluate the wrong pair --
        the mapping moved into the library and is no longer done here.
        Sub-DRT for (csfs[1], csfs[0]) swaps the two rows in canonical order;
        the assertion below keeps the fixture honest.
        """
        bra, ket = 1, 0
        sub = hexpr.Drt.from_csfset(self.csfs[[bra, ket]], self.csfs.shape[1])
        assert not np.array_equal(sub.csfs(), sub.walk_csfs()), \
            "fixture no longer reorders; pick a different (bra,ket)"
        ij_hi, pqrs_hi, coef_hi = hexpr.csf_pair_eval(self.csfs[bra], self.csfs[ket])
        ij_lo, pqrs_lo, coef_lo = self.drt.pair_eval(bra, ket)
        assert len(ij_hi) == len(ij_lo), \
            f"term count mismatch: high={len(ij_hi)}, low={len(ij_lo)}"
        order_hi = np.argsort(pqrs_hi)
        order_lo = np.argsort(pqrs_lo)
        np.testing.assert_array_equal(pqrs_hi[order_hi], pqrs_lo[order_lo])
        np.testing.assert_allclose(coef_hi[order_hi], coef_lo[order_lo], rtol=1e-13)

    def test_csf_block_eval_matches_low_level(self):
        bra_list = [0, 2]
        ket_list = [1, 3, 4]
        bra_csfs = self.csfs[bra_list]
        ket_csfs = self.csfs[ket_list]

        ij_hi, pqrs_hi, coef_hi, bra_pos_hi, ket_pos_hi = \
            hexpr.csf_block_eval(bra_csfs, ket_csfs)
        ij_lo, pqrs_lo, coef_lo, bra_pos_lo, ket_pos_lo = \
            self.drt.block_eval(bra_list, ket_list)

        # Sort both by (bra_pos, ket_pos, pqrs) for comparison
        def sort_key(bp, kp, pq):
            return np.lexsort((pq, kp, bp))

        if len(ij_hi) == 0 and len(ij_lo) == 0:
            return

        assert len(ij_hi) == len(ij_lo), \
            f"term count mismatch: high={len(ij_hi)}, low={len(ij_lo)}"

        ord_hi = sort_key(bra_pos_hi, ket_pos_hi, pqrs_hi)
        ord_lo = sort_key(bra_pos_lo, ket_pos_lo, pqrs_lo)

        np.testing.assert_array_equal(bra_pos_hi[ord_hi], bra_pos_lo[ord_lo])
        np.testing.assert_array_equal(ket_pos_hi[ord_hi], ket_pos_lo[ord_lo])
        np.testing.assert_array_equal(pqrs_hi[ord_hi], pqrs_lo[ord_lo])
        np.testing.assert_allclose(coef_hi[ord_hi], coef_lo[ord_lo], rtol=1e-13)

    def test_csf_subspace_eval_matches_low_level(self):
        idx = [0, 1, 4]
        sub_csfs = self.csfs[idx]

        ij_hi, pqrs_hi, coef_hi, bra_pos_hi, ket_pos_hi = \
            hexpr.csf_subspace_eval(sub_csfs)
        ij_lo, pqrs_lo, coef_lo, bra_pos_lo, ket_pos_lo = \
            self.drt.subspace_eval(idx)

        if len(ij_hi) == 0 and len(ij_lo) == 0:
            return

        assert len(ij_hi) == len(ij_lo), \
            f"term count mismatch: high={len(ij_hi)}, low={len(ij_lo)}"

        def sort_key(bp, kp, pq):
            return np.lexsort((pq, kp, bp))

        ord_hi = sort_key(bra_pos_hi, ket_pos_hi, pqrs_hi)
        ord_lo = sort_key(bra_pos_lo, ket_pos_lo, pqrs_lo)

        np.testing.assert_array_equal(bra_pos_hi[ord_hi], bra_pos_lo[ord_lo])
        np.testing.assert_array_equal(ket_pos_hi[ord_hi], ket_pos_lo[ord_lo])
        np.testing.assert_array_equal(pqrs_hi[ord_hi], pqrs_lo[ord_lo])
        np.testing.assert_allclose(coef_hi[ord_hi], coef_lo[ord_lo], rtol=1e-13)

    def test_csf_pair_eval_returns_correct_dtypes(self):
        ij, pqrs, coef = hexpr.csf_pair_eval(self.csfs[0], self.csfs[1])
        assert ij.dtype   == np.int32
        assert pqrs.dtype == np.int32
        assert coef.dtype == np.float64

    def test_csf_block_eval_bra_pos_within_range(self):
        bra_csfs = self.csfs[[0, 1, 2]]
        ket_csfs = self.csfs[[3, 4]]
        ij, pqrs, coef, bra_pos, ket_pos = hexpr.csf_block_eval(bra_csfs, ket_csfs)
        assert all(0 <= v < 3 for v in bra_pos.tolist())
        assert all(0 <= v < 2 for v in ket_pos.tolist())


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------

class TestErrorPaths:
    def setup_method(self):
        self.drt = hexpr.Drt(kind="soci", nve=4, spin_x2=0,
                              ncore=0, nval=4, next_=0)

    def test_invalid_which_pair_count(self):
        with pytest.raises(ValueError, match="which"):
            self.drt.pair_count(0, 1, which="bogus")

    def test_invalid_which_pair_eval(self):
        with pytest.raises(ValueError, match="which"):
            self.drt.pair_eval(0, 1, which="TWO")

    def test_invalid_which_block_eval(self):
        with pytest.raises(ValueError, match="which"):
            self.drt.block_eval([0], [1], which="")

    def test_invalid_which_subspace_eval(self):
        with pytest.raises(ValueError, match="which"):
            self.drt.subspace_eval([0, 1], which="half")

    def test_out_of_range_bra_pair_count(self):
        with pytest.raises(hexpr.RangeError):
            self.drt.pair_count(9999, 0)

    def test_out_of_range_ket_pair_eval(self):
        with pytest.raises(hexpr.RangeError):
            self.drt.pair_eval(0, 9999)

    def test_out_of_range_index_block_count(self):
        with pytest.raises(hexpr.RangeError):
            self.drt.block_count([0, 9999], [1])

    def test_2d_bra_indices_raises(self):
        with pytest.raises(ValueError, match="1-D"):
            self.drt.block_eval(np.zeros((2, 2), dtype=np.int32), [0])

    def test_2d_ket_indices_raises(self):
        with pytest.raises(ValueError, match="1-D"):
            self.drt.block_eval([0], np.zeros((2, 2), dtype=np.int32))

    def test_2d_subspace_indices_raises(self):
        with pytest.raises(ValueError, match="1-D"):
            self.drt.subspace_eval(np.zeros((2, 2), dtype=np.int32))

    def test_csf_pair_eval_2d_bra_raises(self):
        with pytest.raises(ValueError, match="1-D"):
            hexpr.csf_pair_eval(np.zeros((2, 4), dtype=np.uint8),
                                np.zeros(4, dtype=np.uint8))

    def test_csf_pair_eval_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="length"):
            hexpr.csf_pair_eval(np.zeros(4, dtype=np.uint8),
                                np.zeros(5, dtype=np.uint8))

    def test_csf_block_eval_1d_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            hexpr.csf_block_eval(np.zeros(4, dtype=np.uint8),
                                 np.zeros((1, 4), dtype=np.uint8))

    def test_csf_block_eval_norb_mismatch_raises(self):
        with pytest.raises(ValueError, match="norb"):
            hexpr.csf_block_eval(np.zeros((2, 4), dtype=np.uint8),
                                 np.zeros((2, 5), dtype=np.uint8))

    def test_csf_subspace_eval_1d_raises(self):
        with pytest.raises(ValueError, match="2-D"):
            hexpr.csf_subspace_eval(np.zeros(4, dtype=np.uint8))

    def test_invalid_which_csf_pair_eval(self):
        drt = self.drt
        steps = drt.csfs()[0]
        with pytest.raises(ValueError, match="which"):
            hexpr.csf_pair_eval(steps, steps, which="bad")
