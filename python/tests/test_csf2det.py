"""Tests for Drt.csf_ndets / Drt.csf_to_dets (the csf2det binding).

WHAT IS AND IS NOT ASSERTED HERE
--------------------------------
The overall sign of a CSF is a convention, not a measured quantity: it is
whatever the Eq. (5) coupling order and the block operator order produce,
and nothing inside hexpr constrains it.  <I|J> and U H_det U^T are both
invariant under flipping the sign of a whole CSF.  So these tests assert

  * the number of determinants,
  * the alpha/beta bit masks and their order,
  * coefficient magnitudes and the unit norm,
  * RELATIVE signs within one CSF,

and never the sign of an individual coefficient.  The relative signs are
the ones fixed by validation/csf2det_vs_pyscf.py against PySCF; the same
facts are asserted in C by c/tests/test_csf2det.c.

The two reference DRTs mirror that C test:

  singlet  SOCI nve=2 spin_x2=0 ncore=1 nval=1 next=2   ->  norb=4, CSF "udee"
  triplet  SOCI nve=2 spin_x2=2 ncore=1 nval=1 next=2   ->  norb=4, CSF "uuee"

For "udee" the open shells are orbitals 0 and 1, so the two determinants
are (alpha={0}, beta={1}) and (alpha={1}, beta={0}).  Determinant order is
lexicographically ascending in the set of open-shell positions carrying
alpha, so {0} comes first.
"""
import numpy as np
import pytest
import hexpr

SQRT_HALF = 1.0 / np.sqrt(2.0)


def _singlet_drt():
    return hexpr.Drt(kind='soci', nve=2, spin_x2=0,
                     ncore=1, nval=1, next_=2)


def _triplet_drt():
    return hexpr.Drt(kind='soci', nve=2, spin_x2=2,
                     ncore=1, nval=1, next_=2)


def _index_of(drt, step_string):
    """Canonical index of a CSF given as an 'eudf' step string."""
    return drt.csf_index(hexpr.step_parse(step_string, drt.norb))


# ---------------------------------------------------------------------------
# Shape, dtype and count
# ---------------------------------------------------------------------------

def test_ndets_is_exact_and_matches_to_dets():
    """csf_ndets is a count, not a bound: csf_to_dets writes exactly that many."""
    drt = _singlet_drt()
    i = _index_of(drt, "udee")
    n = drt.csf_ndets(i, 0)
    assert n == 2
    alpha, beta, coef = drt.csf_to_dets(i, 0)
    assert len(coef) == n
    assert alpha.shape[0] == n
    assert beta.shape[0] == n


def test_masks_are_uint64_and_two_dimensional():
    """Masks stay (ndets, nwords) uint64 even when nwords == 1.

    A 1-D result would silently change the caller's indexing at norb > 64,
    and a signed or float dtype would corrupt bit 63.
    """
    drt = _singlet_drt()
    i = _index_of(drt, "udee")
    alpha, beta, coef = drt.csf_to_dets(i, 0)

    nwords = (drt.norb + 63) // 64
    assert nwords == 1
    assert alpha.ndim == 2 and beta.ndim == 2
    assert alpha.shape == (2, nwords)
    assert beta.shape == (2, nwords)
    assert alpha.dtype == np.uint64
    assert beta.dtype == np.uint64
    assert coef.dtype == np.float64
    assert coef.ndim == 1


# ---------------------------------------------------------------------------
# The determinants themselves
# ---------------------------------------------------------------------------

def test_open_shell_singlet_determinants_and_order():
    """udee at Ms=0: two determinants, in ascending alpha-set order."""
    drt = _singlet_drt()
    i = _index_of(drt, "udee")
    alpha, beta, _ = drt.csf_to_dets(i, 0)

    # open shells are orbitals 0 and 1; alpha-set {0} sorts before {1}
    assert int(alpha[0, 0]) == 0b01
    assert int(beta[0, 0]) == 0b10
    assert int(alpha[1, 0]) == 0b10
    assert int(beta[1, 0]) == 0b01


def test_open_shell_singlet_magnitudes_and_norm():
    """udee at Ms=0: both magnitudes 1/sqrt(2), unit norm."""
    drt = _singlet_drt()
    i = _index_of(drt, "udee")
    _, _, coef = drt.csf_to_dets(i, 0)

    assert np.allclose(np.abs(coef), SQRT_HALF)
    assert np.isclose(np.sum(coef ** 2), 1.0)


def test_singlet_and_triplet_relative_signs_are_opposite():
    """The relative sign is the one the PySCF comparison fixed.

    Singlet: the two determinants carry the SAME sign.
    Triplet: OPPOSITE signs.  Neither individual sign is asserted -- the
    overall sign of a CSF is a convention.
    """
    s = _singlet_drt()
    si = _index_of(s, "udee")
    _, _, s_coef = s.csf_to_dets(si, 0)

    t = _triplet_drt()
    ti = _index_of(t, "uuee")
    _, _, t_coef = t.csf_to_dets(ti, 0)

    assert np.allclose(np.abs(t_coef), SQRT_HALF)
    assert s_coef[0] * s_coef[1] > 0.0, "singlet: same relative sign"
    assert t_coef[0] * t_coef[1] < 0.0, "triplet: opposite relative sign"


def test_every_csf_is_normalized_at_every_ms():
    """sum(coef**2) == 1 for every CSF of the triplet DRT at every Ms."""
    drt = _triplet_drt()
    for ms_x2 in (-2, 0, 2):
        for i in range(drt.ncsfs):
            _, _, coef = drt.csf_to_dets(i, ms_x2)
            assert np.isclose(np.sum(coef ** 2), 1.0), \
                f"csf {i} at ms_x2={ms_x2} is not normalized"


# ---------------------------------------------------------------------------
# Multi-word determinants
# ---------------------------------------------------------------------------

_STEP_OCC = {0: 0, 1: 1, 2: 1, 3: 2}   # E, U, D, F -> electrons


def _occupation_from_masks(alpha_row, beta_row, norb):
    """Electrons per spatial orbital, read back out of the bit masks."""
    occ = []
    for p in range(norb):
        w, b = p >> 6, p & 63
        occ.append(((int(alpha_row[w]) >> b) & 1)
                   + ((int(beta_row[w]) >> b) & 1))
    return occ


def test_wide_determinants_span_two_words():
    """norb > 64: masks stay (ndets, 2) and bits land in the right word.

    Everything else in this file runs at norb <= 6, where nwords == 1 and
    the second index is never exercised.  This case walks every CSF of a
    70-orbital DRT and reconstructs the occupation from the masks, so a
    wrong word index or a wrong stride shows up as a mismatch rather than
    as a silently different answer.
    """
    drt = hexpr.Drt(kind='foci', nve=2, spin_x2=0,
                    ncore=1, nval=1, next_=68)
    assert drt.norb == 70
    nwords = (drt.norb + 63) // 64
    assert nwords == 2

    saw_high_word = False
    for i in range(drt.ncsfs):
        alpha, beta, coef = drt.csf_to_dets(i, 0)
        assert alpha.shape == (len(coef), nwords)
        assert beta.shape == (len(coef), nwords)
        assert np.isclose(np.sum(coef ** 2), 1.0)

        want = [_STEP_OCC[s] for s in drt.csf_steps(i).tolist()]
        for k in range(len(coef)):
            got = _occupation_from_masks(alpha[k], beta[k], drt.norb)
            assert got == want, f"csf {i} det {k}: occupation mismatch"
            if int(alpha[k, 1]) or int(beta[k, 1]):
                saw_high_word = True

    assert saw_high_word, \
        "no determinant occupied an orbital >= 64; the wide path was not taken"


# ---------------------------------------------------------------------------
# ms_x2 is a parameter, not an assumption
# ---------------------------------------------------------------------------

def test_negative_ms_is_accepted():
    """Ms is not assumed to equal S; negative Ms is a legal argument.

    At Ms = -S the triplet has a single determinant (both electrons beta).
    """
    drt = _triplet_drt()
    i = _index_of(drt, "uuee")
    alpha, beta, coef = drt.csf_to_dets(i, -2)

    assert len(coef) == 1
    assert np.isclose(abs(coef[0]), 1.0)
    assert int(alpha[0, 0]) == 0
    assert int(beta[0, 0]) == 0b11


def test_ms_x2_of_wrong_parity_raises():
    """2*Ms must have the same parity as spin_x2."""
    drt = _singlet_drt()
    i = _index_of(drt, "udee")
    with pytest.raises(hexpr.HexprError):
        drt.csf_ndets(i, 1)


def test_ms_x2_beyond_spin_raises():
    """|Ms| <= S is required."""
    drt = _singlet_drt()
    i = _index_of(drt, "udee")
    with pytest.raises(hexpr.HexprError):
        drt.csf_ndets(i, 4)


# ---------------------------------------------------------------------------
# Index handling
# ---------------------------------------------------------------------------

def test_index_out_of_range_raises():
    drt = _singlet_drt()
    with pytest.raises(hexpr.HexprError):
        drt.csf_ndets(drt.ncsfs, 0)
    with pytest.raises(hexpr.HexprError):
        drt.csf_ndets(-1, 0)


def test_csfset_index_is_the_input_position():
    """For a CSFSET sub-DRT, index is the INPUT position, as for csf_steps.

    The rows below are the explicit CSF set used by
    validation/csf2det_vs_pyscf.py case (a): fe, ef, ud on two orbitals.
    Input row 2 is the open-shell singlet, so it is the one with two
    determinants -- and rows 0 and 1 are closed shells with one each.
    """
    rows = np.array([[3, 0], [0, 3], [1, 2]], dtype=np.uint8)
    drt = hexpr.Drt.from_csfset(rows, norb=2)
    assert drt.ncsfs == 3

    assert drt.csf_ndets(0, 0) == 1
    assert drt.csf_ndets(1, 0) == 1
    assert drt.csf_ndets(2, 0) == 2

    # row 0 is "fe": orbital 0 doubly occupied, one determinant
    alpha, beta, coef = drt.csf_to_dets(0, 0)
    assert int(alpha[0, 0]) == 0b01
    assert int(beta[0, 0]) == 0b01
    assert np.isclose(abs(coef[0]), 1.0)
