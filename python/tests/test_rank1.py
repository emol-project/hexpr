"""Tests for hexpr.rank1_pair (reduced matrix elements of (a+_p a~_q)^1).

The values were settled against determinants built with csf_to_dets
during development; here the worked examples of the spin
formulation and two identities are asserted, as in c/tests/test_rank1.c.
"""
import itertools
from math import sqrt

import numpy as np
import pytest

import hexpr


def steps(s):
    return np.array(["eudf".index(c) for c in s], dtype=np.uint8)


def R(bra, ket):
    p, q, r = hexpr.rank1_pair(steps(bra), steps(ket))
    return {(int(a), int(b)): float(c) for a, b, c in zip(p, q, r)}


@pytest.mark.parametrize("bra,ket,pq,val", [
    ("udf", "fuu", (2, 0), sqrt(1.5)),     # (a+_3 a~_1)^1, singlet <- triplet
    ("uuu", "udu", (0, 0), sqrt(2.0)),     # spin density, quartet <- doublet
    ("ud", "uu", (0, 0), -sqrt(1.5)),
    ("ud", "uu", (1, 1), sqrt(1.5)),
])
def test_worked_examples(bra, ket, pq, val):
    assert R(bra, ket)[pq] == pytest.approx(val, abs=1e-12)


def _csfs(norb, nel):
    out = []
    for st in itertools.product("eudf", repeat=norb):
        b, ok = 0, True
        for c in st:
            b += {"u": 1, "d": -1}.get(c, 0)
            if b < 0:
                ok = False
                break
        if ok and sum({"e": 0, "u": 1, "d": 1, "f": 2}[c] for c in st) == nel:
            out.append(("".join(st), b))
    return out


@pytest.mark.parametrize("norb,nel", [(3, 2), (3, 3), (4, 4), (4, 3)])
def test_total_spin_and_adjoint(norb, nel):
    cs = _csfs(norb, nel)
    for s, b in cs:
        if b:
            S = b / 2
            tot = sum(v for (p, q), v in R(s, s).items() if p == q)
            assert tot == pytest.approx(-sqrt(2) * sqrt(S * (S + 1) * (2 * S + 1)), abs=1e-12)
    for (bra, sb), (ket, sk) in itertools.product(cs, cs):
        if abs(sb - sk) > 2:
            continue
        fwd, rev = R(bra, ket), R(ket, bra)
        sign = (-1) ** ((sb - sk) // 2)
        for (p, q), v in fwd.items():
            assert rev.get((q, p), 0.0) == pytest.approx(sign * v, abs=1e-12)
        assert len(fwd) == len(rev)


def test_rejects_bad_input():
    with pytest.raises(hexpr.InvalidArg):
        hexpr.rank1_pair(steps("du"), steps("uu"))
    with pytest.raises(ValueError):
        hexpr.rank1_pair(steps("uu"), steps("uuu"))
