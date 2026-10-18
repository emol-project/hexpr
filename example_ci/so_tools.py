"""Shared pieces of the spin-orbit examples (ch2_so, oh_so).

- threej:           Wigner 3j symbol, doubled arguments
- CasSpace:         full CAS of one total spin through hexpr (spin-free H, roots)
- so_blocks:        CSF spin-orbit blocks <I S' M'|H_SO|J S M> from hexpr.rank1_pair
- DetSpace:         determinant space of the CAS, CSF -> determinant matrices,
                    one-electron operators in the determinant basis (the
                    independent check, and L_z for naming Omega)
- so_integrals:     Hermitian h_pq (3, n, n) of the one-electron Breit-Pauli
                    operator from PySCF, true or effective nuclear charges

Conventions (see ch2_so/README.md): PySCF int1e_pnucxp is (p* | nuc x p |),
so h = (alpha^2/2) * 1j * int1e_pnucxp; with effective charges int1e_pnucxp is
replaced by -sum_A Z_eff(A) int1e_prinvxp(A).  H_SO = sum_pq sum_m (-1)^m
h_{pq,-m} s_m(pq), s_m(pq) = -T^1_m(pq)/sqrt2, T^1(pq) = (a+_p a~_q)^1.
"""
import os, sys
from math import factorial, sqrt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import hexpr

ALPHA = 1 / 137.035999084
HARTREE_CM = 219474.6313632


def threej(j1, m1, j2, m2, j3, m3):
    """Wigner 3j, all arguments doubled (2j, 2m)."""
    if m1 + m2 + m3 != 0 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2 or (j1 + j2 + j3) % 2:
        return 0.0
    f = factorial
    tri = (f((j1 + j2 - j3) // 2) * f((j1 - j2 + j3) // 2) * f((-j1 + j2 + j3) // 2)
           / f((j1 + j2 + j3) // 2 + 1))
    pre = sqrt(tri * f((j1 + m1) // 2) * f((j1 - m1) // 2) * f((j2 + m2) // 2)
               * f((j2 - m2) // 2) * f((j3 + m3) // 2) * f((j3 - m3) // 2))
    s = 0.0
    for k in range(0, j1 + j2 + j3 + 2):
        d = [k, (j1 + j2 - j3) // 2 - k, (j1 - m1) // 2 - k, (j2 + m2) // 2 - k,
             (j3 - j2 + m1) // 2 + k, (j3 - j1 - m2) // 2 + k]
        if min(d) < 0:
            continue
        den = 1
        for x in d:
            den *= f(x)
        s += (-1) ** k / den
    return (-1) ** ((j1 - j2 - m3) // 2) * pre * s


def somf_2e(mol, D):
    """Two-electron spin-orbit mean field in the AO basis (real part, as int1e_pnucxp).

    int2e_p1vxp1 is (p* i, cross p j | k l) = int (grad i x grad j)(1) (1/r12)
    k(2) l(2): the same form as int1e_prinvxp with 1/r12 in place of 1/r, so
    the spin-same-orbit operator is (alpha^2/2) 1j int2e_p1vxp1 with a plus sign
    (the nuclear term has -Z). The mean field of Hess et al., with the
    spin-other-orbit term in the -3/2 exchange factors, for a spin-free total
    density D:
        sum_kl D_kl [ g(pq|kl) - 3/2 g(pk|lq) - 3/2 g(lq|pk) ].
    For a closed-shell core density this is exact for the core-active part of
    the Breit-Pauli two-electron operator (GAMESS OPERAT=HSO2P)."""
    G = mol.intor("int2e_p1vxp1", comp=3)                 # (3, i, j, k, l), (ij|kl)
    J = np.einsum("cpqkl,kl->cpq", G, D)
    K1 = np.einsum("cpklq,kl->cpq", G, D)                 # g(pk|lq)
    K2 = np.einsum("clqpk,kl->cpq", G, D)                 # g(lq|pk)
    return J - 1.5 * K1 - 1.5 * K2


def so_integrals(mol, C_act, zeff=None, somf_density=None):
    """Hermitian (3, n, n) spin-orbit integrals in the active MOs.
    zeff: {symbol: Z} effective charges for the one-electron term;
    somf_density: AO total density for the two-electron mean field (somf_2e)."""
    if zeff is None:
        A_ao = mol.intor("int1e_pnucxp", comp=3)
    else:
        A_ao = 0.0
        for a in range(mol.natm):
            with mol.with_rinv_origin(mol.atom_coord(a)):
                A_ao = A_ao - zeff[mol.atom_pure_symbol(a)] * mol.intor("int1e_prinvxp", comp=3)
    if somf_density is not None:
        A_ao = A_ao + somf_2e(mol, somf_density)
    A = np.einsum("ip,cij,jq->cpq", C_act, A_ao, C_act)
    h = 0.5 * ALPHA ** 2 * 1j * A
    assert np.abs(h - h.conj().transpose(0, 2, 1)).max() < 1e-14
    return h


def spherical(h_cart):
    return {+1: -(h_cart[0] + 1j * h_cart[1]) / sqrt(2), 0: h_cart[2],
            -1: (h_cart[0] - 1j * h_cart[1]) / sqrt(2)}


def _tri_inv(t):
    p = np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t.astype(np.float64))) / 2.0).astype(np.int64)
    p = np.where(p * (p - 1) // 2 >= t, p - 1, p)
    return p, t - p * (p - 1) // 2


class CasSpace:
    """Full CAS of total spin s2/2 from hexpr; .e (roots, active energy), .c (columns)."""

    def __init__(self, ncas, nelec, s2, h1, eri):
        self.s2 = s2
        self.drt = hexpr.Drt("soci", nve=nelec, spin_x2=s2, ncore=0, nval=ncas, next_=0)
        self.steps = self.drt.csfs()
        n_oel = ncas * (ncas + 1) // 2
        ij, pqrs, coef = hexpr.Expression(self.drt).terms_full()
        ij, pqrs = ij.astype(np.int64), pqrs.astype(np.int64)
        I = ((-1.0 + np.sqrt(8.0 * ij + 1.0)) / 2.0).astype(np.int64)
        J = ij - I * (I + 1) // 2
        one = pqrs < n_oel
        v = np.empty(len(pqrs))
        p, q = _tri_inv(pqrs[one] + 1)
        v[one] = h1[p - 1, q - 1]
        a, b = _tri_inv(pqrs[~one] - n_oel + 1)
        p, q = _tri_inv(a)
        r, s = _tri_inv(b)
        v[~one] = eri[p - 1, q - 1, r - 1, s - 1]
        H = np.zeros((self.drt.ncsfs, self.drt.ncsfs))
        np.add.at(H, (I, J), coef * v)
        od = I != J
        np.add.at(H, (J[od], I[od]), coef[od] * v[od])
        self.e, self.c = np.linalg.eigh(H)


def so_blocks(bra, ket, h_sph, ncas):
    """{(2M', 2M): <I S' M'|H_SO|J S M>} over all CSFs of two CasSpace objects."""
    sb, sk = bra.s2, ket.s2
    R = np.zeros((len(bra.steps), len(ket.steps), ncas, ncas))
    for i, b in enumerate(bra.steps):
        for j, k in enumerate(ket.steps):
            p, q, r = hexpr.rank1_pair(b, k)
            R[i, j, p, q] = r
    hR = {m: np.einsum("ijpq,pq->ij", R, h_sph[-m]) for m in (-1, 0, 1)}
    out = {}
    for mb in range(-sb, sb + 1, 2):
        for mk in range(-sk, sk + 1, 2):
            X = np.zeros(R.shape[:2], dtype=complex)
            for m in (-1, 0, 1):
                t = threej(sb, -mb, 2, 2 * m, sk, mk)
                if t:
                    X += (-1) ** m * t * hR[m]
            out[mb, mk] = -(1 / sqrt(2)) * (-1) ** ((sb - mb) // 2) * X
    return out


class DetSpace:
    """All determinants of nelec electrons in ncas spatial orbitals; bit p is
    alpha p, bit ncas + p is beta p, ascending operator order (as csf_to_dets)."""

    def __init__(self, ncas, nelec):
        self.n = ncas
        self.dets = [d for d in range(1 << (2 * ncas)) if bin(d).count("1") == nelec]
        self.index = {d: i for i, d in enumerate(self.dets)}

    @staticmethod
    def _sign(d, i):
        return (-1) ** bin(d & ((1 << i) - 1)).count("1")

    def one_body(self, coef):
        """sum_{p q s t} coef[p, q, s, t] a+_{p s} a_{q t}; s, t = 0 alpha, 1 beta."""
        n = self.n
        M = np.zeros((len(self.dets), len(self.dets)), dtype=complex)
        for jd, d in enumerate(self.dets):
            for p, q, s, t in zip(*np.nonzero(coef)):
                i_q, i_p = q + t * n, p + s * n
                if not d >> i_q & 1:
                    continue
                d1 = d ^ (1 << i_q)
                if d1 >> i_p & 1:
                    continue
                d2 = d1 | (1 << i_p)
                M[self.index[d2], jd] += coef[p, q, s, t] * self._sign(d, i_q) * self._sign(d1, i_p)
        return M

    def h_so(self, h_cart):
        s = np.array([[[0, .5], [.5, 0]], [[0, -.5j], [.5j, 0]], [[.5, 0], [0, -.5]]])
        return self.one_body(np.einsum("cpq,cst->pqst", h_cart, s))

    def spin_free(self, o):
        """sum_pq o_pq E_pq for a spatial one-electron matrix o."""
        return self.one_body(np.einsum("pq,st->pqst", o, np.eye(2)))

    def csf_matrix(self, space, ms2):
        U = np.zeros((len(self.dets), space.drt.ncsfs))
        for i in range(space.drt.ncsfs):
            al, be, co = space.drt.csf_to_dets(i, ms2)
            for a, b, x in zip(al[:, 0], be[:, 0], co):
                U[self.index[int(a) | (int(b) << self.n)], i] = x
        return U
