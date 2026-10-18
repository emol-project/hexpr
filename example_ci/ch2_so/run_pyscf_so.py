"""
CH2 / cc-pVDZ, CAS(6,6): spin-orbit matrix from PySCF integrals + hexpr.

The one-electron Breit-Pauli spin-orbit operator
    H_SO = (alpha^2/2) sum_A Z_A (r_A x p)/r_A^3 . s  =  sum_pq sum_m (-1)^m h_{pq,-m} s_m(pq)
is assembled between the CSFs of the singlet and triplet CAS spaces in two
independent ways, for every magnetic component:

  (a) hexpr: reduced matrix elements R_pq = <S'||(a+_p a~_q)^1||S> from
      hexpr.rank1_pair, and
        <S'M'|H_SO|SM> = -(1/sqrt2) (-1)^(S'-M') sum_m (-1)^m (S' 1 S; -M' m M)
                         sum_pq h_{pq,-m} R_pq;
  (b) determinants: H_SO built from spin-orbital operators in the
      determinant basis, and projected on the CSFs expanded by csf_to_dets.

Spatial integrals: PySCF's int1e_pnucxp is (p* | nuc x p |) = int V (grad i x
grad j) with V the nuclear potential energy, so the Hermitian vector operator
is h_pq = (alpha^2/2) * 1j * int1e_pnucxp. (For effective charges:
int1e_pnucxp = -sum_A Z_A int1e_prinvxp(A); replace Z_A by Z_eff.)

The one-electron operator needs no core contribution here: the core is doubly
occupied in every CSF and a CASCI has no core excitation.
Run with the in-tree library:  python run_pyscf_so.py
"""
import sys, os, time
from math import factorial, sqrt
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, fci

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
import hexpr

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--gamess", help="GAMESS SO-CI log to compare with")
ap.add_argument("--zeff", help="effective charges, e.g. C=3.6,H=1.0 (default: true charges)")
args = ap.parse_args()

ALPHA = 1 / 137.035999084
HARTREE_CM = 219474.6313632

# ---------------------------------------------------------------- 3j symbol
def _f(n):
    return factorial(n)

def threej(j1, m1, j2, m2, j3, m3):
    """Wigner 3j, all arguments doubled (2j, 2m)."""
    if m1 + m2 + m3 != 0 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2 or (j1 + j2 + j3) % 2:
        return 0.0
    a = [(j1 + j2 - j3) // 2, (j1 - j2 + j3) // 2, (-j1 + j2 + j3) // 2]
    tri = _f(a[0]) * _f(a[1]) * _f(a[2]) / _f((j1 + j2 + j3) // 2 + 1)
    pre = sqrt(tri * _f((j1 + m1) // 2) * _f((j1 - m1) // 2) * _f((j2 + m2) // 2)
               * _f((j2 - m2) // 2) * _f((j3 + m3) // 2) * _f((j3 - m3) // 2))
    s = 0.0
    for k in range(0, j1 + j2 + j3 + 2):
        d = [k, (j1 + j2 - j3) // 2 - k, (j1 - m1) // 2 - k, (j2 + m2) // 2 - k,
             (j3 - j2 + m1) // 2 + k, (j3 - j1 - m2) // 2 + k]
        if min(d) < 0:
            continue
        den = 1
        for x in d:
            den *= _f(x)
        s += (-1) ** k / den
    return (-1) ** ((j1 - j2 - m3) // 2) * pre * s

# ---------------------------------------------------------------- molecule
# Bohr, so that GAMESS (gamess/ch2_*.inp) sees the identical geometry: the two
# programs convert Angstrom with slightly different constants (1.8e-8 Eh in RHF).
mol = gto.M(atom="""
C   0.0000000000   0.0000000000   0.1899174755
H   0.0000000000   1.8721516716  -1.1400717710
H   0.0000000000  -1.8721516716  -1.1400717710
""", basis="cc-pvdz", unit="Bohr", verbose=0)
mf = scf.RHF(mol)
mf.conv_tol = 1e-12
mf.kernel()
print(f"RHF {mf.e_tot:.10f}")
ncas, nelecas = 6, 6
mc = mcscf.CASCI(mf, ncas, nelecas)
h1eff, ecore = mc.get_h1eff()
C_act = mf.mo_coeff[:, mc.ncore:mc.ncore + ncas]
eri = ao2mo.kernel(mol, C_act, compact=False).reshape((ncas,) * 4)

A_ao = mol.intor("int1e_pnucxp", comp=3)                       # real, antisymmetric
if args.zeff:                                                  # -sum_A Zeff_A prinvxp(A)
    zmap = {k: float(v) for k, v in (kv.split("=") for kv in args.zeff.split(","))}
    zeff = [zmap[mol.atom_pure_symbol(a)] for a in range(mol.natm)]
    A_z = np.zeros_like(A_ao)
    for a in range(mol.natm):
        with mol.with_rinv_origin(mol.atom_coord(a)):
            A_z -= zeff[a] * mol.intor("int1e_prinvxp", comp=3)
    print("Z_eff:", dict(zip([mol.atom_pure_symbol(a) for a in range(mol.natm)], zeff)))
    A_ao = A_z
A = np.einsum("ip,cij,jq->cpq", C_act, A_ao, C_act)
h_cart = 0.5 * ALPHA ** 2 * 1j * A                             # Hermitian (3, n, n)
assert np.abs(h_cart - h_cart.conj().transpose(0, 2, 1)).max() < 1e-14
h_sph = {+1: -(h_cart[0] + 1j * h_cart[1]) / sqrt(2),         # spherical components
          0: h_cart[2],
         -1: (h_cart[0] - 1j * h_cart[1]) / sqrt(2)}

# ---------------------------------------------------------------- spin-free CAS
hexpr.init()
n_oel = ncas * (ncas + 1) // 2

def tri_inv(t):
    p = np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t.astype(np.float64))) / 2.0).astype(np.int64)
    p = np.where(p * (p - 1) // 2 >= t, p - 1, p)
    return p, t - p * (p - 1) // 2

def cas_space(s2):
    drt = hexpr.Drt("soci", nve=nelecas, spin_x2=s2, ncore=0, nval=ncas, next_=0)
    ij, pqrs, coef = hexpr.Expression(drt).terms_full()
    ij, pqrs = ij.astype(np.int64), pqrs.astype(np.int64)
    I = ((-1.0 + np.sqrt(8.0 * ij + 1.0)) / 2.0).astype(np.int64)
    J = ij - I * (I + 1) // 2
    one = pqrs < n_oel
    v = np.empty(len(pqrs))
    p, q = tri_inv(pqrs[one] + 1)
    v[one] = h1eff[p - 1, q - 1]
    a, b = tri_inv(pqrs[~one] - n_oel + 1)
    p, q = tri_inv(a); r, s = tri_inv(b)
    v[~one] = eri[p - 1, q - 1, r - 1, s - 1]
    H = np.zeros((drt.ncsfs, drt.ncsfs))
    np.add.at(H, (I, J), coef * v)
    od = I != J
    np.add.at(H, (J[od], I[od]), coef[od] * v[od])
    e, c = np.linalg.eigh(H)
    return drt, drt.csfs(), e, c

space = {s2: cas_space(s2) for s2 in (0, 2)}
for s2, (drt, steps, e, c) in space.items():
    print(f"2S={s2}: {drt.ncsfs} CSFs, E0 = {ecore + e[0]:.10f}")
# PySCF CASCI for the same spaces (sanity of the spin-free part)
e_s = mcscf.CASCI(mf, ncas, (3, 3)); e_s.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0)
e_t = mcscf.CASCI(mf, ncas, (4, 2))
print(f"PySCF CASCI singlet {e_s.kernel()[0]:.10f}, triplet {e_t.kernel()[0]:.10f}")

# ---------------------------------------------------------------- (a) hexpr
def R_tensor(bra_steps, ket_steps):
    R = np.zeros((len(bra_steps), len(ket_steps), ncas, ncas))
    for i, b in enumerate(bra_steps):
        for j, k in enumerate(ket_steps):
            p, q, r = hexpr.rank1_pair(b, k)
            R[i, j, p, q] = r
    return R

def so_csf_hexpr(sb, sk):
    """{(M', M): <I S' M'|H_SO|J S M>} for all CSFs (2M doubled keys)."""
    R = R_tensor(space[sb][1], space[sk][1])
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

# ---------------------------------------------------------------- (b) determinants
dets = [d for d in range(1 << (2 * ncas)) if bin(d).count("1") == nelecas]
dindex = {d: i for i, d in enumerate(dets)}

def apply_cre(d, i):
    if d >> i & 1:
        return None, 0
    return d | (1 << i), (-1) ** bin(d & ((1 << i) - 1)).count("1")

def apply_ann(d, i):
    if not d >> i & 1:
        return None, 0
    return d ^ (1 << i), (-1) ** bin(d & ((1 << i) - 1)).count("1")

s_mat = {  # (sigma, tau) with sigma, tau = 0 (alpha), 1 (beta)
    0: np.array([[0, 0.5], [0.5, 0]], dtype=complex),
    1: np.array([[0, -0.5j], [0.5j, 0]], dtype=complex),
    2: np.array([[0.5, 0], [0, -0.5]], dtype=complex),
}
coef_pq = np.einsum("cpq,cst->pqst", h_cart, np.array([s_mat[c] for c in range(3)]))
H_det = np.zeros((len(dets), len(dets)), dtype=complex)
for jd, d in enumerate(dets):
    for p in range(ncas):
        for q in range(ncas):
            for s in range(2):
                for t in range(2):
                    c = coef_pq[p, q, s, t]
                    if c == 0:
                        continue
                    d1, s1 = apply_ann(d, q + t * ncas)
                    if d1 is None:
                        continue
                    d2, s2_ = apply_cre(d1, p + s * ncas)
                    if d2 is None:
                        continue
                    H_det[dindex[d2], jd] += c * s1 * s2_
assert np.abs(H_det - H_det.conj().T).max() < 1e-14

def csf_matrix(s2, ms2):
    drt = space[s2][0]
    U = np.zeros((len(dets), drt.ncsfs))
    for i in range(drt.ncsfs):
        al, be, co = drt.csf_to_dets(i, ms2)
        for a, b, x in zip(al[:, 0], be[:, 0], co):
            U[dindex[int(a) | (int(b) << ncas)], i] = x
    return U

U = {(s2, ms): csf_matrix(s2, ms) for s2 in (0, 2) for ms in range(-s2, s2 + 1, 2)}

# ---------------------------------------------------------------- compare
t0 = time.time()
worst, scale = 0.0, 0.0
blocks = {}
for sb in (0, 2):
    for sk in (0, 2):
        if sb + sk < 2:
            continue
        X = so_csf_hexpr(sb, sk)
        blocks[sb, sk] = X
        for (mb, mk), Xa in X.items():
            Xb = U[sb, mb].T @ H_det @ U[sk, mk]
            worst = max(worst, np.abs(Xa - Xb).max())
            scale = max(scale, np.abs(Xb).max())
print(f"\nCSF spin-orbit matrix, hexpr vs determinants, all M: "
      f"max |diff| = {worst:.2e} Eh  (largest element {scale:.3e} Eh)  [{time.time()-t0:.1f}s]")

# ---------------------------------------------------------------- states
nroot = 3
basis = []                      # (2S, root, 2M)
for s2 in (0, 2):
    for n in range(nroot):
        for ms in range(-s2, s2 + 1, 2):
            basis.append((s2, n, ms))
Hsi = np.zeros((len(basis), len(basis)), dtype=complex)
for a, (sb, nb, mb) in enumerate(basis):
    Hsi[a, a] += ecore + space[sb][2][nb]
    for b, (sk, nk, mk) in enumerate(basis):
        if sb + sk < 2:
            continue
        X = blocks[sb, sk][mb, mk]
        Hsi[a, b] += space[sb][3][:, nb] @ X @ space[sk][3][:, nk]
assert np.abs(Hsi - Hsi.conj().T).max() < 1e-12
print("\nsinglet-triplet coupling sqrt(sum_M |<S n|H_SO|T k,M>|^2) (cm-1), rows S n, columns T k:")
for nb in range(nroot):
    row = []
    for nk in range(nroot):
        a = basis.index((0, nb, 0))
        tot = sum(abs(Hsi[a, basis.index((2, nk, ms))]) ** 2 for ms in (-2, 0, 2))
        row.append(f"{sqrt(tot) * HARTREE_CM:10.4f}")
    print("  S%d: " % nb + " ".join(row))
lev = np.linalg.eigvalsh(Hsi)
e0 = min(ecore + space[s][2][0] for s in (0, 2))
print("\nspin-free states (cm-1 above the lowest):")
for s2 in (0, 2):
    print("  2S=%d: " % s2 + ", ".join(f"{(ecore + e - e0) * HARTREE_CM:10.2f}"
                                      for e in space[s2][2][:nroot]))
print("spin-orbit levels (cm-1 above the lowest):")
print("  " + ", ".join(f"{(x - lev[0]) * HARTREE_CM:.4f}" for x in lev))

# ---------------------------------------------------------------- GAMESS
# Compare with a GAMESS SO-CI log (RUNTYP=TRANSITN, OPERAT=HSO1; gamess/).
# GAMESS states are matched to the exact roots by energy (its C1 Davidson can
# skip roots of other irreps of C2v), and every printed <a,M'|HSO|b,M> is
# compared in magnitude: the overall phase of each state is arbitrary.
if args.gamess:
    import re
    txt = open(args.gamess).read()
    blocks_e, cur = [], []
    for n, e in re.findall(r"STATE #\s+(\d+)\s+ENERGY =\s+(-?\d+\.\d+)", txt):
        if int(n) == 1 and cur:
            blocks_e.append(cur); cur = []
        cur.append(float(e))
    blocks_e.append(cur)
    mult_of_drt = {1: 0, 3: 2}           # GAMESS label m-mA -> 2S
    roots = {}                          # (2S, gamess root) -> my root index
    for s2, el in zip((0, 2), blocks_e):
        for i, e in enumerate(el):
            d = np.abs(ecore + space[s2][2] - e)
            k = int(np.argmin(d))
            if d[k] < 1e-7:
                roots[s2, i + 1] = k
    print(f"\nGAMESS: energies of {len(roots)} states matched to exact roots (|dE| < 1e-7):")
    for s2 in (0, 2):
        print(f"  2S={s2}: GAMESS root -> exact root "
              + ", ".join(f"{g}->{k + 1}" for (s, g), k in sorted(roots.items()) if s == s2))
    pat = re.compile(r"<\s*(\d+)-(\d)A\s*\(MS=\s*(-?[\d.]+)\)\|HSO\|\s*(\d+)-(\d)A\s*"
                     r"\(MS=\s*(-?[\d.]+)\)>\s*=\s*(-?[\d.]+)\s+(-?[\d.]+)I")
    worst, n, ratios = 0.0, 0, []
    for g in pat.findall(txt):
        rb, mb_, msb, rk, mk_, msk, re_, im_ = g
        sb, sk = mult_of_drt[int(mb_)], mult_of_drt[int(mk_)]
        if (sb, int(rb)) not in roots or (sk, int(rk)) not in roots or sb + sk < 2:
            continue
        kb, kk = roots[sb, int(rb)], roots[sk, int(rk)]
        Mb, Mk = int(round(2 * float(msb))), int(round(2 * float(msk)))
        mine = space[sb][3][:, kb] @ blocks[sb, sk][Mb, Mk] @ space[sk][3][:, kk]
        gv = complex(float(re_), float(im_))
        d = abs(abs(mine) * HARTREE_CM - abs(gv))
        if d > worst:
            worst, at = d, (g[:6], abs(gv))
        if abs(gv) > 5.0:   # 4 printed decimals: ratio good to 1e-5 or better
            ratios.append(abs(gv) / (abs(mine) * HARTREE_CM))
        n += 1
    print(f"GAMESS HSO elements compared: {n}, max | |hexpr| - |GAMESS| | = {worst:.1e} cm-1"
          " (GAMESS prints 4 decimals)")
    if n:
        print(f"  largest difference at <{at[0][0]}-{at[0][1]}A MS={at[0][2]}|HSO|"
              f"{at[0][3]}-{at[0][4]}A MS={at[0][5]}>, |GAMESS| = {at[1]:.4f} cm-1")
        r = np.array(ratios)
        print(f"  |GAMESS|/|hexpr| over {len(r)} elements > 5 cm-1: mean {r.mean():.8f}, "
              f"min {r.min():.8f}, max {r.max():.8f}")
        # GAMESS's hartree -> cm-1 factor, measured on its own spin-orbit level
        # list against our levels (same Hamiltonian, so only the unit differs)
        g_lev = [float(x) for x in re.findall(r"RELATIVE E=\s+(-?\d+\.\d+) CM-1", txt)]
        g_lev = np.array(g_lev) - min(g_lev)
        my_basis = [(s2, k, ms) for s2 in (0, 2) for (s, g), k in sorted(roots.items())
                    if s == s2 for ms in range(-s2, s2 + 1, 2)]
        Hm = np.zeros((len(my_basis), len(my_basis)), dtype=complex)
        for a, (sb, kb, mb) in enumerate(my_basis):
            Hm[a, a] = ecore + space[sb][2][kb]
            for b, (sk, kk, mk) in enumerate(my_basis):
                if sb + sk >= 2:
                    Hm[a, b] += space[sb][3][:, kb] @ blocks[sb, sk][mb, mk] @ space[sk][3][:, kk]
        my_lev = np.linalg.eigvalsh(Hm)
        my_lev = (my_lev - my_lev[0]) * HARTREE_CM
        sel = my_lev > 10000.0
        f = np.mean(g_lev[sel] / my_lev[sel])
        print(f"  GAMESS hartree->cm-1 relative to CODATA 2018 ({HARTREE_CM}), from its level "
              f"list: {f:.8f} (spread {np.ptp(g_lev[sel] / my_lev[sel]):.1e}); "
              f"levels agree to {np.abs(g_lev - f * my_lev).max():.1e} cm-1 after it")
        worst_f = 0.0
        for g in pat.findall(txt):
            rb, mb_, msb, rk, mk_, msk, re_, im_ = g
            sb, sk = mult_of_drt[int(mb_)], mult_of_drt[int(mk_)]
            if (sb, int(rb)) not in roots or (sk, int(rk)) not in roots or sb + sk < 2:
                continue
            mine = (space[sb][3][:, roots[sb, int(rb)]]
                    @ blocks[sb, sk][int(round(2 * float(msb))), int(round(2 * float(msk)))]
                    @ space[sk][3][:, roots[sk, int(rk)]])
            worst_f = max(worst_f, abs(f * abs(mine) * HARTREE_CM
                                       - abs(complex(float(re_), float(im_)))))
        print(f"  with that unit: max | f|hexpr| - |GAMESS| | = {worst_f:.1e} cm-1 "
              f"over all {n} elements (print rounding 5e-5)")
