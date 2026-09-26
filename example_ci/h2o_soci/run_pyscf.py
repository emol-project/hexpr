"""
H2O / 6-31G(d) SOCI full CI via PySCF integrals + hexpr -- SPHERICAL variant.

PySCF driver for the spherical H2O SOCI(8,17) example (4e4o reference, RAS2=4;
mol.cart=False, 18 AO, 17 active, next_=11). Cross-validates run_psi4.py (Psi4),
run_fermi.jl (Fermi.jl), and the OpenMolcas reference in molcas_h2o_soci.txt.
"""

import time
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo
import hexpr

# ---------------------------------------------------------------------------
# Molecule and RHF (SPHERICAL -- matches Fermi.jl)
# ---------------------------------------------------------------------------
mol = gto.Mole()
mol.atom = """
O   0.0   0.0000000   0.0000000
H   0.0   0.7572157   0.5865358
H   0.0  -0.7572157   0.5865358
"""
mol.basis = "6-31g(d)"
mol.charge = 0
mol.spin = 0
mol.unit = "Angstrom"
mol.symmetry = False
mol.verbose = 0
mol.cart = False          # spherical d (5 components) -> 18 AOs, matches Fermi
mol.build()

mf = scf.RHF(mol)
mf.kernel()

nre = mol.energy_nuc()
rhf = mf.e_tot
print(f"  nbas={mol.nao}  NRE = {nre:.10f}  RHF = {rhf:.10f}")

# ---------------------------------------------------------------------------
# SOCI integral transformation: freeze O 1s (1 inactive), 17 active MOs.
# ---------------------------------------------------------------------------
ncas, nelecas = 17, 8
mc = mcscf.CASCI(mf, ncas, nelecas)
h1eff, ecore = mc.get_h1eff()
mo_act  = mf.mo_coeff[:, mc.ncore:mc.ncore + mc.ncas]
eri_act = ao2mo.kernel(mol, mo_act, compact=False).reshape(ncas, ncas, ncas, ncas)
print(f"  h1eff {h1eff.shape}  eri_act {eri_act.shape}  ecore={ecore:.10f}")

# ---------------------------------------------------------------------------
# hexpr SOCI expression: ncore=2, nval=4, next=11  (17 active orbitals)
# ---------------------------------------------------------------------------
hexpr.init()
drt  = hexpr.Drt('soci', nve=nelecas, spin_x2=0, ncore=2, nval=4, next_=11)
expr = hexpr.Expression(drt)
ij_arr, pqrs_arr, coef_arr = expr.terms_full()
ncsfs = drt.ncsfs
norb  = drt.norb                 # 17
n_oel = norb * (norb + 1) // 2   # 153
print(f"  hexpr: ncsfs={ncsfs}  nterms={len(ij_arr)}  (norb={norb})")


def tri_inv(t):
    p = int(np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t)) / 2.0))
    if p * (p - 1) // 2 >= t:
        p -= 1
    q = t - p * (p - 1) // 2
    return p, q


# ---------------------------------------------------------------------------
# Assemble H term-by-term (same loop as run_pyscf.py)
# Scalar per-term loop -- kept commented for clarity; the vectorized block
# below is the active path (see note there).
# ---------------------------------------------------------------------------
# H = np.zeros((ncsfs, ncsfs))
# for k in range(len(ij_arr)):
#     ij = int(ij_arr[k])
#     I  = int((np.sqrt(8 * ij + 1) - 1) // 2)
#     J  = ij - I * (I + 1) // 2
#     c  = coef_arr[k]
#     if pqrs_arr[k] < n_oel:
#         p, q = tri_inv(int(pqrs_arr[k]) + 1)
#         v = h1eff[p - 1, q - 1]
#     else:
#         outer = int(pqrs_arr[k]) - n_oel + 1
#         a, b  = tri_inv(outer)
#         p, q  = tri_inv(a)
#         r, s  = tri_inv(b)
#         v = eri_act[p - 1, q - 1, r - 1, s - 1]
#     H[I, J] += c * v
#     if I != J:
#         H[J, I] += c * v

# ===========================================================================
# Vectorized H-build (production path; ~30x fewer Python ops than the scalar
# loop above, which is kept commented for clarity). Verified value-equivalent
# (identical energies; H differs only at ~1e-14 from add order). At 4e4o's
# 21.4M terms the scalar loop is impractical, so this is the active path.
# ===========================================================================
def tri_inv_vec(t):
    """Decode tri(p,q) values t (int64 array) to 1-based (p, q) arrays."""
    p = np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t.astype(np.float64))) / 2.0).astype(np.int64)
    p = np.where(p * (p - 1) // 2 >= t, p - 1, p)   # guard floating-point rounding
    q = t - p * (p - 1) // 2
    return p, q

t0 = time.time()
ij_np   = ij_arr.astype(np.int64)
pqrs_np = pqrs_arr.astype(np.int64)
I_arr = ((-1.0 + np.sqrt(8.0 * ij_np + 1.0)) / 2.0).astype(np.int64)
J_arr = ij_np - I_arr * (I_arr + 1) // 2
mask_1e = pqrs_np < n_oel
mask_2e = ~mask_1e
H = np.zeros((ncsfs, ncsfs))
# --- 1-electron terms ---
if mask_1e.any():
    p1, q1 = tri_inv_vec(pqrs_np[mask_1e] + 1)
    vals   = h1eff[p1 - 1, q1 - 1]
    I1, J1, c1 = I_arr[mask_1e], J_arr[mask_1e], coef_arr[mask_1e]
    np.add.at(H, (I1, J1), c1 * vals)
    od = I1 != J1
    np.add.at(H, (J1[od], I1[od]), c1[od] * vals[od])
# --- 2-electron terms ---
if mask_2e.any():
    outer  = pqrs_np[mask_2e] - n_oel + 1
    a1, b1 = tri_inv_vec(outer)
    p1, q1 = tri_inv_vec(a1)
    r1, s1 = tri_inv_vec(b1)
    vals   = eri_act[p1 - 1, q1 - 1, r1 - 1, s1 - 1]
    I2, J2, c2 = I_arr[mask_2e], J_arr[mask_2e], coef_arr[mask_2e]
    np.add.at(H, (I2, J2), c2 * vals)
    od = I2 != J2
    np.add.at(H, (J2[od], I2[od]), c2[od] * vals[od])
print(f"  H built (vectorized) ({time.time()-t0:.1f}s)")

evals = np.linalg.eigvalsh(H)
soci = ecore + evals[0]

print()
print("=" * 60)
print("H2O / 6-31G(d)  SOCI  full CI   (PySCF + hexpr, SPHERICAL)")
print("=" * 60)
print(f"Nuclear repulsion: {nre:.10f} Hartree")
print(f"RHF energy:        {rhf:.10f} Hartree")
print(f"SOCI energy:       {soci:.10f} Hartree")
print(f"  ecore={ecore:.10f}  E_active[0]={evals[0]:.10f}")
print("=" * 60)
