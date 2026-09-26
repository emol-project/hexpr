"""
Validation: H2O / 6-31G(d) RHF and SOCI energies vs OpenMolcas RASSCF.

All three quantities (NRE, RHF, SOCI) agree with OpenMolcas to 1e-8 Hartree.
The SOCI total energy printed by OpenMolcas has 8 digits of precision, which
sets the achievable tolerance for cross-code comparison.

OpenMolcas partition (RAS1=2, RAS2=2, RAS3=14, nActEl=8 2 2, CIONLY):
  Inactive (1 MO):  O 1s
  Active (18 MOs): MO #2-19 -- hexpr ncore=2, nval=2, next=14

mol.cart=True: OpenMolcas with 6-31G* uses Cartesian d functions (19 AOs);
spherical d gives 18 AOs and a ~1.4e-3 Hartree error on RHF.

Index decode: see docs/PQRS_INDEX_SPEC.md.
Reference:    validation/h2o_soci/README.md
Exit code:    0 (all PASS) or 1 (any FAIL).
"""

import sys
import time
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo
import hexpr

NRE_TARGET  = 9.1892096196
RHF_TARGET  = -76.0105034691
SOCI_TARGET = -76.19806885
NRE_TOL  = 1e-8
RHF_TOL  = 1e-8
SOCI_TOL = 1e-8

# ---------------------------------------------------------------------------
# Molecule and RHF
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
mol.verbose = 3
# OpenMolcas with 6-31G* uses Cartesian d functions (6 components, including
# the s-type x^2+y^2+z^2 mixture), giving 19 AOs. PySCF defaults to spherical
# d (5 components, 18 AOs), which shifts the RHF energy by ~1.4e-3 Hartree.
# cart=True matches OpenMolcas.
mol.cart = True
mol.build()

mf = scf.RHF(mol)
mf.kernel()

nre = mol.energy_nuc()
rhf = mf.e_tot

# ---------------------------------------------------------------------------
# SOCI integrals via PySCF
# ---------------------------------------------------------------------------
# OpenMolcas Inactive=1 -> PySCF ncore=1 (freeze O 1s);
# active = all remaining 18 MOs
ncas, nelecas, ncore_pyscf = 18, 8, 1
mc = mcscf.CASCI(mf, ncas, nelecas)
t0 = time.time()
h1eff, ecore = mc.get_h1eff()
mo_act  = mf.mo_coeff[:, mc.ncore:mc.ncore + mc.ncas]
eri_act = ao2mo.kernel(mol, mo_act, compact=False).reshape(ncas, ncas, ncas, ncas)
print(f"  Integrals ready ({time.time()-t0:.1f}s)  "
      f"h1eff {h1eff.shape}  eri_act {eri_act.shape}")

# ---------------------------------------------------------------------------
# hexpr Brooks expression
# hexpr partition: ncore=2 (MO2-3), nval=2 (MO4-5), next=14 (MO6-19)
# pqrs indices run over the 18-orbital active space (0-based 0..17)
# ---------------------------------------------------------------------------
hexpr.init()
drt  = hexpr.Drt('soci', nve=nelecas, spin_x2=0, ncore=2, nval=2, next_=14)
expr = hexpr.Expression(drt)
t0 = time.time()
ij_arr, pqrs_arr, coef_arr = expr.terms_full()
ncsfs  = drt.ncsfs
norb   = drt.norb    # 18
n_oel  = norb * (norb + 1) // 2   # 171
print(f"  hexpr: ncsfs={ncsfs}  nterms={len(ij_arr)}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Vectorised tri inverse (see docs/PQRS_INDEX_SPEC.md)
# ---------------------------------------------------------------------------
def tri_inv_vec(t):
    """Decode tri(p,q) values t (int64 array) to 1-based (p, q) arrays."""
    p = np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t.astype(np.float64))) / 2.0).astype(np.int64)
    p = np.where(p * (p - 1) // 2 >= t, p - 1, p)   # guard floating-point rounding
    q = t - p * (p - 1) // 2
    return p, q

# ---------------------------------------------------------------------------
# Vectorised H matrix construction
# ---------------------------------------------------------------------------
t0 = time.time()
ij_np   = ij_arr.astype(np.int64)
pqrs_np = pqrs_arr.astype(np.int64)

# Decode ij -> (I, J) with I >= J
I_arr = ((-1.0 + np.sqrt(8.0 * ij_np + 1.0)) / 2.0).astype(np.int64)
J_arr = ij_np - I_arr * (I_arr + 1) // 2

# Split 1-electron / 2-electron terms
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

print(f"  H built ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Diagonalise and compute total energy
# ---------------------------------------------------------------------------
t0 = time.time()
E_active = np.linalg.eigh(H)[0]
soci = ecore + E_active[0]
print(f"  eigh ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
nre_diff  = abs(nre  - NRE_TARGET)
rhf_diff  = abs(rhf  - RHF_TARGET)
soci_diff = abs(soci - SOCI_TARGET)
nre_ok  = nre_diff  <= NRE_TOL
rhf_ok  = rhf_diff  <= RHF_TOL
soci_ok = soci_diff <= SOCI_TOL
passed  = nre_ok and rhf_ok and soci_ok

print()
print("=" * 70)
print("H2O / 6-31G(d)  SOCI validation")
print("=" * 70)
print(f"Nuclear repulsion: {nre:.10f}  target {NRE_TARGET:.10f}"
      f"  diff {nre_diff:.2e}  {'OK' if nre_ok else 'FAIL'}")
print(f"RHF energy:        {rhf:.10f}  target {RHF_TARGET:.10f}"
      f"  diff {rhf_diff:.2e}  {'OK' if rhf_ok else 'FAIL'}")
print(f"SOCI energy:       {soci:.10f}  target {SOCI_TARGET:.10f}"
      f"  diff {soci_diff:.2e}  {'OK' if soci_ok else 'FAIL'}")
print(f"  ecore={ecore:.10f}  E_active[0]={E_active[0]:.10f}")
print(f"Status: {'PASS' if passed else 'FAIL'}")
print("=" * 70)

sys.exit(0 if passed else 1)
