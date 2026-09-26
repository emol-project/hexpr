"""
H2O / 6-31G(d) SOCI full CI via Psi4 integrals + hexpr -- SPHERICAL variant.

Psi4 integral backend for the spherical H2O SOCI(8,17) example (4e4o reference,
RAS2=4) -- the Psi4 sibling of run_pyscf.py; only the integral-generation block
differs. The hexpr evaluation, Hamiltonian build, diagonalization, and CI-vector
report are shared with run_pyscf.py. Freezing the O 1s core leaves 17 active
orbitals, so the hexpr SOCI partition uses next_=11 (4-orbital reference RAS2).

Two Psi4 details:
  - 6-31G* with spherical d (puream True) -> 18 AOs (Pople sets default to
    Cartesian in Psi4, so puream=True is required to force 5D).
  - scf_type pk (exact ERIs, no density fitting).
Psi4 has no one-call active-space effective Hamiltonian, so ecore / h1eff
are built by hand following PySCF's h1e_for_cas (freeze O 1s):
    core_dm = 2 C_core C_core^T
    corevhf = J(core_dm) - 1/2 K(core_dm)
    ecore   = E_nuc + tr(core_dm . hcore) + 1/2 tr(core_dm . corevhf)
    h1eff   = C_act^T (hcore + corevhf) C_act

Reproduces the spherical SOCI = -76.19711710 Ha (4e4o reference; matches PySCF /
Fermi.jl / OpenMolcas to ~1e-8).

Requires Psi4. Note: psi4 1.9.1 needs libint2.so.2; install with libint
pinned below 2.10 (e.g. `conda install "psi4=1.9.1" "libint<2.10"` ->
libint 2.9.0). See example_ci/README.
"""

import os
import time
import numpy as np
import psi4
import hexpr

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Molecule and RHF (Psi4)
# Geometry copied verbatim from run_pyscf.py. no_com / no_reorient keep Psi4
# in the same frame PySCF uses (PySCF does not recenter/reorient at
# symmetry=False). 6-31G* + puream True = spherical d (5 comp) -> nbf 18,
# matching OpenMolcas's "Spherical d" card and run_pyscf.py's
# mol.cart=False.
# ---------------------------------------------------------------------------
psi4.core.set_output_file(os.path.join(SCRIPT_DIR, "psi4_h2o_soci.out"), False)
mol = psi4.geometry("""
O   0.0   0.0000000   0.0000000
H   0.0   0.7572157   0.5865358
H   0.0  -0.7572157   0.5865358
units angstrom
symmetry c1
no_com
no_reorient
""")
psi4.set_options({
    "basis": "6-31G*",
    "puream": True,            # spherical d (5 components) -> nbf 18
    "scf_type": "pk",          # exact ERIs (no density fitting)
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
})
rhf, wfn = psi4.energy("scf", return_wfn=True)
nre = mol.nuclear_repulsion_energy()
assert wfn.basisset().nbf() == 18, \
    f"expected nbf=18 (spherical d), got {wfn.basisset().nbf()}"
print(f"  NRE = {nre:.10f}  RHF = {rhf:.10f}")

# ---------------------------------------------------------------------------
# SOCI integral transformation (Psi4)
# Inactive=1 (freeze O 1s); active = all remaining 17 MOs.
# ncore=1, ncas=17, nelecas=8.  ecore / h1eff built per PySCF h1e_for_cas.
# ---------------------------------------------------------------------------
ncas, nelecas, ncore_pyscf = 17, 8, 1
t0 = time.time()

mints  = psi4.core.MintsHelper(wfn.basisset())
Ca     = np.asarray(wfn.Ca_subset("AO", "ALL"))         # (18, 18), energy-ordered
core   = Ca[:, 0:ncore_pyscf]                           # frozen O 1s
active = Ca[:, ncore_pyscf:ncore_pyscf + ncas]          # 17 active MOs
hcore  = np.asarray(wfn.H())                            # T + V (AO basis)
I_ao   = np.asarray(mints.ao_eri())                     # (18,18,18,18) chemists (pq|rs)

# Effective active-space 1e Hamiltonian + core energy (PySCF h1e_for_cas).
core_dm = 2.0 * core @ core.T
J = np.einsum("pqrs,rs->pq", I_ao, core_dm)
K = np.einsum("prqs,rs->pq", I_ao, core_dm)
corevhf = J - 0.5 * K
ecore = (nre
         + np.einsum("ij,ji", core_dm, hcore)
         + 0.5 * np.einsum("ij,ji", core_dm, corevhf))
h1eff = active.T @ (hcore + corevhf) @ active           # (17, 17)

Cact    = psi4.core.Matrix.from_array(active)
eri_act = np.asarray(
    mints.mo_eri(Cact, Cact, Cact, Cact)).reshape(ncas, ncas, ncas, ncas)
print(f"  Integrals ready ({time.time()-t0:.1f}s)  "
      f"h1eff {h1eff.shape}  eri_act {eri_act.shape}")

# ===========================================================================
# From the "hexpr SOCI expression" section onward, this file is byte-identical
# to run_pyscf.py (only the integral block above differs; run_pyscf.py also
# writes integrals.h5 at the end, which this script omits).
# ===========================================================================

# ---------------------------------------------------------------------------
# hexpr SOCI expression
# hexpr partition: ncore=2 (MO2-3), nval=4 (MO4-7), next=11 (MO8-18)
# ---------------------------------------------------------------------------
hexpr.init()
drt  = hexpr.Drt('soci', nve=nelecas, spin_x2=0, ncore=2, nval=4, next_=11)
expr = hexpr.Expression(drt)
t0 = time.time()
ij_arr, pqrs_arr, coef_arr = expr.terms_full()
ncsfs  = drt.ncsfs
norb   = drt.norb    # 17
n_oel  = norb * (norb + 1) // 2   # 153
print(f"  hexpr: ncsfs={ncsfs}  nterms={len(ij_arr)}  ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Scalar tri inverse (see docs/PQRS_INDEX_SPEC.md)
# ---------------------------------------------------------------------------
def tri_inv(t):
    """Decode tri(p, q) = p*(p-1)/2 + q to (p, q) (1-based)."""
    p = int(np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t)) / 2.0))
    if p * (p - 1) // 2 >= t:
        p -= 1
    q = t - p * (p - 1) // 2
    return p, q


def _csf_label(row):
    """Format step codes (0..3) as both e/u/d/f symbols and raw integers."""
    codes = [int(c) for c in row]
    sym  = "".join("eudf"[c] for c in codes)
    nums = " ".join(str(c) for c in codes)
    return sym, nums

# ---------------------------------------------------------------------------
# H matrix construction: one explicit pass over the (ij, pqrs, coef) terms.
# This is the same per-term loop the H2 example uses -- it shows plainly how to
# consume the expression: decode the packed CSF-pair index ij -> (I, J), decode
# the orbital index from pqrs, look up the MO integral, accumulate into H.
# A vectorized version of this same construction is given (commented) below.
# ---------------------------------------------------------------------------
# --- Scalar per-term loop (kept commented for clarity; the vectorized block
#     below is the active path -- see note there). -------------------------
# t0 = time.time()
# H = np.zeros((ncsfs, ncsfs))
# for k in range(len(ij_arr)):
#     # Decode the packed CSF-pair index ij -> (I, J), I >= J.
#     ij = int(ij_arr[k])
#     I  = int((np.sqrt(8 * ij + 1) - 1) // 2)
#     J  = ij - I * (I + 1) // 2
#     c  = coef_arr[k]
#
#     if pqrs_arr[k] < n_oel:
#         # 1-electron term: pqrs = tri(p, q) - 1
#         p, q = tri_inv(int(pqrs_arr[k]) + 1)
#         v = h1eff[p - 1, q - 1]
#     else:
#         # 2-electron term: pqrs = tri(tri(p,q), tri(r,s)) - 1 + n_oel
#         outer = int(pqrs_arr[k]) - n_oel + 1
#         a, b  = tri_inv(outer)
#         p, q  = tri_inv(a)
#         r, s  = tri_inv(b)
#         v = eri_act[p - 1, q - 1, r - 1, s - 1]
#
#     H[I, J] += c * v
#     if I != J:
#         H[J, I] += c * v
#
# print(f"  H built ({time.time()-t0:.1f}s)")

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

# ---------------------------------------------------------------------------
# Diagonalise
# ---------------------------------------------------------------------------
t0 = time.time()
evals, evecs = np.linalg.eigh(H)
E_active = evals
soci = ecore + E_active[0]
civec = evecs[:, 0]   # ground state (eigh returns ascending eigenvalues)
print(f"  eigh ({time.time()-t0:.1f}s)")

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print()
print("=" * 60)
print("H2O / 6-31G(d)  SOCI  full CI   (Psi4 + hexpr, SPHERICAL)")
print("=" * 60)
print(f"Nuclear repulsion: {nre:.10f} Hartree")
print(f"RHF energy:        {rhf:.10f} Hartree")
print(f"SOCI energy:       {soci:.10f} Hartree")
print(f"  ecore={ecore:.10f}  E_active[0]={E_active[0]:.10f}")
print("=" * 60)

# ---------------------------------------------------------------------------
# Ground-state CI vector (dominant configurations only).
# Ordering anchor: H was filled by decoding the packed pair index ij into DRT
# CSF indices (I, J), so civec[i] corresponds 1:1 to drt.csfs()[i] (DRT order).
# ---------------------------------------------------------------------------
step_rows = drt.csfs()   # uint8 (ncsfs, norb), DRT order -- matches civec
assert step_rows.shape[0] == len(civec), \
    f"drt.csfs() rows {step_rows.shape[0]} != civec length {len(civec)}"

WEIGHT_CUTOFF = 0.001   # show CSFs whose weight (c**2) is at least this
print()
print(f"Dominant CSFs (weight = c**2 >= {WEIGHT_CUTOFF}):")
w = civec**2
order = [i for i in np.argsort(-np.abs(civec)) if w[i] >= WEIGHT_CUTOFF]
for i in order:
    sym, nums = _csf_label(step_rows[i])
    print(f"  {sym}  ({nums})   c = {civec[i]:+.6f}   w = {w[i]:.6f}")
print(f"  ({len(order)} of {len(civec)} CSFs shown; weight cutoff {WEIGHT_CUTOFF})")
print("=" * 60)
