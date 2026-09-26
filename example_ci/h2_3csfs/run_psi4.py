"""
Example: H2 / STO-3G CAS(2,2) full CI via direct CSF set input,
with integrals computed by Psi4.

Psi4 sibling of run_pyscf.py: identical calculation, but ENUC / ECORE /
h1eff / eri_act come from Psi4 instead of PySCF. Everything from the CSF
set onward (the hexpr csf_subspace_eval evaluation, H build, diagonalize,
CI-vector report) is reused verbatim from run_pyscf.py -- only the
integral-generation block differs.

CAS(2,2) has no frozen core (ncore = 0), so ECORE == ENUC and h1eff is just
the bare active-space one-electron Hamiltonian (no core J/K term).

Requires Psi4. Note: psi4 1.9.1 needs libint2.so.2; install with libint
pinned below 2.10 (e.g. `conda install "psi4=1.9.1" "libint<2.10"` ->
libint 2.9.0). See example_ci/README.
"""

import numpy as np
import psi4
import hexpr

# ---------------------------------------------------------------------------
# Molecule and RHF (Psi4)
# Geometry copied verbatim from run_pyscf.py. no_com / no_reorient keep Psi4
# in the same frame PySCF uses (PySCF does not recenter/reorient at
# symmetry=False); immaterial to the energy, kept for a clean comparison.
# ---------------------------------------------------------------------------
psi4.core.be_quiet()
mol = psi4.geometry("""
H   0.00000000   0.00000000   0.00000000
H   0.00000000   0.00000000   0.74140000
units angstrom
symmetry c1
no_com
no_reorient
""")
psi4.set_options({
    "basis": "sto-3g",
    "scf_type": "pk",          # exact ERIs (no density fitting)
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
})
rhf, wfn = psi4.energy("scf", return_wfn=True)
assert wfn.basisset().nbf() == 2, f"expected nbf=2, got {wfn.basisset().nbf()}"

# ---------------------------------------------------------------------------
# CAS(2,2) integrals (no frozen core: ncore = 0, all MOs active)
# ---------------------------------------------------------------------------
mints   = psi4.core.MintsHelper(wfn.basisset())
Ca      = np.asarray(wfn.Ca_subset("AO", "ALL"))      # (2, 2), energy-ordered
active  = Ca[:, 0:2]                                   # all MOs active (ncore = 0)
hcore   = np.asarray(wfn.H())                          # T + V (AO basis)

ENUC  = mol.nuclear_repulsion_energy()
ECORE = ENUC                                           # ncore = 0  ->  no core, corevhf = 0
h1eff = active.T @ hcore @ active                      # (2, 2)

Cact    = psi4.core.Matrix.from_array(active)
eri_act = np.asarray(mints.mo_eri(Cact, Cact, Cact, Cact)).reshape(2, 2, 2, 2)
print(f"  RHF = {rhf:.10f}  ENUC = {ENUC:.10f}")

# ===========================================================================
# Everything below is reused verbatim from run_pyscf.py
# ===========================================================================

# ---------------------------------------------------------------------------
# CSF set: three singlet CSFs for H2 CAS(2,2)
# Step codes: 0=E (empty), 1=U (spin-up), 2=D (spin-down), 3=F (doubly occupied)
# ---------------------------------------------------------------------------
csfs = np.array([
    [3, 0],   # F E:  sigma^2
    [0, 3],   # E F:  sigma*^2
    [1, 2],   # U D:  sigma alpha, sigma* beta (open-shell singlet)
], dtype=np.uint8)

norb  = 2
ncsfs = 3

# ---------------------------------------------------------------------------
# Evaluate the CI Hamiltonian over the 3-CSF subspace
# ---------------------------------------------------------------------------
hexpr.init()

ij, pqrs, coef, bra_pos, ket_pos = hexpr.csf_subspace_eval(
    csfs, which='full')

print(f"hexpr.csf_subspace_eval: {len(ij)} terms over {ncsfs} CSFs")

# ---------------------------------------------------------------------------
# Build the CI Hamiltonian matrix H[I, J]
# Decode each term's (p, q, r, s) from the packed pqrs index.
# See docs/PQRS_INDEX_SPEC.md for the full encoding.
# ---------------------------------------------------------------------------
n_oel = norb * (norb + 1) // 2   # 3, number of 1e index pairs

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

H = np.zeros((ncsfs, ncsfs))
for k in range(len(ij)):
    I, J = bra_pos[k], ket_pos[k]
    c = coef[k]
    if pqrs[k] < n_oel:
        # 1-electron term
        p, q = tri_inv(int(pqrs[k]) + 1)
        v = h1eff[p - 1, q - 1]
    else:
        # 2-electron term
        outer = int(pqrs[k]) - n_oel + 1
        a, b  = tri_inv(outer)
        p, q  = tri_inv(a)
        r, s  = tri_inv(b)
        v = eri_act[p - 1, q - 1, r - 1, s - 1]
    H[I, J] += c * v
    if I != J:
        H[J, I] += c * v

# ---------------------------------------------------------------------------
# Diagonalize and report the three total energies (ecore + active energy)
# ---------------------------------------------------------------------------
evals, evecs = np.linalg.eigh(H)
E_active = evals
E_total  = ECORE + E_active   # eigh returns ascending eigenvalues

print()
print("=" * 60)
print("H2 / STO-3G  CAS(2,2)  full CI  (3 singlet CSFs)")
print("=" * 60)
for k, e in enumerate(E_total, start=1):
    print(f"Root {k}: {e:+.10f} Hartree")
print("=" * 60)

# ---------------------------------------------------------------------------
# Per-root CI vectors.
# Ordering anchor: H was filled by bra_pos/ket_pos (positions into the input
# `csfs` array), so evecs[:, k][i] corresponds to csfs[i] for every root k
# -- NOT a DRT ordering (there is no Drt object in this script). The
# amplitude sign is an arbitrary global phase per root.
# ---------------------------------------------------------------------------
assert evecs.shape[0] == csfs.shape[0], \
    f"evecs rows {evecs.shape[0]} != csfs rows {csfs.shape[0]}"

print()
print("Per-root CI vectors (all CSFs):")
for k in range(ncsfs):
    civec_k = evecs[:, k]
    print(f"\nRoot {k + 1}:  E = {E_total[k]:+.10f} Hartree")
    for i in np.argsort(-np.abs(civec_k)):
        sym, nums = _csf_label(csfs[i])
        print(f"  {sym}  ({nums})   c = {civec_k[i]:+.6f}   w = {civec_k[i]**2:.6f}")
print("=" * 60)
