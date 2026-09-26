"""
Example: H2 / STO-3G CAS(2,2) full CI via direct CSF set input,
with integrals computed by PySCF.

Same calculation as run.py, but builds the molecule and computes
ENUC / ECORE / h1eff / eri_act from PySCF rather than using
hard-coded values. Modify the molecule definition below to study
any geometry / basis / charge / spin.

Requires PySCF (`pip install pyscf`).
"""

import sys
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo
import hexpr

# ---------------------------------------------------------------------------
# Molecule
# Modify atom coordinates / basis / charge / spin to study different cases.
# ---------------------------------------------------------------------------
mol = gto.Mole()
mol.atom = """
H   0.00000000   0.00000000   0.00000000
H   0.00000000   0.00000000   0.74140000
"""
mol.basis = "sto-3g"
mol.charge = 0
mol.spin = 0
mol.unit = "Angstrom"
mol.symmetry = False
mol.build()

# ---------------------------------------------------------------------------
# RHF and CAS(2,2) integral transformation
# ---------------------------------------------------------------------------
mf = scf.RHF(mol)
mf.kernel()

ncas, nelecas = 2, 2
mc = mcscf.CASCI(mf, ncas, nelecas)
h1eff, ecore = mc.get_h1eff()
mo_act  = mf.mo_coeff[:, mc.ncore:mc.ncore + mc.ncas]
eri_act = ao2mo.kernel(mol, mo_act, compact=False).reshape(
    ncas, ncas, ncas, ncas)
ENUC  = mol.energy_nuc()
ECORE = ecore   # for CAS(2,2) with no core, ECORE == ENUC

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
