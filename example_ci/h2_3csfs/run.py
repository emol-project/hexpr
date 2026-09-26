"""
Example: H2 / STO-3G CAS(2,2) full CI via direct CSF set input.

Demonstrates the high-level CSF-set workflow:
  - Build a sub-DRT directly from an explicit list of CSF step-code
    sequences (no need to declare nve / spin_x2 / orbital partition).
  - Use hexpr.csf_subspace_eval to compute the CI Hamiltonian over the
    square CSF subspace in one call.
  - Diagonalize to obtain the three CASCI total energies.

The three singlet CSFs for H2 CAS(2,2) are:
  CSF 0:  F E  -- sigma^2          (closed-shell, both electrons in MO 0)
  CSF 1:  E F  -- sigma*^2         (closed-shell, both electrons in MO 1)
  CSF 2:  U D  -- sigma alpha,
                  sigma* beta      (open-shell singlet coupling)

This script uses hard-coded integrals for H2 / STO-3G at
R = 0.7414 Angstrom, so it has no PySCF dependency. For an
arbitrary geometry or basis, see the companion run_pyscf.py
which computes the integrals from scratch with PySCF.

The file molcas_h2.ci.txt in this directory lists the OpenMolcas
reference numbers for the default case if you want to spot-check
the output.
"""

import sys
import numpy as np
import hexpr

# ---------------------------------------------------------------------------
# Integrals -- hard-coded values from PySCF
#
# These integrals were computed by PySCF for H2 / STO-3G CAS(2,2) at the
# equilibrium geometry:
#
#     H   0.00000000   0.00000000   0.00000000
#     H   0.00000000   0.00000000   0.74140000   (Angstrom)
#
# For an arbitrary geometry or basis, use run_pyscf.py instead, which
# computes ENUC / ECORE / h1eff / eri_act from PySCF.
# ---------------------------------------------------------------------------
ENUC  = 0.713753993687618
ECORE = 0.713753993687618   # equals enuc (no doubly-occupied core MOs)

h1eff = np.array([
    [-1.252463573564898,  0.0               ],
    [ 0.0,               -0.475948715220964 ],
])

eri_act = np.zeros((2, 2, 2, 2))
eri_act[0, 0, 0, 0] = 0.674488766356838
eri_act[1, 1, 1, 1] = 0.697393767423027
eri_act[0, 0, 1, 1] = eri_act[1, 1, 0, 0] = 0.663468096423568
eri_act[0, 1, 0, 1] = eri_act[0, 1, 1, 0] = \
eri_act[1, 0, 0, 1] = eri_act[1, 0, 1, 0] = 0.181288808211496

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
E_active = np.linalg.eigh(H)[0]
E_total  = ECORE + E_active

print()
print("=" * 60)
print("H2 / STO-3G  CAS(2,2)  full CI  (3 singlet CSFs)")
print("=" * 60)
for k, e in enumerate(E_total, start=1):
    print(f"Root {k}: {e:+.10f} Hartree")
print("=" * 60)
