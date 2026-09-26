# Example: H2 / STO-3G CAS(2,2) full CI via direct CSF set input.
#
# Demonstrates the high-level CSF-set workflow:
#   - Build a sub-DRT directly from an explicit list of CSF step-code
#     sequences (no need to declare nve / spin_x2 / orbital partition).
#   - Use HExpr.csf_subspace_eval to compute the CI Hamiltonian over the
#     square CSF subspace in one call.
#   - Diagonalize to obtain the three CASCI total energies.
#
# The three singlet CSFs for H2 CAS(2,2) are:
#   CSF 0:  F E  -- sigma^2          (closed-shell, both electrons in MO 0)
#   CSF 1:  E F  -- sigma*^2         (closed-shell, both electrons in MO 1)
#   CSF 2:  U D  -- sigma alpha,
#                   sigma* beta      (open-shell singlet coupling)
#
# This is a self-contained quantum-chemistry mini-program. With the default
# hard-coded integrals it runs H2 / STO-3G at R = 0.7414 Angstrom; uncomment
# the HDF5 block to read integrals computed elsewhere (e.g. by PySCF) for
# any geometry or basis. The file molcas_h2.ci.txt in this directory lists
# the OpenMolcas reference numbers for the default case if you want to
# spot-check the output.

using HExpr
using LinearAlgebra
using Printf

# ---------------------------------------------------------------------------
# Integrals -- hard-coded values from PySCF
#
# These integrals were computed by PySCF for H2 / STO-3G CAS(2,2) at the
# equilibrium geometry:
#
#     H   0.00000000   0.00000000   0.00000000
#     H   0.00000000   0.00000000   0.74140000   (Angstrom)
#
# Hard-coding keeps this example dependency-free (no PySCF or HDF5 required
# just to run it). To use a different geometry or basis, generate the
# integrals with PySCF (or another QC code), dump them to an HDF5 file,
# and use the HDF5 block below instead.
# ---------------------------------------------------------------------------
const ENUC  = 0.713753993687618
const ECORE = 0.713753993687618    # equals enuc (no doubly-occupied core MOs)

h1eff = [
    -1.252463573564898   0.0;
     0.0                -0.475948715220964
]

eri_act = zeros(Float64, 2, 2, 2, 2)
eri_act[1, 1, 1, 1] = 0.674488766356838
eri_act[2, 2, 2, 2] = 0.697393767423027
eri_act[1, 1, 2, 2] = eri_act[2, 2, 1, 1] = 0.663468096423568
eri_act[1, 2, 1, 2] = eri_act[1, 2, 2, 1] =
eri_act[2, 1, 1, 2] = eri_act[2, 1, 2, 1] = 0.181288808211496

# ---------------------------------------------------------------------------
# Alternative: load integrals from an HDF5 file
#
# Uncomment this block (and comment out the hard-coded values above) to
# read precomputed integrals from an HDF5 file. The file is expected to
# contain three datasets:
#   enuc    : scalar Float64
#   h1eff   : 2 x 2  Float64 array
#   eri_act : 2 x 2 x 2 x 2 Float64 array
# Generate such a file with PySCF: see the alternative block in run.py
# for the PySCF code, then dump h1eff, eri_act, and enuc with h5py.
#
# using HDF5
# ENUC, ECORE, h1eff, eri_act = h5open("integrals.h5", "r") do f
#     enuc = read(f["enuc"])
#     h1  = read(f["h1eff"])
#     eri = read(f["eri_act"])
#     (enuc, enuc, h1, eri)   # ECORE == ENUC for CAS(2,2) with no core
# end
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# CSF set: three singlet CSFs for H2 CAS(2,2)
# Step codes: 0=E (empty), 1=U (spin-up), 2=D (spin-down), 3=F (doubly occupied)
# Rows are CSFs.
# ---------------------------------------------------------------------------
csfs = UInt8[
    3 0;    # F E:  sigma^2
    0 3;    # E F:  sigma*^2
    1 2     # U D:  sigma alpha, sigma* beta (open-shell singlet)
]

norb  = 2
ncsfs = 3

# ---------------------------------------------------------------------------
# Evaluate the CI Hamiltonian over the 3-CSF subspace
# ---------------------------------------------------------------------------
HExpr.init()

ij, pqrs, coef, bra_pos, ket_pos = HExpr.csf_subspace_eval(csfs; which=:full)

println("HExpr.csf_subspace_eval: $(length(ij)) terms over $ncsfs CSFs")

# ---------------------------------------------------------------------------
# Build the CI Hamiltonian matrix H[I, J]
# Decode each term's (p, q, r, s) from the packed pqrs index.
# See docs/PQRS_INDEX_SPEC.md for the full encoding.
# ---------------------------------------------------------------------------
const n_oel = norb * (norb + 1) ÷ 2   # 3, number of 1e index pairs

"Decode tri(p, q) = p*(p-1)/2 + q to (p, q) (1-based)."
function tri_inv(t::Integer)
    p = Int(ceil((-1.0 + sqrt(1.0 + 8.0 * t)) / 2.0))
    if p * (p - 1) ÷ 2 >= t
        p -= 1
    end
    q = t - p * (p - 1) ÷ 2
    return p, q
end

H = zeros(Float64, ncsfs, ncsfs)
for k in 1:length(ij)
    # bra_pos / ket_pos are 0-based positions from hexpr; add 1 for Julia.
    I = Int(bra_pos[k]) + 1
    J = Int(ket_pos[k]) + 1
    c = coef[k]
    if pqrs[k] < n_oel
        # 1-electron term
        p, q = tri_inv(Int(pqrs[k]) + 1)
        v = h1eff[p, q]
    else
        # 2-electron term
        outer = Int(pqrs[k]) - n_oel + 1
        a, b  = tri_inv(outer)
        p, q  = tri_inv(a)
        r, s  = tri_inv(b)
        v = eri_act[p, q, r, s]
    end
    H[I, J] += c * v
    if I != J
        H[J, I] += c * v
    end
end

# ---------------------------------------------------------------------------
# Diagonalize and report the three total energies (ecore + active energy)
# ---------------------------------------------------------------------------
E_active = eigvals(Symmetric(H))
E_total  = ECORE .+ E_active

println()
println("=" ^ 60)
println("H2 / STO-3G  CAS(2,2)  full CI  (3 singlet CSFs)")
println("=" ^ 60)
for k in 1:ncsfs
    @printf("Root %d: %+.10f Hartree\n", k, E_total[k])
end
println("=" ^ 60)
