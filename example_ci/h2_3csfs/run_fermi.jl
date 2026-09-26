# Example: H2 / STO-3G CAS(2,2) full CI via direct CSF set input,
# with integrals from Fermi.jl.
#
# Julia analogue of run_pyscf.py / run_psi4.py for the h2_3csfs example, and
# the integral-driven counterpart to run.jl (which uses hard-coded integrals).
# Here the RHF reference and the CAS(2,2) MO integrals are computed by Fermi.jl
# -- which builds its one- and two-electron integrals through the
# GaussianBasis.jl library -- and the CI Hamiltonian over the explicit 3-CSF
# set is built with hexpr's csf_subspace_eval (no DRT; the CSF set is given
# directly).
#
# The three singlet CSFs for H2 CAS(2,2) are:
#   CSF 0:  F E  -- sigma^2     (closed shell, both electrons in MO 0)
#   CSF 1:  E F  -- sigma*^2    (closed shell, both electrons in MO 1)
#   CSF 2:  U D  -- open-shell singlet (sigma alpha, sigma* beta)
#
# STO-3G has only s functions on H (nbf = 2), so "spherical" vs "cartesian" is
# immaterial here, and CAS(2,2) has no frozen core (ECORE == ENUC; h1eff is the
# bare active-space 1e Hamiltonian, T+V in the MO basis).
#
# This engine needs an external Fermi.jl checkout and its depot (it is not
# covered by example_ci/Project.toml); see example_ci/README.md for the
# Fermi.jl env and launch.

using Fermi
using Fermi.Integrals: IntegralHelper, compute!
import Molecules
using LinearAlgebra
using HExpr

# ---------------------------------------------------------------------------
# Molecule and RHF (STO-3G)
# ---------------------------------------------------------------------------
Fermi.Options.set("printstyle", "none")
Fermi.Options.set("basis", "sto-3g")
# DIIS / ODA convergence accelerators divide by ~0 on this 2-AO / 1-occupied
# system and poison the Fock matrix with NaNs; plain SCF converges fine here.
Fermi.Options.set("diis", false)
Fermi.Options.set("oda", false)
@molecule {
  H   0.00000000   0.00000000   0.00000000
  H   0.00000000   0.00000000   0.74140000
}

aoints = IntegralHelper{Float64}()
rhf    = Fermi.HartreeFock.RHF(aoints)

mol = aoints.molecule
nre = Molecules.nuclear_repulsion(mol.atoms)
println("  RHF = ", round(rhf.energy; digits=10), "   ENUC = ", round(nre; digits=10))

# ---------------------------------------------------------------------------
# CAS(2,2) MO integrals in the RHF basis (no frozen core -> all 2 MOs active).
#   h1eff = T + V (MO basis),   eri_act[p,q,r,s] = (pq|rs) chemist
# ---------------------------------------------------------------------------
moints = IntegralHelper(orbitals=rhf.orbitals)
h1eff  = compute!(moints, aoints, "T")
h1eff += compute!(moints, aoints, "V")
eri_act = compute!(moints, "ERI")

norb  = size(h1eff, 1)            # 2
ECORE = nre                       # CAS(2,2): no core, ECORE == ENUC

# ---------------------------------------------------------------------------
# CSF set: three singlet CSFs for H2 CAS(2,2)
# Step codes: 0=E (empty), 1=U (spin-up), 2=D (spin-down), 3=F (doubly occupied)
# ---------------------------------------------------------------------------
csfs = UInt8[3 0;     # F E:  sigma^2
             0 3;     # E F:  sigma*^2
             1 2]     # U D:  open-shell singlet
ncsfs = size(csfs, 1)

# ---------------------------------------------------------------------------
# hexpr: evaluate the CI Hamiltonian over the 3-CSF subspace.
# Returns (ij, pqrs, coef, bra_pos, ket_pos); bra_pos/ket_pos are 0-based
# positions into `csfs`. (ij is unused here -- we index H by bra_pos/ket_pos.)
# ---------------------------------------------------------------------------
HExpr.init()
ij, pqrs, coef, bra_pos, ket_pos = HExpr.csf_subspace_eval(csfs; which=:full)
println("  hexpr.csf_subspace_eval: ", length(ij), " terms over ", ncsfs, " CSFs")

n_oel = norb * (norb + 1) ÷ 2     # 3

# tri(p,q) = p(p-1)/2 + q  inverse  (1-based, p >= q)
function tri_inv(t::Int)
    p = Int(ceil((-1.0 + sqrt(1.0 + 8.0 * t)) / 2.0))
    if p * (p - 1) ÷ 2 >= t
        p -= 1
    end
    q = t - p * (p - 1) ÷ 2
    return p, q
end

# ---------------------------------------------------------------------------
# Build H[I, J] (hexpr orbital indices 0-based; bra_pos/ket_pos 0-based).
# ---------------------------------------------------------------------------
H = zeros(ncsfs, ncsfs)
for k in eachindex(ij)
    I = Int(bra_pos[k])
    J = Int(ket_pos[k])
    c = coef[k]
    pq = Int(pqrs[k])
    if pq < n_oel
        p, q = tri_inv(pq + 1)
        v = h1eff[p, q]
    else
        outer = pq - n_oel + 1
        a, b = tri_inv(outer)
        p, q = tri_inv(a)
        r, s = tri_inv(b)
        v = eri_act[p, q, r, s]
    end
    H[I + 1, J + 1] += c * v
    if I != J
        H[J + 1, I + 1] += c * v
    end
end

# ---------------------------------------------------------------------------
# Diagonalize and report the three total energies (ECORE + active energy).
# ---------------------------------------------------------------------------
F = eigen(Symmetric(H))
E_total = ECORE .+ F.values

println()
println("="^60)
println("H2 / STO-3G  CAS(2,2)  full CI  (3 singlet CSFs, Fermi.jl + hexpr)")
println("="^60)
for k in 1:ncsfs
    println("Root ", k, ": ", (F.values[k] >= 0 ? "+" : ""),
            round(E_total[k]; digits=10), " Hartree")
end
println("="^60)

# ---------------------------------------------------------------------------
# Per-root CI vectors. evecs[:,k][i] corresponds to csfs[i,:] (input order).
# The amplitude sign is an arbitrary global phase per root.
# ---------------------------------------------------------------------------
sym = ('e', 'u', 'd', 'f')
println()
println("Per-root CI vectors (all CSFs):")
for k in 1:ncsfs
    civec = F.vectors[:, k]
    println()
    println("Root ", k, ":  E = ", round(E_total[k]; digits=10), " Hartree")
    for i in sortperm(abs.(civec); rev=true)
        codes = Int.(csfs[i, :])
        label = join(sym[c + 1] for c in codes)
        nums  = join(string.(codes), " ")
        println("  ", label, "  (", nums, ")   c = ", round(civec[i]; digits=6),
                "   w = ", round(civec[i]^2; digits=6))
    end
end
println("="^60)

HExpr.shutdown()
