#!/usr/bin/env julia
# =============================================================================
# H2O / 6-31G(d)  SOCI full CI  via  Fermi.jl integrals + hexpr.
#
# Julia analogue of run_pyscf.py, but the RHF reference and MO integrals come
# from Fermi.jl instead of PySCF. Fermi/GaussianBasis only supports SPHERICAL
# basis functions, so this uses spherical d (5 components) -> 18 AOs/MOs, vs.
# the cartesian (6 components) -> 19 in the PySCF/OpenMolcas reference. Freezing
# the O 1s core leaves 17 active orbitals (one fewer than the cartesian 18), so
# the hexpr SOCI partition uses next=11 (4-orbital reference RAS2):
#     ncore=2, nval=4, next=11   (norb = 17),   nve = 8
#
# Pipeline:
#   - Fermi.jl computes RHF and the MO 1e (T+V) and 2e (chemist (pq|rs)) ints.
#   - the O 1s core (MO 1) is folded into an effective 1e Hamiltonian + ecore.
#   - hexpr builds the SOCI DRT and emits the full CI expression (ij,pqrs,coef).
#   - the CI Hamiltonian is assembled term-by-term, diagonalized, lowest state
#     reported.
#
# Run with both the Fermi project and the hexpr Julia binding on the load path:
#   export JULIA_DEPOT_PATH="$PWD/local:$HOME/.julia"
#   export HEXPR_PREFIX=<your hexpr install prefix>
#   export HEXPR_LIBRARY=$HEXPR_PREFIX/lib/libhexpr.so
#   export LD_LIBRARY_PATH=$HEXPR_PREFIX/lib:$LD_LIBRARY_PATH
#   export JULIA_LOAD_PATH=$HEXPR_PREFIX/share/hexpr/julia:
#   julia --project=<your Fermi.jl checkout> run_fermi.jl
# =============================================================================

using Fermi
using Fermi.Integrals: IntegralHelper, compute!
import Molecules
using LinearAlgebra
using HExpr

# ---------------------------------------------------------------------------
# Molecule and RHF (spherical 6-31G*)
# ---------------------------------------------------------------------------
Fermi.Options.set("printstyle", "none")
Fermi.Options.set("basis", "6-31g*")
@molecule {
  O   0.0   0.0000000   0.0000000
  H   0.0   0.7572157   0.5865358
  H   0.0  -0.7572157   0.5865358
}

aoints = IntegralHelper{Float64}()
rhf    = Fermi.HartreeFock.RHF(aoints)

mol  = aoints.molecule
nre  = Molecules.nuclear_repulsion(mol.atoms)
rhf_E = rhf.energy
println("  NRE = ", round(nre;   digits=10), "   RHF = ", round(rhf_E; digits=10))

# ---------------------------------------------------------------------------
# MO integrals in the RHF basis  (hp = T+V, eri[p,q,r,s] = (pq|rs) chemist)
# ---------------------------------------------------------------------------
moints = IntegralHelper(orbitals=rhf.orbitals)
hp  = compute!(moints, aoints, "T")
hp += compute!(moints, aoints, "V")
eri = compute!(moints, "ERI")
nbas = size(hp, 1)                       # 18 (spherical)

# ---------------------------------------------------------------------------
# Freeze the O 1s core (MO 1) -> effective 1e Hamiltonian + core energy.
# chemist notation: eri[p,q,r,s] = (pq|rs)
#   Ecore   = Vnuc + Σ_i 2 h_ii + Σ_ij [ 2(ii|jj) - (ij|ji) ]     (i,j ∈ core)
#   h1eff_pq = h_pq + Σ_i [ 2(pq|ii) - (pi|iq) ]                   (p,q ∈ active)
# ---------------------------------------------------------------------------
nfrozen = 1
core = 1:nfrozen                          # MO 1
act  = (nfrozen + 1):nbas                 # MO 2..18  -> 17 active orbitals
nact = length(act)

Ecore = nre + 2 * sum(hp[i, i] for i in core) +
        sum(2 * eri[i, i, j, j] - eri[i, j, j, i] for i in core, j in core)

h1eff = zeros(nact, nact)
for (pp, p) in enumerate(act), (qq, q) in enumerate(act)
    v = hp[p, q]
    for i in core
        v += 2 * eri[p, q, i, i] - eri[p, i, i, q]
    end
    h1eff[pp, qq] = v
end

eri_act = eri[act, act, act, act]         # nact^4, active-local indices

# ---------------------------------------------------------------------------
# hexpr SOCI expression
#   hexpr partition over the 17 active orbitals: ncore=2, nval=4, next=11
# ---------------------------------------------------------------------------
HExpr.init()
nve = 8
drt = HExpr.Drt("soci"; nve=nve, spin_x2=0, ncore=2, nval=4, next=nact - 6)
@assert HExpr.norb(drt) == nact "hexpr norb $(HExpr.norb(drt)) != nact $nact"

expr = HExpr.Expression(drt)
ij_arr, pqrs_arr, coef_arr = HExpr.read_full(expr)
ncsfs = Int(HExpr.ncsfs(drt))
norb  = Int(HExpr.norb(drt))              # 17
n_oel = norb * (norb + 1) ÷ 2             # 153
println("  hexpr: ncsfs=", ncsfs, "  nterms=", length(ij_arr),
        "  (norb=", norb, ", n_oel=", n_oel, ")")

# ---------------------------------------------------------------------------
# tri(p,q) = p(p-1)/2 + q  inverse  (1-based, p >= q)
# ---------------------------------------------------------------------------
function tri_inv(t::Int)
    p = Int(ceil((-1.0 + sqrt(1.0 + 8.0 * t)) / 2.0))
    if p * (p - 1) ÷ 2 >= t
        p -= 1
    end
    q = t - p * (p - 1) ÷ 2
    return p, q
end

# ---------------------------------------------------------------------------
# Assemble H: one pass over the (ij, pqrs, coef) terms (hexpr indices 0-based).
# ---------------------------------------------------------------------------
H = zeros(ncsfs, ncsfs)
for k in eachindex(ij_arr)
    ij = Int(ij_arr[k])
    I  = Int(floor((sqrt(8.0 * ij + 1.0) - 1.0) / 2.0))   # 0-based bra
    J  = ij - I * (I + 1) ÷ 2                              # 0-based ket
    c  = coef_arr[k]
    pq = Int(pqrs_arr[k])

    if pq < n_oel
        # 1-electron term: pqrs = tri(p,q) - 1
        p, q = tri_inv(pq + 1)
        v = h1eff[p, q]
    else
        # 2-electron term: pqrs = tri(tri(p,q), tri(r,s)) - 1 + n_oel
        outer = pq - n_oel + 1
        a, b = tri_inv(outer)
        p, q = tri_inv(a)
        r, s = tri_inv(b)
        v = eri_act[p, q, r, s]
    end

    H[I + 1, J + 1] += c * v               # +1 for Julia 1-based CSF index
    if I != J
        H[J + 1, I + 1] += c * v
    end
end

# ---------------------------------------------------------------------------
# Validation: the reference CSF (4 active orbitals doubly occupied: ffff0...)
# is the RHF determinant, so  ecore + <ref|H|ref>  must equal the RHF energy.
# ---------------------------------------------------------------------------
ref_steps = vcat(fill(UInt8(3), 4), fill(UInt8(0), norb - 4))
ref = HExpr.csf_index(drt, ref_steps) + 1     # 0-based -> 1-based
E_ref = Ecore + H[ref, ref]
println("  check: ecore + <ref|H|ref> = ", round(E_ref; digits=10),
        "   (RHF = ", round(rhf_E; digits=10), ",  Δ = ", E_ref - rhf_E, ")")
@assert isapprox(E_ref, rhf_E; atol=1e-8) "reference CSF energy != RHF"

# ---------------------------------------------------------------------------
# Diagonalize
# ---------------------------------------------------------------------------
F = eigen(Symmetric(H))
E_active = F.values
soci  = Ecore + E_active[1]
civec = F.vectors[:, 1]

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
println()
println("="^60)
println("H2O / 6-31G(d)  SOCI  full CI   (Fermi.jl + hexpr, spherical)")
println("="^60)
println("Nuclear repulsion: ", round(nre;   digits=10), " Hartree")
println("RHF energy:        ", round(rhf_E;  digits=10), " Hartree")
println("SOCI energy:       ", round(soci;   digits=10), " Hartree")
println("  ecore=", round(Ecore; digits=10), "  E_active[1]=", round(E_active[1]; digits=10))
println("="^60)

# Dominant CSFs (weight = c^2). csfs(drt) rows are in DRT order, matching civec.
step_rows = HExpr.csfs(drt)               # (ncsfs, norb) UInt8
@assert size(step_rows, 1) == length(civec)
sym = ('e', 'u', 'd', 'f')
wcut = 0.001
w = civec .^ 2
order = [i for i in sortperm(abs.(civec); rev=true) if w[i] >= wcut]
println()
println("Dominant CSFs (weight = c^2 >= ", wcut, "):")
for i in order
    codes = Int.(step_rows[i, :])
    label = join(sym[c + 1] for c in codes)
    nums  = join(string.(codes), " ")
    println("  ", label, "  (", nums, ")   c = ",
            round(civec[i]; digits=6), "   w = ", round(w[i]; digits=6))
end
println("  (", length(order), " of ", length(civec), " CSFs shown; cutoff ", wcut, ")")
println("="^60)

HExpr.shutdown()
