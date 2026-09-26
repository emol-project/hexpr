# Example: H2O / 6-31G(d) SOCI -- generating the hexpr expression in Ruby.
#
# Demonstrates how to drive hexpr from Ruby for a realistic SOCI workflow:
#   - Build the SOCI(8,17) DRT for H2O / 6-31G(d).
#   - Construct the expression and query the term counts.
#
# Unlike run_pyscf.py and run_fermi.jl, this script does NOT build the CI Hamiltonian
# or diagonalize -- it just shows how to drive hexpr from Ruby and reports
# how many expression terms come out. The CI evaluation requires MO
# integrals (one- and two-electron in the active space), which we don't
# attempt to compute in Ruby. To run a full CI, follow the same pattern
# as run_pyscf.py: get integrals from PySCF, decode the (ij, pqrs, coef)
# arrays, assemble H,
# and diagonalize.

require "hexpr"

HExpr.init

# ---------------------------------------------------------------------------
# Build the SOCI(8,17) DRT and the full expression
# hexpr partition: ncore=2 (MO2-3), nval=4 (MO4-7), next=11 (MO8-18)
# ---------------------------------------------------------------------------
drt  = HExpr::Drt.new(kind: "soci", nve: 8, spin_x2: 0,
                       ncore: 2, nval: 4, next: 11)
expr = HExpr::Expression.new(drt)

puts "=" * 60
puts "H2O / 6-31G(d)  SOCI  hexpr expression"
puts "=" * 60
puts "norb         = #{drt.norb}"
puts "ncsfs        = #{drt.ncsfs}"
puts "nnodes       = #{drt.nnodes}"
puts "nterms_full  = #{expr.nterms_full}"
puts "nterms_one   = #{expr.nterms_one}"
puts "=" * 60

# ---------------------------------------------------------------------------
# Dumping the full expression
#
# expr.read_full returns three parallel arrays (ij, pqrs, coef) of length
# nterms_full. For SOCI(8,17) with the 4e4o reference that's about 21 million
# terms -- DO NOT print them all unless you really want to (minutes of runtime,
# many MB of output). This block is left commented out by design.
#
# ij, pqrs, coef = expr.read_full
# ij.each_with_index do |ij_k, k|
#   puts "  k=#{k}  ij=#{ij_k}  pqrs=#{pqrs[k]}  coef=#{coef[k]}"
# end
# ---------------------------------------------------------------------------

HExpr.shutdown
