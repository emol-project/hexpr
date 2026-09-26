# Example: H2 / STO-3G CAS(2,2) -- generating the hexpr expression in Ruby.
#
# Demonstrates the high-level CSF-set workflow from Ruby:
#   - Build a sub-DRT directly from an explicit list of CSF step-code
#     sequences via HExpr::Drt.from_csfset.
#   - Construct the expression with HExpr::Expression.new(drt) and query
#     the term counts.
#
# Unlike run.py and run.jl, this script does NOT build the CI Hamiltonian
# or diagonalize -- it just shows how to drive hexpr from Ruby and reports
# how many expression terms come out. To do a full CI in Ruby, decode the
# (ij, pqrs, coef) arrays the same way run.py does. See the alternative
# block at the bottom for a sketch.
#
# The three singlet CSFs for H2 CAS(2,2) are:
#   CSF 0:  F E  -- sigma^2          (closed-shell, both electrons in MO 0)
#   CSF 1:  E F  -- sigma*^2         (closed-shell, both electrons in MO 1)
#   CSF 2:  U D  -- sigma alpha,
#                   sigma* beta      (open-shell singlet coupling)

require "hexpr"

HExpr.init

# ---------------------------------------------------------------------------
# CSF set: three singlet CSFs for H2 CAS(2,2)
# Step codes: 0=E (empty), 1=U (spin-up), 2=D (spin-down), 3=F (doubly occupied)
# ---------------------------------------------------------------------------
csfs = [
  [3, 0],   # F E:  sigma^2
  [0, 3],   # E F:  sigma*^2
  [1, 2],   # U D:  sigma alpha, sigma* beta (open-shell singlet)
]
norb = 2

# ---------------------------------------------------------------------------
# Build the sub-DRT and the expression
# ---------------------------------------------------------------------------
drt  = HExpr::Drt.from_csfset(csfs, norb)
expr = HExpr::Expression.new(drt)

puts "=" * 60
puts "H2 / STO-3G  CAS(2,2)  hexpr expression  (3 singlet CSFs)"
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
# nterms_full. For H2 / CAS(2,2) that's 15 terms, harmless to print -- but
# for a typical SOCI calculation it can easily be 10^5-10^6 terms, so don't
# enable this carelessly on a big system.
#
 ij, pqrs, coef = expr.read_full
 ij.each_with_index do |ij_k, k|
   puts "  k=#{k}  ij=#{ij_k}  pqrs=#{pqrs[k]}  coef=#{coef[k]}"
 end
# ---------------------------------------------------------------------------

HExpr.shutdown
