using Test
using HExpr

@testset "HExpr" begin
    include("test_smoke.jl")
    include("test_drt.jl")
    include("test_expr.jl")
    include("test_csfset.jl")
    include("test_pair.jl")
    include("test_completeness.jl")
end
