@testset "drt" begin
    HExpr.init()
    d = HExpr.Drt("soci"; nve=4, spin_x2=0, ncore=0, nval=4, next=0)
    @test HExpr.norb(d)      == 4
    @test HExpr.ncsfs(d)     == 20
    @test HExpr.nnodes(d)    == 14
    @test HExpr.kind(d)      == HExpr.HEXPR_DRT_SOCI
    @test HExpr.nve_state(d) == 4
    @test HExpr.spin_x2(d)   == 0
    finalize(d)
    HExpr.shutdown()
end
