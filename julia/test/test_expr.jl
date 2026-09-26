@testset "expr" begin
    HExpr.init()
    d = HExpr.Drt("soci"; nve=4, spin_x2=0, ncore=0, nval=4, next=0)
    e = HExpr.Expression(d)

    @test HExpr.nterms_full(e) == 638
    @test HExpr.nterms_one(e)  == 124

    ij, pqrs, coef = HExpr.read_full(e)
    @test length(ij)   == 638
    @test length(pqrs) == 638
    @test length(coef) == 638
    @test ij[1]   == 0
    @test pqrs[1] == 5
    @test isapprox(coef[1], 2.0)

    ij1, pqrs1, coef1 = HExpr.read_one(e)
    @test length(ij1)   == 124
    @test length(pqrs1) == 124
    @test length(coef1) == 124

    finalize(e)
    finalize(d)
    HExpr.shutdown()
end
