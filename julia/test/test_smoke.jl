@testset "smoke" begin
    HExpr.init()
    s = HExpr.version_string()
    @test !isempty(s)
    @test HExpr.version_major() == 1
    @test HExpr.version_minor() == 0
    @test HExpr.version_patch() == 0
    HExpr.shutdown()
end
