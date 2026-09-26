@testset "Drt.from_csfset" begin
    HExpr.init()

    full = HExpr.Drt("soci"; nve=4, spin_x2=0, ncore=0, nval=4, next=0)
    @test HExpr.norb(full)   == 4
    @test HExpr.ncsfs(full)  == 20
    @test HExpr.nnodes(full) == 14

    csfs_mat = HExpr.csfs(full)
    @test size(csfs_mat) == (20, 4)

    # Round-trip: sub-DRT from full CSF set matches full DRT shape
    sub = HExpr.Drt(csfs_mat, 4)
    @test HExpr.norb(sub)   == 4
    @test HExpr.ncsfs(sub)  == 20
    @test HExpr.nnodes(sub) == 14

    # 2-CSF subset has fewer nodes than full
    sub2 = HExpr.Drt(csfs_mat[1:2, :], 4)
    @test HExpr.norb(sub2)   == 4
    @test HExpr.ncsfs(sub2)  >= 2
    @test HExpr.nnodes(sub2) < HExpr.nnodes(full)

    # CSFSET: csfs()/ncsfs speak input order; walk_* exposes the superset.
    # Non-canonical input triggering a canonical superset: N=2, M=4.
    inp = UInt8[1 0 3 2;   # uefd
                0 1 2 3]   # eudf
    subs = HExpr.Drt(inp, 4)
    @test HExpr.ncsfs(subs) == 2
    @test HExpr.csfs(subs) == inp                 # input rows, input order
    @test HExpr.walk_ncsfs(subs) == 4             # canonical superset M
    @test HExpr.walk_ncsfs(subs) >= HExpr.ncsfs(subs)
    walk = HExpr.walk_csfs(subs)
    @test size(walk) == (4, 4)
    walk_rows = Set(eachrow(walk))
    @test all(r -> r in walk_rows, eachrow(inp))

    # SOCI: walk_csfs/walk_ncsfs equal csfs/ncsfs (no input set)
    @test HExpr.walk_ncsfs(full) == HExpr.ncsfs(full)
    @test HExpr.walk_csfs(full)  == HExpr.csfs(full)

    # Invalid step value caught by C library
    bad = UInt8[0 1 2 4]  # 1x4 matrix, step=4 is invalid
    @test_throws Exception HExpr.Drt(bad, 4)

    # Shape mismatch caught by Julia validation
    wrong = zeros(UInt8, 2, 3)
    @test_throws Exception HExpr.Drt(wrong, 4)

    HExpr.shutdown()
end
