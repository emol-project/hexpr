# Binding-completeness wrappers: csf_steps, step_parse, drt_node_info, and
# the *_show* functions. Reference case: s_full_4o4e (norb=4, ncsfs=20).
@testset "binding completeness" begin
    HExpr.init()

    drt = HExpr.Drt("soci"; nve=4, spin_x2=0, ncore=0, nval=4, next=0)
    @test HExpr.ncsfs(drt) == 20

    # step_parse("fude", 4) -> [3,1,2,0]
    @test HExpr.step_parse("fude", 4) == UInt8[3, 1, 2, 0]

    # csf_steps matches the csfs matrix and round-trips through csf_index
    rows = HExpr.csfs(drt)
    for i in 0:(Int(HExpr.ncsfs(drt)) - 1)
        steps = HExpr.csf_steps(drt, i)
        @test steps == rows[i + 1, :]
        @test HExpr.csf_index(drt, steps) == i
    end

    # drt_node_info against known nodes
    @test HExpr.drt_node_info(drt, 0) == (orb = 0, nve = 0, is = 0)
    @test HExpr.drt_node_info(drt, Int(HExpr.nnodes(drt)) - 1) ==
          (orb = 4, nve = 4, is = 0)

    # *_show*: callable, return nothing, do not throw
    expr = HExpr.Expression(drt)
    @test HExpr.drt_show(drt) === nothing
    @test HExpr.expr_show_full(expr) === nothing
    @test HExpr.expr_show_one(expr) === nothing

    # ---- csf2det: the determinant expansion ----
    #
    # What is and is not asserted: the overall sign of a CSF is a
    # convention, not a measured quantity, and both <I|J> and
    # U H_det U^T are invariant under flipping it.  So the count, the
    # masks, the magnitudes and the RELATIVE signs are asserted here,
    # and no individual coefficient sign ever is.

    # Closed-shell 'ffee': one determinant, orbitals 0 and 1 doubly
    # occupied.  With no open shells there is no coupling and no sign
    # convention in play, so the masks are forced.
    idx = HExpr.csf_index(drt, HExpr.step_parse("ffee", 4))
    @test HExpr.csf_ndets(drt, idx, 0) == 1
    a, b, c = HExpr.csf_to_dets(drt, idx, 0)
    @test eltype(a) === UInt64
    @test size(a) == (1, 1)
    @test size(b) == (1, 1)
    @test a[1, 1] == 0x03
    @test b[1, 1] == 0x03
    @test isapprox(abs(c[1]), 1.0)

    # Open-shell singlet 'udfe': two determinants, equal magnitudes, the
    # SAME relative sign, and each the other's spin flip.
    idx = HExpr.csf_index(drt, HExpr.step_parse("udfe", 4))
    @test HExpr.csf_ndets(drt, idx, 0) == 2
    a, b, c = HExpr.csf_to_dets(drt, idx, 0)
    @test size(a) == (2, 1)
    @test all(isapprox.(abs.(c), 1 / sqrt(2)))
    @test isapprox(sum(c .^ 2), 1.0)
    @test c[1] * c[2] > 0
    @test a[1, 1] == b[2, 1]
    @test b[1, 1] == a[2, 1]
    @test a[1, 1] != a[2, 1]

    # ms_x2 is a parameter, not an assumption: wrong parity and |Ms| > S
    # are both errors.
    @test_throws ErrorException HExpr.csf_ndets(drt, idx, 1)
    @test_throws ErrorException HExpr.csf_ndets(drt, idx, 2)

    # Multi-word determinants.  Everything above runs at norb = 4, where
    # nwords == 1 and the second mask index is never used.  This walks a
    # 70-orbital DRT and reconstructs the occupation from the masks, so a
    # wrong word index or stride fails instead of answering differently.
    wide = HExpr.Drt("foci"; nve=2, spin_x2=0, ncore=1, nval=1, next=68)
    @test HExpr.norb(wide) == 70
    wide_norb = Int(HExpr.norb(wide))
    nwords = (wide_norb + 63) >> 6
    @test nwords == 2

    step_occ = Dict(UInt8(0) => 0, UInt8(1) => 1, UInt8(2) => 1, UInt8(3) => 2)
    saw_high_word = false
    for i in 0:(Int(HExpr.ncsfs(wide)) - 1)
        a, b, c = HExpr.csf_to_dets(wide, i, 0)
        @test size(a) == (length(c), nwords)
        @test isapprox(sum(c .^ 2), 1.0)
        want = [step_occ[s] for s in HExpr.csf_steps(wide, i)]
        for k in 1:length(c)
            got = [Int((a[k, (p >> 6) + 1] >> (p & 63)) & 0x01) +
                   Int((b[k, (p >> 6) + 1] >> (p & 63)) & 0x01)
                   for p in 0:(wide_norb - 1)]
            @test got == want
            if a[k, 2] != 0 || b[k, 2] != 0
                saw_high_word = true
            end
        end
    end
    @test saw_high_word

    HExpr.shutdown()
end
