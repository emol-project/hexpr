@testset "pair / block / subspace" begin
    HExpr.init()

    drt = HExpr.Drt("soci"; nve=4, spin_x2=0, ncore=0, nval=4, next=0)
    nc  = Int(HExpr.ncsfs(drt))   # 20
    csfs_mat = HExpr.csfs(drt)    # (20, 4) UInt8

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    function triangle(a, b)
        hi, lo = max(a, b), min(a, b)
        hi * (hi + 1) ÷ 2 + lo
    end

    # Aggregate pair_eval over all (bra, ket), first-wins dedup on (ij, pqrs).
    function aggregate_all_pairs(d; which=:full)
        n = Int(HExpr.ncsfs(d))
        seen = Dict{Tuple{Int32,Int32}, Float64}()
        for bra in 0:n-1, ket in 0:n-1
            ij_v, pq_v, co_v = HExpr.pair_eval(d, bra, ket; which=which)
            for k in eachindex(ij_v)
                key = (ij_v[k], pq_v[k])
                haskey(seen, key) || (seen[key] = co_v[k])
            end
        end
        ks = sort!(collect(keys(seen)))
        Int32[k[1] for k in ks], Int32[k[2] for k in ks], Float64[seen[k] for k in ks]
    end

    # Sort block/subspace result by (bra_pos, ket_pos, pqrs).
    sort_block(bp, kp, pq) = sortperm(collect(zip(bp, kp, pq)))

    # ------------------------------------------------------------------
    # pair_count / pair_max_terms / pair_eval - basics
    # ------------------------------------------------------------------

    @testset "pair basics" begin
        @test HExpr.pair_count(drt, 0, 1) >= 0
        @test HExpr.pair_max_terms(drt; which=:full) > 0
        @test HExpr.pair_max_terms(drt; which=:one) <= HExpr.pair_max_terms(drt; which=:full)

        ij, pqrs, coef = HExpr.pair_eval(drt, 0, 1)
        @test eltype(ij)   == Int32
        @test eltype(pqrs) == Int32
        @test eltype(coef) == Float64

        cnt = HExpr.pair_count(drt, 0, 1)
        @test length(ij) == length(pqrs) == length(coef) == cnt

        # diagonal should produce some terms in total
        total_diag = sum(length(HExpr.pair_eval(drt, i, i)[1]) for i in 0:nc-1)
        @test total_diag > 0

        # :one is a subset of :full  (same pqrs keys)
        ij_f, pq_f, _ = HExpr.pair_eval(drt, 0, 1; which=:full)
        ij_o, pq_o, _ = HExpr.pair_eval(drt, 0, 1; which=:one)
        @test length(ij_o) <= length(ij_f)
        full_keys = Set(zip(ij_f, pq_f))
        for k in eachindex(ij_o)
            @test (ij_o[k], pq_o[k]) in full_keys
        end

        # ij == triangle(bra, ket) for all returned terms
        for bra in 0:2, ket in 0:2
            ij_v, _, _ = HExpr.pair_eval(drt, bra, ket)
            expected = Int32(triangle(bra, ket))
            @test all(==(expected), ij_v)
        end
    end

    # ------------------------------------------------------------------
    # Aggregate test: pair_eval over all pairs matches read_full
    # ------------------------------------------------------------------

    @testset "pair_eval aggregate matches Expression" begin
        ij_agg, pq_agg, co_agg = aggregate_all_pairs(drt)

        expr = HExpr.Expression(drt)
        ij_ref, pq_ref, co_ref = HExpr.read_full(expr)

        # sort reference by (ij, pqrs)
        ord = sortperm(collect(zip(ij_ref, pq_ref)))
        ij_ref_s = ij_ref[ord]
        pq_ref_s = pq_ref[ord]
        co_ref_s = co_ref[ord]

        @test length(ij_agg) == length(ij_ref_s)
        @test ij_agg  == ij_ref_s
        @test pq_agg  == pq_ref_s
        @test all(abs.(co_agg .- co_ref_s) .< 1e-10)

        finalize(expr)
    end

    # ------------------------------------------------------------------
    # block_count / block_max_terms / block_eval
    # ------------------------------------------------------------------

    @testset "block basics" begin
        bra_list = Int32[0, 1]
        ket_list = Int32[2, 3]

        @test HExpr.block_count(drt, bra_list, ket_list) >= 0
        @test HExpr.block_max_terms(drt, 2, 2) >= HExpr.block_count(drt, bra_list, ket_list)

        ij, pqrs, coef, bp, kp = HExpr.block_eval(drt, bra_list, ket_list)
        @test eltype(ij)   == Int32
        @test eltype(pqrs) == Int32
        @test eltype(coef) == Float64
        @test eltype(bp)   == Int32
        @test eltype(kp)   == Int32

        cnt = HExpr.block_count(drt, bra_list, ket_list)
        @test length(ij) == length(pqrs) == length(coef) == length(bp) == length(kp) == cnt
    end

    @testset "block_eval bra_pos / ket_pos semantics" begin
        bra_list = Int32[5, 9, 12]
        ket_list = Int32[2, 7]
        ij, pqrs, coef, bp, kp = HExpr.block_eval(drt, bra_list, ket_list)

        # positions are 0-based indices into the input arrays, not raw CSF indices
        @test all(0 .<= bp .< length(bra_list))
        @test all(0 .<= kp .< length(ket_list))
        # specifically: values must be in {0,1,2}, not in {5,9,12}
        if !isempty(bp)
            @test Set(bp) ⊆ Set(Int32[0, 1, 2])
        end
    end

    @testset "block_eval matches per-pair pair_eval" begin
        bra_list = Int32[0, 1]
        ket_list = Int32[2, 3]
        ij_b, pq_b, co_b, bp, kp = HExpr.block_eval(drt, bra_list, ket_list)

        for (ip, bra) in enumerate(bra_list), (jp, ket) in enumerate(ket_list)
            ij_p, pq_p, co_p = HExpr.pair_eval(drt, bra, ket)
            mask = (bp .== Int32(ip - 1)) .& (kp .== Int32(jp - 1))
            # sort both by pqrs
            ord_b = sortperm(pq_b[mask])
            ord_p = sortperm(pq_p)
            @test pq_b[mask][ord_b] == pq_p[ord_p]
            @test ij_b[mask][ord_b] == ij_p[ord_p]
            @test all(abs.(co_b[mask][ord_b] .- co_p[ord_p]) .< 1e-12)
        end
    end

    @testset "block_eval empty arrays" begin
        ij, pqrs, coef, bp, kp = HExpr.block_eval(drt, Int32[], Int32[0, 1])
        @test length(ij) == 0
        ij2, pqrs2, coef2, bp2, kp2 = HExpr.block_eval(drt, Int32[0, 1], Int32[])
        @test length(ij2) == 0
    end

    # ------------------------------------------------------------------
    # subspace_count / subspace_max_terms / subspace_eval
    # ------------------------------------------------------------------

    @testset "subspace basics" begin
        idx = Int32[0, 1, 2]
        @test HExpr.subspace_count(drt, idx) >= 0
        @test HExpr.subspace_max_terms(drt, length(idx)) >= HExpr.subspace_count(drt, idx)
    end

    @testset "subspace_eval matches block_eval" begin
        idx = Int32[0, 3, 7]
        ij_s, pq_s, co_s, bp_s, kp_s = HExpr.subspace_eval(drt, idx)
        ij_b, pq_b, co_b, bp_b, kp_b = HExpr.block_eval(drt, idx, idx)
        @test ij_s == ij_b
        @test pq_s == pq_b
        @test co_s == co_b
        @test bp_s == bp_b
        @test kp_s == kp_b
    end

    @testset "subspace_eval pos within range" begin
        idx = Int32[0, 1, 2, 5]
        _, _, _, bp, kp = HExpr.subspace_eval(drt, idx)
        @test all(0 .<= bp .< length(idx))
        @test all(0 .<= kp .< length(idx))
    end

    # ------------------------------------------------------------------
    # High-level CSF step-code functions
    # ------------------------------------------------------------------

    @testset "csf_pair_eval matches low-level" begin
        bra, ket = 0, 3
        bra_steps = csfs_mat[bra+1, :]   # Julia is 1-indexed
        ket_steps = csfs_mat[ket+1, :]

        ij_hi, pq_hi, co_hi = HExpr.csf_pair_eval(bra_steps, ket_steps)
        ij_lo, pq_lo, co_lo = HExpr.pair_eval(drt, bra, ket)

        ord_hi = sortperm(pq_hi)
        ord_lo = sortperm(pq_lo)
        @test length(ij_hi) == length(ij_lo)
        @test pq_hi[ord_hi] == pq_lo[ord_lo]
        @test all(abs.(co_hi[ord_hi] .- co_lo[ord_lo]) .< 1e-12)
    end

    @testset "csf_pair_eval diagonal" begin
        idx = 5
        steps = csfs_mat[idx+1, :]
        ij_hi, pq_hi, co_hi = HExpr.csf_pair_eval(steps, steps)
        ij_lo, pq_lo, co_lo = HExpr.pair_eval(drt, idx, idx)
        ord_hi = sortperm(pq_hi)
        ord_lo = sortperm(pq_lo)
        @test length(ij_hi) == length(ij_lo)
        @test pq_hi[ord_hi] == pq_lo[ord_lo]
        @test all(abs.(co_hi[ord_hi] .- co_lo[ord_lo]) .< 1e-12)
    end

    @testset "csf_pair_eval reordering" begin
        # bra,ket whose 2-CSF sub-DRT reorders (canonical order != input).
        # The eval API takes input-order indices, so csf_pair_eval passes rows
        # 0 and 1 and the library maps them onto the reordered canonical walk.
        # (csfs[1],csfs[0]) swaps; the assertion below keeps the fixture honest.
        bra, ket = 1, 0
        sub = HExpr.Drt(csfs_mat[[bra+1, ket+1], :], size(csfs_mat, 2))
        @test HExpr.csfs(sub) != HExpr.walk_csfs(sub)   # fixture really reorders

        ij_hi, pq_hi, co_hi = HExpr.csf_pair_eval(csfs_mat[bra+1, :], csfs_mat[ket+1, :])
        ij_lo, pq_lo, co_lo = HExpr.pair_eval(drt, bra, ket)
        ord_hi = sortperm(pq_hi)
        ord_lo = sortperm(pq_lo)
        @test length(ij_hi) == length(ij_lo)
        @test pq_hi[ord_hi] == pq_lo[ord_lo]
        @test all(abs.(co_hi[ord_hi] .- co_lo[ord_lo]) .< 1e-12)
    end

    @testset "csf_pair_eval return types" begin
        ij, pqrs, coef = HExpr.csf_pair_eval(csfs_mat[1, :], csfs_mat[2, :])
        @test eltype(ij)   == Int32
        @test eltype(pqrs) == Int32
        @test eltype(coef) == Float64
    end

    @testset "csf_block_eval matches low-level" begin
        bra_list = Int32[0, 2]
        ket_list = Int32[1, 3, 4]
        bra_csfs = csfs_mat[bra_list.+1, :]   # 0-based -> 1-based
        ket_csfs = csfs_mat[ket_list.+1, :]

        ij_hi, pq_hi, co_hi, bp_hi, kp_hi = HExpr.csf_block_eval(bra_csfs, ket_csfs)
        ij_lo, pq_lo, co_lo, bp_lo, kp_lo = HExpr.block_eval(drt, bra_list, ket_list)

        @test length(ij_hi) == length(ij_lo)
        if !isempty(ij_hi)
            ord_hi = sort_block(bp_hi, kp_hi, pq_hi)
            ord_lo = sort_block(bp_lo, kp_lo, pq_lo)
            @test bp_hi[ord_hi] == bp_lo[ord_lo]
            @test kp_hi[ord_hi] == kp_lo[ord_lo]
            @test pq_hi[ord_hi] == pq_lo[ord_lo]
            @test all(abs.(co_hi[ord_hi] .- co_lo[ord_lo]) .< 1e-12)
        end
    end

    @testset "csf_block_eval bra_pos within range" begin
        bra_csfs = csfs_mat[1:3, :]
        ket_csfs = csfs_mat[4:5, :]
        _, _, _, bp, kp = HExpr.csf_block_eval(bra_csfs, ket_csfs)
        @test all(0 .<= bp .< 3)
        @test all(0 .<= kp .< 2)
    end

    @testset "csf_subspace_eval matches low-level" begin
        idx = Int32[0, 1, 4]
        sub_csfs = csfs_mat[idx.+1, :]

        ij_hi, pq_hi, co_hi, bp_hi, kp_hi = HExpr.csf_subspace_eval(sub_csfs)
        ij_lo, pq_lo, co_lo, bp_lo, kp_lo = HExpr.subspace_eval(drt, idx)

        @test length(ij_hi) == length(ij_lo)
        if !isempty(ij_hi)
            ord_hi = sort_block(bp_hi, kp_hi, pq_hi)
            ord_lo = sort_block(bp_lo, kp_lo, pq_lo)
            @test bp_hi[ord_hi] == bp_lo[ord_lo]
            @test kp_hi[ord_hi] == kp_lo[ord_lo]
            @test pq_hi[ord_hi] == pq_lo[ord_lo]
            @test all(abs.(co_hi[ord_hi] .- co_lo[ord_lo]) .< 1e-12)
        end
    end

    # ------------------------------------------------------------------
    # Error paths
    # ------------------------------------------------------------------

    @testset "error paths" begin
        # invalid which
        @test_throws ArgumentError HExpr.pair_count(drt, 0, 1; which=:bogus)
        @test_throws ArgumentError HExpr.pair_eval(drt, 0, 1; which=:two)
        @test_throws ArgumentError HExpr.block_eval(drt, [0], [1]; which=:all)
        @test_throws ArgumentError HExpr.subspace_eval(drt, [0]; which=:half)
        @test_throws ArgumentError HExpr.csf_pair_eval(csfs_mat[1,:], csfs_mat[2,:]; which=:bad)

        # out-of-range indices
        @test_throws Exception HExpr.pair_count(drt, 9999, 0)
        @test_throws Exception HExpr.pair_eval(drt, 0, 9999)
        @test_throws Exception HExpr.block_count(drt, [0, 9999], [1])

        # csf_pair_eval length mismatch
        @test_throws Exception HExpr.csf_pair_eval(UInt8[0,1,2,3], UInt8[0,1,2])

        # csf_block_eval norb mismatch
        @test_throws Exception HExpr.csf_block_eval(
            zeros(UInt8, 2, 4), zeros(UInt8, 2, 5))
    end

    finalize(drt)
    HExpr.shutdown()
end
