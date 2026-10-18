# rank1_pair (reduced matrix elements of (a+_p a~_q)^1, since 1.2): the worked
# examples of the spin formulation and two identities, as in
# c/tests/test_rank1.c and python/tests/test_rank1.py.
@testset "rank1_pair" begin
    HExpr.init()

    steps(s) = UInt8[findfirst(==(c), "eudf") - 1 for c in s]

    function R(bra, ket)
        p, q, r = HExpr.rank1_pair(steps(bra), steps(ket))
        Dict((Int(p[k]), Int(q[k])) => r[k] for k in eachindex(p))
    end

    @testset "worked examples" begin
        @test R("udf", "fuu")[(2, 0)] ≈ sqrt(1.5) atol=1e-12   # singlet <- triplet
        @test R("uuu", "udu")[(0, 0)] ≈ sqrt(2.0) atol=1e-12   # quartet <- doublet
        @test R("ud", "uu")[(0, 0)] ≈ -sqrt(1.5) atol=1e-12
        @test R("ud", "uu")[(1, 1)] ≈ sqrt(1.5) atol=1e-12
    end

    # All CSFs (step string, 2S) of norb orbitals and nel electrons.
    function all_csfs(norb, nel)
        out = Tuple{String,Int}[]
        for t in Iterators.product(ntuple(_ -> "eudf", norb)...)
            b, ok = 0, true
            for c in t
                b += c == 'u' ? 1 : c == 'd' ? -1 : 0
                b < 0 && (ok = false; break)
            end
            n = sum(c == 'f' ? 2 : c == 'e' ? 0 : 1 for c in t)
            ok && n == nel && push!(out, (String(collect(t)), b))
        end
        out
    end

    @testset "total spin and adjoint ($norb, $nel)" for (norb, nel) in
            [(3, 2), (3, 3), (4, 4), (4, 3)]
        cs = all_csfs(norb, nel)
        for (s, b) in cs
            b == 0 && continue
            S = b / 2
            tot = sum(v for ((p, q), v) in R(s, s) if p == q)
            @test tot ≈ -sqrt(2) * sqrt(S * (S + 1) * (2S + 1)) atol=1e-12
        end
        for (bra, sb) in cs, (ket, sk) in cs
            abs(sb - sk) > 2 && continue
            fwd, rev = R(bra, ket), R(ket, bra)
            sign = iseven((sb - sk) ÷ 2) ? 1.0 : -1.0
            @test length(fwd) == length(rev)
            for ((p, q), v) in fwd
                @test get(rev, (q, p), 0.0) ≈ sign * v atol=1e-12
            end
        end
    end

    @testset "bad input" begin
        @test_throws ErrorException HExpr.rank1_pair(steps("du"), steps("uu"))
        @test_throws ErrorException HExpr.rank1_pair(steps("uu"), steps("uuu"))
    end
end
