require_relative "helper"

# HExpr.rank1_pair (reduced matrix elements of (a+_p a~_q)^1, since 1.2): the
# worked examples of the spin formulation and two identities, as in
# c/tests/test_rank1.c and python/tests/test_rank1.py.
class TestRank1 < Minitest::Test
  def setup
    HExpr.init
  end

  def teardown
    HExpr.shutdown
  end

  def steps(s)
    s.chars.map { |c| "eudf".index(c) }
  end

  def r(bra, ket)
    p, q, v = HExpr.rank1_pair(steps(bra), steps(ket))
    p.zip(q, v).to_h { |a, b, c| [[a, b], c] }
  end

  def test_worked_examples
    assert_in_delta Math.sqrt(1.5), r("udf", "fuu")[[2, 0]], 1e-12  # singlet <- triplet
    assert_in_delta Math.sqrt(2.0), r("uuu", "udu")[[0, 0]], 1e-12  # quartet <- doublet
    assert_in_delta(-Math.sqrt(1.5), r("ud", "uu")[[0, 0]], 1e-12)
    assert_in_delta Math.sqrt(1.5), r("ud", "uu")[[1, 1]], 1e-12
  end

  # All CSFs [step string, 2S] of norb orbitals and nel electrons.
  def all_csfs(norb, nel)
    %w[e u d f].repeated_permutation(norb).filter_map do |t|
      b = 0
      ok = t.all? { |c| b += { "u" => 1, "d" => -1 }.fetch(c, 0); b >= 0 }
      n = t.sum { |c| { "e" => 0, "u" => 1, "d" => 1, "f" => 2 }[c] }
      [t.join, b] if ok && n == nel
    end
  end

  def test_total_spin_and_adjoint
    [[3, 2], [3, 3], [4, 4], [4, 3]].each do |norb, nel|
      cs = all_csfs(norb, nel)
      cs.each do |s, b|
        next if b.zero?
        sp = b / 2.0
        tot = r(s, s).sum { |(p, q), v| p == q ? v : 0.0 }
        assert_in_delta(-Math.sqrt(2) * Math.sqrt(sp * (sp + 1) * (2 * sp + 1)), tot, 1e-12)
      end
      cs.product(cs).each do |(bra, sb), (ket, sk)|
        next if (sb - sk).abs > 2
        fwd = r(bra, ket)
        rev = r(ket, bra)
        sign = ((sb - sk) / 2).even? ? 1.0 : -1.0
        assert_equal fwd.size, rev.size
        fwd.each { |(p, q), v| assert_in_delta sign * v, rev.fetch([q, p], 0.0), 1e-12 }
      end
    end
  end

  def test_rejects_bad_input
    assert_raises(RuntimeError) { HExpr.rank1_pair(steps("du"), steps("uu")) }
    assert_raises(ArgumentError) { HExpr.rank1_pair(steps("uu"), steps("uuu")) }
  end
end
