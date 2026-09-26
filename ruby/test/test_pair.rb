require_relative "helper"
require "set"

class TestPair < Minitest::Test
  CSFS_4O4E = [
    [0, 0, 3, 3], [0, 1, 2, 3], [0, 1, 3, 2], [0, 3, 0, 3],
    [0, 3, 1, 2], [0, 3, 3, 0], [1, 0, 2, 3], [1, 0, 3, 2],
    [1, 1, 2, 2], [1, 2, 0, 3], [1, 2, 1, 2], [1, 2, 3, 0],
    [1, 3, 0, 2], [1, 3, 2, 0], [3, 0, 0, 3], [3, 0, 1, 2],
    [3, 0, 3, 0], [3, 1, 0, 2], [3, 1, 2, 0], [3, 3, 0, 0],
  ].freeze

  def setup
    HExpr.init
    @drt = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                           ncore: 0, nval: 4, next: 0)
    @nc = @drt.ncsfs          # 20
    @csfs = CSFS_4O4E
  end

  def teardown
    HExpr.shutdown
  end

  # ----------------------------------------------------------------
  # Helpers
  # ----------------------------------------------------------------

  def triangle(a, b)
    hi, lo = [a, b].max, [a, b].min
    hi * (hi + 1) / 2 + lo
  end

  # Sweep all (bra, ket) pairs; collect (ij, pqrs, coef) with first-wins dedup on (ij, pqrs).
  def aggregate_all_pairs(drt, which: :full)
    n = drt.ncsfs
    seen = {}
    (0...n).each do |bra|
      (0...n).each do |ket|
        ij_v, pq_v, co_v = drt.pair_eval(bra, ket, which: which)
        ij_v.zip(pq_v, co_v).each do |ij, pq, co|
          key = [ij, pq]
          seen[key] ||= co
        end
      end
    end
    ks = seen.keys.sort
    [ks.map { |k| k[0] }, ks.map { |k| k[1] }, ks.map { |k| seen[k] }]
  end

  def sort_block(bp, kp, pq)
    bp.zip(kp, pq).each_with_index.sort_by { |(b, k, p), _i| [b, k, p] }.map(&:last)
  end

  # ----------------------------------------------------------------
  # pair_count / pair_max_terms / pair_eval basics
  # ----------------------------------------------------------------

  def test_pair_count_non_negative
    assert @drt.pair_count(0, 1) >= 0
  end

  def test_pair_max_terms_positive
    assert @drt.pair_max_terms(which: :full) > 0
  end

  def test_pair_max_terms_one_le_full
    assert @drt.pair_max_terms(which: :one) <= @drt.pair_max_terms(which: :full)
  end

  def test_pair_eval_types
    ij, pqrs, coef = @drt.pair_eval(0, 1)
    assert_instance_of Array, ij
    assert_instance_of Array, pqrs
    assert_instance_of Array, coef
    assert ij.all?   { |v| v.is_a?(Integer) }
    assert pqrs.all? { |v| v.is_a?(Integer) }
    assert coef.all? { |v| v.is_a?(Float) }
  end

  def test_pair_eval_length_matches_count
    cnt = @drt.pair_count(0, 1)
    ij, pqrs, coef = @drt.pair_eval(0, 1)
    assert_equal cnt, ij.length
    assert_equal cnt, pqrs.length
    assert_equal cnt, coef.length
  end

  def test_pair_eval_diagonal_has_terms
    total = (0...@nc).sum { |i| @drt.pair_eval(i, i).first.length }
    assert total > 0
  end

  def test_pair_eval_one_subset_of_full
    ij_f, pq_f, _ = @drt.pair_eval(0, 1, which: :full)
    ij_o, pq_o, _ = @drt.pair_eval(0, 1, which: :one)
    assert ij_o.length <= ij_f.length
    full_keys = Set.new(ij_f.zip(pq_f))
    ij_o.zip(pq_o).each { |k| assert full_keys.include?(k) }
  end

  def test_pair_eval_ij_equals_triangle
    (0..2).each do |bra|
      (0..2).each do |ket|
        ij_v, _, _ = @drt.pair_eval(bra, ket)
        expected = triangle(bra, ket)
        ij_v.each { |v| assert_equal expected, v }
      end
    end
  end

  # ----------------------------------------------------------------
  # Aggregate: pair_eval over all pairs matches read_full
  # ----------------------------------------------------------------

  def test_pair_eval_aggregate_matches_expression
    ij_agg, pq_agg, co_agg = aggregate_all_pairs(@drt)

    expr = HExpr::Expression.new(@drt)
    ij_ref, pq_ref, co_ref = expr.read_full

    ord = ij_ref.zip(pq_ref).each_with_index.sort_by { |(ij, pq), _| [ij, pq] }.map(&:last)
    ij_ref_s  = ord.map { |i| ij_ref[i] }
    pq_ref_s  = ord.map { |i| pq_ref[i] }
    co_ref_s  = ord.map { |i| co_ref[i] }

    assert_equal ij_ref_s.length, ij_agg.length
    assert_equal ij_ref_s,  ij_agg
    assert_equal pq_ref_s,  pq_agg
    co_agg.zip(co_ref_s).each { |a, b| assert_in_delta b, a, 1e-10 }
  end

  # ----------------------------------------------------------------
  # block_count / block_max_terms / block_eval basics
  # ----------------------------------------------------------------

  def test_block_count_non_negative
    assert @drt.block_count([0, 1], [2, 3]) >= 0
  end

  def test_block_max_terms_ge_count
    cnt = @drt.block_count([0, 1], [2, 3])
    assert @drt.block_max_terms(2, 2) >= cnt
  end

  def test_block_eval_types
    ij, pqrs, coef, bp, kp = @drt.block_eval([0, 1], [2, 3])
    assert_instance_of Array, ij
    assert_instance_of Array, pqrs
    assert_instance_of Array, coef
    assert_instance_of Array, bp
    assert_instance_of Array, kp
    assert ij.all?   { |v| v.is_a?(Integer) }
    assert pqrs.all? { |v| v.is_a?(Integer) }
    assert coef.all? { |v| v.is_a?(Float) }
    assert bp.all?   { |v| v.is_a?(Integer) }
    assert kp.all?   { |v| v.is_a?(Integer) }
  end

  def test_block_eval_length_matches_count
    cnt = @drt.block_count([0, 1], [2, 3])
    ij, pqrs, coef, bp, kp = @drt.block_eval([0, 1], [2, 3])
    [ij, pqrs, coef, bp, kp].each { |a| assert_equal cnt, a.length }
  end

  def test_block_eval_bra_pos_ket_pos_semantics
    bra_list = [5, 9, 12]
    ket_list = [2, 7]
    _, _, _, bp, kp = @drt.block_eval(bra_list, ket_list)
    assert bp.all? { |v| v >= 0 && v < bra_list.length }
    assert kp.all? { |v| v >= 0 && v < ket_list.length }
    # positions must be in {0,1,2}, not in {5,9,12}
    unless bp.empty?
      assert (Set.new(bp) - Set.new([0, 1, 2])).empty?
    end
  end

  def test_block_eval_matches_per_pair_pair_eval
    bra_list = [0, 1]
    ket_list = [2, 3]
    ij_b, pq_b, co_b, bp, kp = @drt.block_eval(bra_list, ket_list)

    bra_list.each_with_index do |bra, ip|
      ket_list.each_with_index do |ket, jp|
        ij_p, pq_p, co_p = @drt.pair_eval(bra, ket)
        mask = bp.zip(kp).each_with_index.select { |(b, k), _| b == ip && k == jp }.map(&:last)
        pq_b_sel = mask.map { |i| pq_b[i] }
        ij_b_sel = mask.map { |i| ij_b[i] }
        co_b_sel = mask.map { |i| co_b[i] }

        ord_b = pq_b_sel.each_with_index.sort_by { |v, _| v }.map(&:last)
        ord_p = pq_p.each_with_index.sort_by { |v, _| v }.map(&:last)

        assert_equal ord_p.map { |i| pq_p[i] }, ord_b.map { |i| pq_b_sel[i] }
        assert_equal ord_p.map { |i| ij_p[i] }, ord_b.map { |i| ij_b_sel[i] }
        ord_b.zip(ord_p) { |ib, ip2|
          assert_in_delta co_p[ip2], co_b_sel[ib], 1e-12
        }
      end
    end
  end

  def test_block_eval_empty_bra
    ij, pqrs, coef, bp, kp = @drt.block_eval([], [0, 1])
    assert_equal 0, ij.length
  end

  def test_block_eval_empty_ket
    ij, pqrs, coef, bp, kp = @drt.block_eval([0, 1], [])
    assert_equal 0, ij.length
  end

  # ----------------------------------------------------------------
  # subspace_count / subspace_max_terms / subspace_eval
  # ----------------------------------------------------------------

  def test_subspace_count_non_negative
    assert @drt.subspace_count([0, 1, 2]) >= 0
  end

  def test_subspace_max_terms_ge_count
    idx = [0, 1, 2]
    assert @drt.subspace_max_terms(idx.length) >= @drt.subspace_count(idx)
  end

  def test_subspace_eval_matches_block_eval
    idx = [0, 3, 7]
    ij_s, pq_s, co_s, bp_s, kp_s = @drt.subspace_eval(idx)
    ij_b, pq_b, co_b, bp_b, kp_b = @drt.block_eval(idx, idx)
    assert_equal ij_b, ij_s
    assert_equal pq_b, pq_s
    assert_equal co_b, co_s
    assert_equal bp_b, bp_s
    assert_equal kp_b, kp_s
  end

  def test_subspace_eval_pos_within_range
    idx = [0, 1, 2, 5]
    _, _, _, bp, kp = @drt.subspace_eval(idx)
    assert bp.all? { |v| v >= 0 && v < idx.length }
    assert kp.all? { |v| v >= 0 && v < idx.length }
  end

  # ----------------------------------------------------------------
  # High-level csf_pair_eval
  # ----------------------------------------------------------------

  def test_csf_pair_eval_matches_low_level
    bra, ket = 0, 3
    bra_steps = @csfs[bra]
    ket_steps = @csfs[ket]

    ij_hi, pq_hi, co_hi = HExpr.csf_pair_eval(bra_steps, ket_steps)
    ij_lo, pq_lo, co_lo = @drt.pair_eval(bra, ket)

    ord_hi = pq_hi.each_with_index.sort_by { |v, _| v }.map(&:last)
    ord_lo = pq_lo.each_with_index.sort_by { |v, _| v }.map(&:last)

    assert_equal ij_lo.length, ij_hi.length
    assert_equal ord_lo.map { |i| pq_lo[i] }, ord_hi.map { |i| pq_hi[i] }
    ord_hi.zip(ord_lo) { |ih, il|
      assert_in_delta co_lo[il], co_hi[ih], 1e-12
    }
  end

  def test_csf_pair_eval_diagonal
    idx = 5
    steps = @csfs[idx]
    ij_hi, pq_hi, co_hi = HExpr.csf_pair_eval(steps, steps)
    ij_lo, pq_lo, co_lo = @drt.pair_eval(idx, idx)

    ord_hi = pq_hi.each_with_index.sort_by { |v, _| v }.map(&:last)
    ord_lo = pq_lo.each_with_index.sort_by { |v, _| v }.map(&:last)

    assert_equal ij_lo.length, ij_hi.length
    assert_equal ord_lo.map { |i| pq_lo[i] }, ord_hi.map { |i| pq_hi[i] }
    ord_hi.zip(ord_lo) { |ih, il|
      assert_in_delta co_lo[il], co_hi[ih], 1e-12
    }
  end

  def test_csf_pair_eval_reordering
    # bra,ket whose 2-CSF sub-DRT reorders (canonical order != input order).
    # The eval API takes input-order indices, so csf_pair_eval passes rows 0
    # and 1 and the library maps them onto the reordered canonical walk.
    # (csfs[1], csfs[0]) swaps; the assertion below keeps the fixture honest.
    bra, ket = 1, 0
    sub = HExpr::Drt.from_csfset([@csfs[bra], @csfs[ket]], @csfs[bra].size)
    refute_equal sub.csfs, sub.walk_csfs,
                 "fixture no longer reorders; pick a different (bra,ket)"

    ij_hi, pq_hi, co_hi = HExpr.csf_pair_eval(@csfs[bra], @csfs[ket])
    ij_lo, pq_lo, co_lo = @drt.pair_eval(bra, ket)

    ord_hi = pq_hi.each_with_index.sort_by { |v, _| v }.map(&:last)
    ord_lo = pq_lo.each_with_index.sort_by { |v, _| v }.map(&:last)

    assert_equal ij_lo.length, ij_hi.length
    assert_equal ord_lo.map { |i| pq_lo[i] }, ord_hi.map { |i| pq_hi[i] }
    ord_hi.zip(ord_lo) { |ih, il|
      assert_in_delta co_lo[il], co_hi[ih], 1e-12
    }
  end

  def test_csf_pair_eval_return_types
    ij, pqrs, coef = HExpr.csf_pair_eval(@csfs[0], @csfs[1])
    assert ij.all?   { |v| v.is_a?(Integer) }
    assert pqrs.all? { |v| v.is_a?(Integer) }
    assert coef.all? { |v| v.is_a?(Float) }
  end

  # ----------------------------------------------------------------
  # High-level csf_block_eval
  # ----------------------------------------------------------------

  def test_csf_block_eval_matches_low_level
    bra_list = [0, 2]
    ket_list = [1, 3, 4]
    bra_csfs = bra_list.map { |i| @csfs[i] }
    ket_csfs = ket_list.map { |i| @csfs[i] }

    ij_hi, pq_hi, co_hi, bp_hi, kp_hi = HExpr.csf_block_eval(bra_csfs, ket_csfs)
    ij_lo, pq_lo, co_lo, bp_lo, kp_lo = @drt.block_eval(bra_list, ket_list)

    assert_equal ij_lo.length, ij_hi.length
    unless ij_hi.empty?
      ord_hi = sort_block(bp_hi, kp_hi, pq_hi)
      ord_lo = sort_block(bp_lo, kp_lo, pq_lo)
      assert_equal ord_lo.map { |i| bp_lo[i] }, ord_hi.map { |i| bp_hi[i] }
      assert_equal ord_lo.map { |i| kp_lo[i] }, ord_hi.map { |i| kp_hi[i] }
      assert_equal ord_lo.map { |i| pq_lo[i] }, ord_hi.map { |i| pq_hi[i] }
      ord_hi.zip(ord_lo) { |ih, il|
        assert_in_delta co_lo[il], co_hi[ih], 1e-12
      }
    end
  end

  def test_csf_block_eval_bra_pos_within_range
    bra_csfs = @csfs[0, 3]
    ket_csfs = @csfs[3, 2]
    _, _, _, bp, kp = HExpr.csf_block_eval(bra_csfs, ket_csfs)
    assert bp.all? { |v| v >= 0 && v < 3 }
    assert kp.all? { |v| v >= 0 && v < 2 }
  end

  # ----------------------------------------------------------------
  # High-level csf_subspace_eval
  # ----------------------------------------------------------------

  def test_csf_subspace_eval_matches_low_level
    idx = [0, 1, 4]
    sub_csfs = idx.map { |i| @csfs[i] }

    ij_hi, pq_hi, co_hi, bp_hi, kp_hi = HExpr.csf_subspace_eval(sub_csfs)
    ij_lo, pq_lo, co_lo, bp_lo, kp_lo = @drt.subspace_eval(idx)

    assert_equal ij_lo.length, ij_hi.length
    unless ij_hi.empty?
      ord_hi = sort_block(bp_hi, kp_hi, pq_hi)
      ord_lo = sort_block(bp_lo, kp_lo, pq_lo)
      assert_equal ord_lo.map { |i| bp_lo[i] }, ord_hi.map { |i| bp_hi[i] }
      assert_equal ord_lo.map { |i| kp_lo[i] }, ord_hi.map { |i| kp_hi[i] }
      assert_equal ord_lo.map { |i| pq_lo[i] }, ord_hi.map { |i| pq_hi[i] }
      ord_hi.zip(ord_lo) { |ih, il|
        assert_in_delta co_lo[il], co_hi[ih], 1e-12
      }
    end
  end

  # ----------------------------------------------------------------
  # Error paths
  # ----------------------------------------------------------------

  def test_pair_count_invalid_which
    assert_raises(ArgumentError) { @drt.pair_count(0, 1, which: :bogus) }
  end

  def test_pair_eval_invalid_which
    assert_raises(ArgumentError) { @drt.pair_eval(0, 1, which: :two) }
  end

  def test_block_eval_invalid_which
    assert_raises(ArgumentError) { @drt.block_eval([0], [1], which: :all) }
  end

  def test_subspace_eval_invalid_which
    assert_raises(ArgumentError) { @drt.subspace_eval([0], which: :half) }
  end

  def test_csf_pair_eval_invalid_which
    assert_raises(ArgumentError) { HExpr.csf_pair_eval(@csfs[0], @csfs[1], which: :bad) }
  end

  def test_pair_count_out_of_range
    assert_raises(RuntimeError) { @drt.pair_count(9999, 0) }
  end

  def test_pair_eval_out_of_range
    assert_raises(RuntimeError) { @drt.pair_eval(0, 9999) }
  end

  def test_block_count_out_of_range
    assert_raises(RuntimeError) { @drt.block_count([0, 9999], [1]) }
  end

  def test_csf_pair_eval_length_mismatch
    assert_raises(ArgumentError) { HExpr.csf_pair_eval([0, 1, 2, 3], [0, 1, 2]) }
  end

  def test_csf_block_eval_norb_mismatch
    bra_csfs = [[0, 1, 2, 3], [1, 0, 2, 3]]
    ket_csfs = [[0, 1, 2], [1, 0, 2]]
    assert_raises(ArgumentError) { HExpr.csf_block_eval(bra_csfs, ket_csfs) }
  end
end
