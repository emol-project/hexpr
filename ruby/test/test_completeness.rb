require_relative "helper"

# Binding-completeness wrappers: csf_steps, step_parse, node_info, and the
# *_show* methods. Reference case: s_full_4o4e (norb=4, ncsfs=20).
class TestCompleteness < Minitest::Test
  def setup    = HExpr.init
  def teardown = HExpr.shutdown

  def ref_drt
    HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0, ncore: 0, nval: 4, next: 0)
  end

  def test_step_parse_fude
    assert_equal [3, 1, 2, 0], HExpr.step_parse("fude", 4)
  end

  def test_csf_steps_round_trip
    drt  = ref_drt
    rows = drt.csfs
    assert_equal 20, drt.ncsfs
    drt.ncsfs.times do |i|
      steps = drt.csf_steps(i)
      assert_equal rows[i], steps
      assert_equal i, drt.csf_index(steps)
    end
  end

  def test_node_info_known_nodes
    drt = ref_drt
    assert_equal [0, 0, 0], drt.node_info(0)
    assert_equal [4, 4, 0], drt.node_info(drt.nnodes - 1)
  end

  def test_show_methods_callable
    drt  = ref_drt
    expr = HExpr::Expression.new(drt)
    assert_nil drt.show
    assert_nil expr.show_full
    assert_nil expr.show_one
  end

  # --- csf2det: the determinant expansion ---
  #
  # What is and is not asserted: the overall sign of a CSF is a
  # convention, not a measured quantity, and both <I|J> and U H_det U^T
  # are invariant under flipping it. So the count, the masks, the
  # magnitudes and the RELATIVE signs are asserted here, and no
  # individual coefficient sign ever is.

  # Closed-shell "ffee": one determinant, orbitals 0 and 1 doubly
  # occupied. With no open shells there is no coupling and no sign
  # convention in play, so the masks are forced.
  def test_csf2det_closed_shell
    drt = ref_drt
    idx = drt.csf_index(HExpr.step_parse("ffee", 4))
    assert_equal 1, drt.csf_ndets(idx, 0)

    alpha, beta, coef = drt.csf_to_dets(idx, 0)
    assert_equal 1, alpha.size
    assert_equal 1, alpha[0].size          # nwords == 1 at norb == 4
    assert_equal 0b11, alpha[0][0]
    assert_equal 0b11, beta[0][0]
    assert_in_delta 1.0, coef[0].abs, 1e-12
  end

  # Open-shell singlet "udfe": two determinants, equal magnitudes, the
  # SAME relative sign, and each the other's spin flip.
  def test_csf2det_open_shell_singlet
    drt = ref_drt
    idx = drt.csf_index(HExpr.step_parse("udfe", 4))
    assert_equal 2, drt.csf_ndets(idx, 0)

    alpha, beta, coef = drt.csf_to_dets(idx, 0)
    assert_equal 2, alpha.size
    assert_equal 2, coef.size
    rhalf = 1.0 / Math.sqrt(2.0)
    assert_in_delta rhalf, coef[0].abs, 1e-12
    assert_in_delta rhalf, coef[1].abs, 1e-12
    assert_in_delta 1.0, coef.sum { |c| c * c }, 1e-12
    assert coef[0] * coef[1] > 0, "singlet: expected the same relative sign"
    assert_equal alpha[0], beta[1]
    assert_equal beta[0], alpha[1]
    refute_equal alpha[0], alpha[1]
  end

  # ms_x2 is a parameter, not an assumption: wrong parity and |Ms| > S
  # are both errors.
  def test_csf2det_rejects_impossible_ms
    drt = ref_drt
    idx = drt.csf_index(HExpr.step_parse("udfe", 4))
    assert_raises(RuntimeError) { drt.csf_ndets(idx, 1) }
    assert_raises(RuntimeError) { drt.csf_ndets(idx, 2) }
  end

  # Multi-word determinants. Every case above runs at norb == 4, where
  # nwords == 1 and the second index of a mask row is never used. This
  # walks a 70-orbital DRT and reconstructs the occupation from the
  # masks, so a wrong word index or stride fails instead of answering
  # differently.
  def test_csf2det_multi_word
    drt = HExpr::Drt.new(kind: "foci", nve: 2, spin_x2: 0,
                         ncore: 1, nval: 1, next: 68)
    assert_equal 70, drt.norb
    nwords = (drt.norb + 63) / 64
    assert_equal 2, nwords

    step_occ = { 0 => 0, 1 => 1, 2 => 1, 3 => 2 }
    saw_high_word = false

    drt.ncsfs.times do |i|
      alpha, beta, coef = drt.csf_to_dets(i, 0)
      assert_equal coef.size, alpha.size
      assert_equal nwords, alpha[0].size
      assert_in_delta 1.0, coef.sum { |c| c * c }, 1e-12

      want = drt.csf_steps(i).map { |s| step_occ.fetch(s) }
      coef.each_index do |k|
        got = (0...drt.norb).map do |p|
          w = p / 64
          b = p % 64
          ((alpha[k][w] >> b) & 1) + ((beta[k][w] >> b) & 1)
        end
        assert_equal want, got, "csf #{i} det #{k}: occupation mismatch"
        saw_high_word = true if alpha[k][1] != 0 || beta[k][1] != 0
      end
    end

    assert saw_high_word,
           "no determinant occupied an orbital >= 64; wide path not taken"
  end
end
