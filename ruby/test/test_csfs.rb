require_relative "helper"

class TestCsfs < Minitest::Test
  def setup    = HExpr.init
  def teardown = HExpr.shutdown

  def test_csfs_shape_and_values_s_full_4o4e
    drt = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                          ncore: 0, nval: 4, next: 0)
    rows = drt.csfs

    assert_kind_of Array, rows
    assert_equal drt.ncsfs, rows.size
    assert_equal 20, rows.size

    rows.each do |row|
      assert_kind_of Array, row
      assert_equal drt.norb, row.size
      assert_equal 4, row.size
      row.each do |s|
        assert_kind_of Integer, s
        assert_includes 0..3, s
      end
    end
  end

  def test_csfs_round_trips_into_from_csfset
    drt  = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                           ncore: 0, nval: 4, next: 0)
    rows = drt.csfs

    sub = HExpr::Drt.from_csfset(rows, drt.norb)
    assert_equal drt.norb,  sub.norb
    assert_equal drt.ncsfs, sub.ncsfs
    assert_equal HExpr::HEXPR_DRT_CSFSET, sub.kind

    # The step-code matrix read back out of the sub-DRT matches the input.
    assert_equal rows, sub.csfs
  end

  def test_csfs_on_small_subdrt
    # Two CSFs that share a top node; the sub-DRT exposes a non-empty matrix.
    sub = HExpr::Drt.from_csfset([[0, 0, 3, 3], [0, 1, 2, 3]], 4)
    refute_empty sub.csfs
    sub.csfs.each { |row| assert_equal 4, row.size }
  end
end
