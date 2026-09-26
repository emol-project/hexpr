require_relative "helper"

class TestCsfsetDrt < Minitest::Test
  def setup    = HExpr.init
  def teardown = HExpr.shutdown

  # CSF step codes for s_full_4o4e (4 electrons in 4 orbitals, S=0).
  # From reference/s_full_4o4e.drt.txt "CSF List" section, e=0 u=1 d=2 f=3.
  CSFS_4O4E = [
    [0, 0, 3, 3], [0, 1, 2, 3], [0, 1, 3, 2], [0, 3, 0, 3],
    [0, 3, 1, 2], [0, 3, 3, 0], [1, 0, 2, 3], [1, 0, 3, 2],
    [1, 1, 2, 2], [1, 2, 0, 3], [1, 2, 1, 2], [1, 2, 3, 0],
    [1, 3, 0, 2], [1, 3, 2, 0], [3, 0, 0, 3], [3, 0, 1, 2],
    [3, 0, 3, 0], [3, 1, 0, 2], [3, 1, 2, 0], [3, 3, 0, 0],
  ].freeze

  def test_round_trip_full_4o4e
    full = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                           ncore: 0, nval: 4, next: 0)
    assert_equal 4,  full.norb
    assert_equal 20, full.ncsfs
    assert_equal 14, full.nnodes

    sub = HExpr::Drt.from_csfset(CSFS_4O4E, 4)
    assert_equal 4,  sub.norb
    assert_equal 20, sub.ncsfs
    assert_equal 14, sub.nnodes
    assert_equal HExpr::HEXPR_DRT_CSFSET, sub.kind
  end

  def test_two_csf_subset
    csfs = [[0, 0, 3, 3], [0, 1, 2, 3]]
    sub = HExpr::Drt.from_csfset(csfs, 4)
    assert_equal 4, sub.norb
    assert sub.ncsfs >= 2
    assert sub.nnodes < 14
  end

  def test_input_order_and_walk_superset
    # Non-canonical input triggering a canonical superset: N=2, M=4.
    inp = [[1, 0, 3, 2],   # uefd
           [0, 1, 2, 3]]   # eudf
    sub = HExpr::Drt.from_csfset(inp, 4)

    # csfs/ncsfs == the input rows in input order
    assert_equal 2, sub.ncsfs
    assert_equal inp, sub.csfs

    # walk_* == canonical superset M >= N
    assert_equal 4, sub.walk_ncsfs
    assert sub.walk_ncsfs >= sub.ncsfs
    walk = sub.walk_csfs
    assert_equal 4, walk.length
    inp.each { |row| assert_includes walk, row }
  end

  def test_walk_equals_csfs_for_soci
    soci = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                          ncore: 0, nval: 4, next: 0)
    assert_equal soci.ncsfs, soci.walk_ncsfs
    assert_equal soci.csfs,  soci.walk_csfs
  end

  def test_invalid_step_raises
    csfs = [[0, 1, 2, 4]]  # step=4 is invalid; C library returns NULL
    assert_raises(RuntimeError) {
      HExpr::Drt.from_csfset(csfs, 4)
    }
  end

  def test_shape_mismatch_raises
    csfs = [[0, 0, 3]]  # row length 3, norb=4
    assert_raises(ArgumentError) {
      HExpr::Drt.from_csfset(csfs, 4)
    }
  end

  def test_non_array_raises
    assert_raises(ArgumentError) {
      HExpr::Drt.from_csfset("not an array", 4)
    }
  end
end
