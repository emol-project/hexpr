require_relative "helper"

class TestDrt < Minitest::Test
  def setup    = HExpr.init
  def teardown = HExpr.shutdown

  def test_s_full_4o4e
    drt = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                          ncore: 0, nval: 4, next: 0)
    assert_equal 4,  drt.norb
    assert_equal 20, drt.ncsfs
    assert_equal 14, drt.nnodes
    assert_equal HExpr::HEXPR_DRT_SOCI, drt.kind
    assert_equal 4,  drt.nve_state
    assert_equal 0,  drt.spin_x2
  end
end
