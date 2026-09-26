require_relative "helper"

class TestExpr < Minitest::Test
  def setup    = HExpr.init
  def teardown = HExpr.shutdown

  def test_s_full_4o4e
    drt  = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0,
                           ncore: 0, nval: 4, next: 0)
    expr = HExpr::Expression.new(drt)

    assert_equal 638, expr.nterms_full
    assert_equal 124, expr.nterms_one

    ij, pqrs, coef = expr.read_full
    assert_equal 638, ij.length
    assert_equal 638, pqrs.length
    assert_equal 638, coef.length
    assert_equal 0,   ij.first
    assert_equal 5,   pqrs.first
    assert_in_delta 2.0, coef.first, 1e-12

    ij1, pqrs1, coef1 = expr.read_one
    assert_equal 124, ij1.length
    assert_equal 124, pqrs1.length
    assert_equal 124, coef1.length
  end
end
