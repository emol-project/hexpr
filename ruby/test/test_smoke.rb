require_relative "helper"

class TestSmoke < Minitest::Test
  def setup    = HExpr.init
  def teardown = HExpr.shutdown

  def test_version_string_nonempty
    refute_empty HExpr.version_string
  end

  def test_version_components
    assert_equal 1, HExpr.version_major
    assert_equal 0, HExpr.version_minor
    assert_equal 0, HExpr.version_patch
  end
end
