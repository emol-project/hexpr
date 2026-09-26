"""Tests for the binding-completeness wrappers.

Covers csf_steps, step_parse, and the *_show* methods added so that every
public C function is callable from the Python binding.  Reference case:
s_full_4o4e (norb=4, ncsfs=20).
"""
import io

import pytest
import hexpr


def _ref_drt():
    return hexpr.Drt(kind="soci", nve=4, spin_x2=0, ncore=0, nval=4, next_=0)


def test_step_parse_fude():
    assert list(hexpr.step_parse("fude", 4)) == [3, 1, 2, 0]


def test_csf_steps_round_trip():
    """csf_steps matches the csfs matrix and round-trips through csf_index."""
    drt = _ref_drt()
    csfs = drt.csfs()
    assert drt.ncsfs == 20
    for i in range(drt.ncsfs):
        steps = drt.csf_steps(i)
        assert list(steps) == list(csfs[i])
        assert drt.csf_index(steps) == i


def test_node_info_known_nodes():
    drt = _ref_drt()
    assert drt.node_info(0) == (0, 0, 0)
    assert drt.node_info(drt.nnodes - 1) == (4, 4, 0)


def test_show_methods_stdout_no_raise():
    """show* with default stdout must be callable and not raise."""
    drt = _ref_drt()
    expr = hexpr.Expression(drt)
    drt.show()
    expr.show_full()
    expr.show_one()


def test_show_non_stdout_not_implemented():
    drt = _ref_drt()
    with pytest.raises(NotImplementedError):
        drt.show(file=io.StringIO())
