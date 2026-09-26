"""Shared fixtures and case definitions for hexpr pytest suite."""
from pathlib import Path
import hexpr

# Ensure the library is initialized before any test runs.
hexpr.init()

REPO_ROOT = Path(__file__).parent.parent.parent  # .../hexpr/
REF_DIR   = REPO_ROOT / "reference"

# Each entry: (case_name, kwargs for hexpr.Drt constructor)
# Mirrors run_reference.rb's CASES list.
CASES = [
    ("s_full_4o4e",  dict(kind="soci", nve=4,  spin_x2=0, ncore=0, nval=4, next_=0)),
    ("t_full_4o4e",  dict(kind="soci", nve=4,  spin_x2=2, ncore=0, nval=4, next_=0)),
    ("s_full_6o6e",  dict(kind="soci", nve=6,  spin_x2=0, ncore=0, nval=6, next_=0)),
    ("t_full_6o6e",  dict(kind="soci", nve=6,  spin_x2=2, ncore=0, nval=6, next_=0)),
    ("s_soci_small", dict(kind="soci", nve=6,  spin_x2=0, ncore=1, nval=4, next_=1)),
    ("s_foci_small", dict(kind="foci", nve=6,  spin_x2=0, ncore=1, nval=4, next_=1)),
    ("t_soci_small", dict(kind="soci", nve=6,  spin_x2=2, ncore=1, nval=4, next_=1)),
    ("t_foci_small", dict(kind="foci", nve=6,  spin_x2=2, ncore=1, nval=4, next_=1)),
    ("s_soci_med",   dict(kind="soci", nve=10, spin_x2=0, ncore=2, nval=4, next_=4)),
    ("s_foci_med",   dict(kind="foci", nve=10, spin_x2=0, ncore=2, nval=4, next_=4)),
    ("t_soci_med",   dict(kind="soci", nve=10, spin_x2=2, ncore=2, nval=4, next_=4)),
    ("t_foci_med",   dict(kind="foci", nve=10, spin_x2=2, ncore=2, nval=4, next_=4)),
]
