# example_ci/

This directory contains worked examples that drive **hexpr** end-to-end
to produce CI energies. They are written to be read, not to validate
hexpr -- for validation, see [`../validation/`](../validation/).

Each example is a self-contained quantum-chemistry mini-program. The
default geometry and basis are hard-coded for a small, fast case, but
the scripts are written so a curious user can swap in a different
geometry or basis and re-run.

## Examples

| Directory | Method | System | API focus |
|---|---|---|---|
| [`h2_3csfs/`](h2_3csfs/) | full CI from 3 explicit CSFs | H2 / STO-3G / CAS(2,2) | `csf_subspace_eval` (high-level CSF-set API) |
| [`h2o_soci/`](h2o_soci/) | SOCI(8,17) full CI | H2O / 6-31G(d), spherical | `Drt('soci', ...)` + `Expression.terms_full` |

Each directory contains, for the same example:

- Python entry point(s):
  - `h2_3csfs/` has three: `run.py` (hard-coded integrals, no PySCF
    needed), `run_pyscf.py` (PySCF-based, full CI with arbitrary
    geometry/basis), and `run_psi4.py` (the same calculation with the
    integrals computed by Psi4 instead of PySCF).
  - `h2o_soci/` has three: `run_pyscf.py` (PySCF-based, the main
    SOCI driver), `run_psi4.py` (the same SOCI calculation with the
    integrals computed by Psi4), and `run_pyscf_foci.py` (a PySCF-based
    FOCI(8,17) variant of the same H2O system).
- Julia versions of the same calculation: `h2_3csfs/` has two -- `run.jl`
  (hard-coded integrals, optional HDF5 load; dependency-free) and
  `run_fermi.jl` (Fermi.jl computes its own integrals, via the
  GaussianBasis.jl library); `h2o_soci/` has `run_fermi.jl` (Fermi.jl
  computes its own integrals -- no `integrals.h5`; see `h2o_soci/README.md`
  for its Fermi.jl setup).
- `run_hexpr.rb` -- Ruby version that drives hexpr through the expression
  step only (no Hamiltonian assembly or diagonalization)
- `molcas_*.txt` -- the OpenMolcas input and reference energies if you
  want to spot-check the output against an independent QC code

Read all three language versions side by side for a complete picture of
how hexpr is used from each binding.

## Running the Python examples

```sh
python example_ci/h2_3csfs/run.py           # hard-coded integrals, no PySCF
python example_ci/h2_3csfs/run_pyscf.py     # PySCF, arbitrary geometry/basis
python example_ci/h2_3csfs/run_psi4.py      # Psi4 instead of PySCF
python example_ci/h2o_soci/run_pyscf.py     # PySCF, the main H2O driver
python example_ci/h2o_soci/run_psi4.py      # Psi4 instead of PySCF
```

The H2 `run.py` is self-contained (integrals hard-coded). The H2
`run_pyscf.py` and the H2O `run_pyscf.py` call PySCF to compute the
MO integrals. (The H2O Julia engine, run_fermi.jl, computes its own
integrals through Fermi.jl, so no integrals.h5 is exchanged for H2O.)

## Running with Psi4 instead of PySCF

Both cases also ship a `run_psi4.py` that computes the same active-space
integrals with **Psi4** rather than PySCF. Only the integral-generation
block differs; the hexpr evaluation, the Hamiltonian build, the
diagonalization, and the CI-vector report are shared with `run_pyscf.py`.
Running the Psi4 variant is an independent check that the hexpr evaluation
is not tied to one integral backend. Both backends reproduce the
OpenMolcas reference energies to ~1e-8 Ha:

```
H2O SOCI : -76.1971171057
H2 roots :  -1.1372701747   -0.1699013905   +0.4798361182
```

```sh
python example_ci/h2_3csfs/run_psi4.py
python example_ci/h2o_soci/run_psi4.py
```

`run_psi4.py` requires `psi4` (PySCF is not needed). Note that
`psi4 1.9.1` links `libint2.so.2`; pin libint **below 2.10** so that
SONAME is present, otherwise Psi4 fails to import:

```sh
conda create -n psi4 "psi4=1.9.1" "libint<2.10" numpy   # -> libint 2.9.0
```

Two Psi4 settings make the integrals line up with OpenMolcas and the
PySCF script: 6-31G* with **spherical d** functions (`puream True`, 18
AOs for water -- the H2O script asserts `nbf == 18`), and `scf_type pk`
(exact ERIs, no density fitting). Psi4 has no one-call active-space
effective Hamiltonian, so `run_psi4.py` builds `ecore` / `h1eff` by hand
following PySCF's `h1e_for_cas` (for the H2 case, which has no frozen
core, this reduces to `ecore = ENUC` and a bare `h1eff`). For H2O,
neither run_pyscf.py nor run_psi4.py writes integrals.h5; the Julia
engine run_fermi.jl computes its own integrals via Fermi.jl.

## Running the Julia examples

There are two Julia engines, with different environment needs.

**`run.jl`** (`h2_3csfs/` only) is dependency-free: hard-coded integrals,
no PySCF / Psi4 / Fermi.jl. It runs under the separate Julia environment
defined by `example_ci/Project.toml` (kept separate from `julia/Project.toml`
so that the HExpr binding itself has no HDF5 dependency).

First time only:

```sh
julia --project=example_ci -e 'using Pkg; Pkg.instantiate()'
```

Then:

```sh
julia --project=example_ci example_ci/h2_3csfs/run.jl
```

**`run_fermi.jl`** (`h2_3csfs/` and `h2o_soci/`) uses Fermi.jl to compute its
own RHF and integrals (spherical only). It needs the hexpr env (for the Julia
HExpr binding) **plus an external Fermi.jl checkout and its depot**. This
script is **not** covered by `example_ci/Project.toml` (which declares only
HDF5 + HExpr); point Julia at your Fermi.jl project and its local depot, e.g.

```sh
export JULIA_DEPOT_PATH=<Fermi.jl checkout>/local:$HOME/.julia
julia --project=<Fermi.jl checkout> example_ci/h2_3csfs/run_fermi.jl
julia --project=<Fermi.jl checkout> example_ci/h2o_soci/run_fermi.jl
```

See `example_ci/h2o_soci/README.md` for the full Fermi.jl checkout, depot, and
env details.

## Running the Ruby examples

```sh
cd ruby
bundle exec ruby -Ilib ../example_ci/h2_3csfs/run_hexpr.rb
bundle exec ruby -Ilib ../example_ci/h2o_soci/run_hexpr.rb
```

The Ruby examples only build the hexpr expression and print the term counts.

## Customizing the calculations

For `h2_3csfs/`, use `run.py` for the hard-coded default case (no
PySCF required), or `run_pyscf.py` to recompute the integrals from
scratch via PySCF at any geometry or basis. For `h2o_soci/`,
`run_pyscf.py` is the only PySCF driver; edit the molecule
definition at the top to study a different geometry or basis. The
`run_psi4.py` scripts mirror their `run_pyscf.py` siblings -- edit the
molecule block at the top the same way (for the H2O case, keep
`puream True` so the spherical-d basis still matches OpenMolcas).

`run.jl` has a hard-coded integrals block at the top for the default
case (so the example runs out of the box). For `h2_3csfs/`, `run.jl`
also has a commented-out HDF5-load alternative -- comment out the
hard-coded block and uncomment the alternative block to switch modes.

## See also

- Project root: [`../README.md`](../README.md)
- Validation: [`../validation/`](../validation/)
- Python binding: [`../python/README.md`](../python/README.md)
- Julia binding: [`../julia/README.md`](../julia/README.md)
- Ruby binding: [`../ruby/README.md`](../ruby/README.md)
- Index encoding: [`../docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
