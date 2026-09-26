# h2o_soci -- SOCI(8,17) for water, four engines

This example builds the second-order CI (SOCI) wavefunction for water in a
6-31G* basis with hexpr, and cross-checks the energy across four independent
quantum-chemistry engines. It is a worked example meant to be read; for
validation of hexpr itself see `../../validation/`.

## The calculation

- Molecule / basis: H2O, 6-31G* with **spherical** d functions -- the d shell
  is taken in 5 real spherical-harmonic components (not 6 Cartesian) -> **18 AO**.
- Active space: the O 1s is frozen (1 inactive), leaving **17 active** orbitals
  with 8 electrons. The SOCI DRT partitions them RAS1=2 / RAS2=4 / RAS3=11
  (hexpr `ncore=2, nval=4, next=11`).
- All four engines use the same geometry and active space; only the integral
  source and the basis-function convention differ.

The earlier Cartesian-d form of this example used 19 AO / 18 active / RAS3=14
and gave a SOCI of -76.19806885. The ~2.7e-3 Ha gap to the spherical value
below is a genuine basis-set difference (Cartesian d carries one extra
function), not an error. The spherical numbers are the current reference.

An earlier RAS2=2 reference (still spherical, 17 active) gave -76.19536401;
widening the full-CI reference to RAS2=4 lowers SOCI by ~1.75 mHa to
-76.19711710 and enlarges the CI from 1,431 to 17,685 CSFs.

## Cross-engine reference energies (Hartree)

| Engine            | RHF            | SOCI           |
|-------------------|----------------|----------------|
| Psi4 + hexpr      | -76.0091065449 | -76.1971171057 |
| PySCF + hexpr     | -76.0091065449 | -76.1971171060 |
| Fermi.jl + hexpr  | -76.0091065429 | -76.1971171013 |
| OpenMolcas        | -76.0091065429 | -76.19711710   |

All four agree to ~1e-8 Ha. OpenMolcas and Fermi.jl share an identical SCF;
Psi4 and PySCF share another.

## Files

| File                  | Engine     | Integrals                  | hexpr binding |
|-----------------------|------------|----------------------------|---------------|
| `run_pyscf.py`        | PySCF      | PySCF (`mol.cart=False`)   | Python        |
| `run_psi4.py`         | Psi4       | Psi4 (`puream=True`)       | Python        |
| `run_fermi.jl`        | Fermi.jl   | Fermi.jl / GaussianBasis   | Julia (`HExpr`) |
| `run_hexpr.rb`        | --         | none (expression only)     | Ruby          |
| `molcas_h2o_soci.txt` | OpenMolcas | native RASSCF              | -- (reference deck) |

`run_pyscf.py` and `run_psi4.py` are siblings: only the integral-generation
block differs; the hexpr evaluation, Hamiltonian build, diagonalization, and
CI-vector report are shared. `run_hexpr.rb` only builds the hexpr expression
and prints term counts (no Hamiltonian, no energy). `molcas_h2o_soci.txt` is a
self-contained OpenMolcas deck (native basis blocks + Spherical d card,
coordinates in bohr) plus the reference energies; copy the deck into a file and
run it with pymolcas.

## Compute note

The 4e4o CI has 17,685 CSFs and the hexpr SOCI expression has ~21.4M terms.
The Python/Julia drivers build the dense Hamiltonian via a vectorized
NumPy/array path (the scalar per-term loop is kept commented for clarity) and
diagonalize it densely. The dense 17,685x17,685 matrix needs ~2.5 GB; full
diagonalization with eigenvectors (run_psi4.py's dominant-CSF report) peaks
near ~14 GB RAM, so run it on a machine with >=16 GB. PySCF (eigenvalues only)
and the OpenMolcas Davidson solver are lighter.

## Environment

hexpr is reached through its language bindings, located via environment
variables (the scripts hard-code no path):

    export HEXPR_PREFIX=<your hexpr install prefix>
    export HEXPR_LIBRARY=$HEXPR_PREFIX/lib/libhexpr.so
    export LD_LIBRARY_PATH=$HEXPR_PREFIX/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$HEXPR_PREFIX/share/hexpr/python:$PYTHONPATH   # Python binding
    export JULIA_LOAD_PATH=$HEXPR_PREFIX/share/hexpr/julia:          # Julia binding

Per engine, from a clean shell:

- **PySCF** (`run_pyscf.py`): the hexpr env above + a Python environment with
  numpy and pyscf. Then `python run_pyscf.py`.
- **Psi4** (`run_psi4.py`): the hexpr env above + a conda environment with
  psi4 1.9.1 pinned against libint < 2.10 so libint2.so.2 is present:
  `conda create -n psi4 -c conda-forge "psi4=1.9.1" "libint<2.10" numpy`.
  Then `python run_psi4.py`.
- **Fermi.jl** (`run_fermi.jl`): the hexpr env above (for the Julia HExpr
  binding) **plus an external Fermi.jl checkout and its depot**. This script is
  not covered by `../Project.toml` (which declares only HDF5 + HExpr); point
  Julia at your Fermi.jl project and its local depot, e.g.

      export JULIA_DEPOT_PATH=<Fermi.jl checkout>/local:$HOME/.julia
      julia --project=<Fermi.jl checkout> run_fermi.jl

- **OpenMolcas** (`molcas_h2o_soci.txt`): a built OpenMolcas (`MOLCAS`) with
  `pymolcas` on PATH (pymolcas needs pyparsing). Copy the deck out of the .txt
  and run it.

## See also

- Parent index: `../README.md`
- Validation: `../../validation/`
