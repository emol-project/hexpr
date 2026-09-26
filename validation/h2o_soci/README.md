# H2O / 6-31G(d) / SOCI

## System

| Property | Value |
|---|---|
| Molecule | H2O |
| Charge / Multiplicity | 0 / 1 (singlet) |
| Symmetry | C1 (none) |
| Basis | 6-31G(d) = 6-31G* (Cartesian d functions, 19 AOs) |
| Method | SOCI (Second-Order CI) |
| Total orbitals | 19 |
| Total electrons | 10 |

## Geometry (Angstrom)
O   0.0000000   0.0000000   0.0000000
H   0.0000000   0.7572157   0.5865358
H   0.0000000  -0.7572157   0.5865358

## Reference values (from OpenMolcas)

| Quantity | Value (Hartree) |
|---|---|
| Nuclear repulsion energy | 9.1892096196 |
| RHF energy | -76.0105034691 |
| SOCI CI energy (root 1) | -76.19806885 |

> **Note.** These numbers use the earlier Cartesian-d form of the basis
> (19 AO / 18 active / RAS3=14). The current reference is the spherical-d
> form (18 AO / 17 active / RAS3=11), giving SOCI = -76.19711710; see
> `../../example_ci/h2o_soci/README.md`. The ~2.7e-3 Ha gap is a genuine
> basis-set difference, not a discrepancy.

## OpenMolcas partition

The RASSCF block in `h2o_rasscf.inp` declares:

| Region | MOs | RASSCF keyword | Description |
|---|---|---|---|
| Inactive | #1 | `Inactive = 1` | O 1s, doubly occupied, not in CI |
| RAS1 (hole space) | #2-3 | `RAS1 = 2`, max 2 holes | hexpr core |
| RAS2 (fully active) | #4-5 | `RAS2 = 2` | hexpr active (valence) |
| RAS3 (particle space) | #6-19 | `RAS3 = 14`, max 2 particles | hexpr external |
| Active electrons | -- | `nActEl = 8 2 2` | 8 active electrons |

`CIONLY` disables orbital optimization, so the calculation is a single CI
diagonalization with the SCF orbitals.

## hexpr partition

```python
hexpr.Drt('soci', nve=8, spin_x2=0, ncore=2, nval=2, next_=14)
```

The `pqrs` integers returned by hexpr index into the 18-orbital active space
(everything except the single inactive O 1s, MO #1).

## Tolerances

| Quantity | Tolerance |
|---|---|
| Nuclear repulsion energy | 1e-8 Hartree |
| RHF energy | 1e-8 Hartree |
| CI energy | 1e-8 Hartree |

## How to run

```sh
python validation/h2o_soci/pyscf_run.py
```

Run from the repository root. The script is self-contained: no arguments
needed. It computes RHF and SOCI energies via PySCF + hexpr, prints
computed vs. target values for NRE, RHF, and SOCI, and exits 0 (PASS) or
1 (FAIL).

## OpenMolcas reference files

- `h2o_rasscf.inp` -- OpenMolcas input (preserved as evidence).
- `h2o.xyz` -- Geometry referenced by the input.

The full OpenMolcas output is not committed; the reference numerical values
above (NRE, RHF, SOCI) suffice to validate hexpr. Reproduce by running
`pymolcas h2o_rasscf.inp` in this directory.
