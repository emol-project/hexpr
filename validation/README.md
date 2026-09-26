# validation/

## Overview

This directory contains an independent reference calculation used to validate
the PySCF + hexpr integration pipeline.

The reference system is H2O / 6-31G(d) / SOCI. For this system:

- **OpenMolcas** is the external ground truth -- it was run once, and its
  input file is preserved.
- `pyscf_run.py` is a self-contained PySCF script that validates the RHF
  energy against OpenMolcas and then computes the CI energy using hexpr,
  comparing it against the OpenMolcas RASSCF reference.

## Reference system

| Directory | Basis | Method | RHF target (Hartree) | CI target (Hartree) |
|---|---|---|---|---|
| `h2o_soci/` | 6-31G(d) | SOCI | -76.0105034691 | -76.19806885 |

Nuclear repulsion: **9.1892096196 Hartree**.

## Methodology

`pyscf_run.py` uses PySCF as a computational bridge between OpenMolcas and
hexpr in three steps:

1. **RHF** -- PySCF computes the restricted Hartree-Fock energy and verifies
   it matches OpenMolcas to within the RHF tolerance. This confirms that the
   two codes agree on geometry, basis set, and SCF procedure before any CI
   code is involved.
2. **Integral transformation** -- PySCF transforms AO integrals to the MO
   basis for the chosen active space (one-electron and two-electron
   integrals).
3. **CI via hexpr** -- hexpr builds the DRT and evaluates the Brooks
   expression to produce CI Hamiltonian matrix elements. PySCF diagonalizes
   the CI matrix. The resulting CI energy is compared against the OpenMolcas
   CI reference.

This three-step chain confirms that hexpr correctly evaluates the CI
expression for the given active space and integral data.

## Tolerances

| Quantity | Tolerance |
|---|---|
| Nuclear repulsion energy | 1e-8 Hartree |
| RHF energy | 1e-8 Hartree |
| CI energy | 1e-8 Hartree |

All three quantities agree with OpenMolcas to within 1e-8 Hartree (the
OpenMolcas printed precision is 8 digits for the SOCI total energy).

## Reproducing the validation

Run from the repository root (requires `libhexpr.so`, or `libhexpr.dylib` on
macOS -- see [c/README.md](../c/README.md)):

```sh
python validation/h2o_soci/pyscf_run.py
```

The script prints computed vs. target values and exits 0 (PASS) or 1 (FAIL).

## Notes

- The OpenMolcas input file (`h2o_rasscf.inp`) and geometry (`h2o.xyz`) are
  preserved in `h2o_soci/`. The full OpenMolcas output is **not** committed;
  the reference numerical values above (NRE, RHF, SOCI) suffice to validate
  hexpr. Anyone with OpenMolcas installed can reproduce the reference by
  running `pymolcas h2o_rasscf.inp` in that directory.
- 6-31G(d) note: OpenMolcas with `Basis = 6-31G*` uses Cartesian d functions
  (19 AOs); PySCF defaults to spherical d (18 AOs, energy differs by
  ~1.4e-3 Hartree). `pyscf_run.py` sets `mol.cart = True` to match
  OpenMolcas.

## See also

- Project root: [`../README.md`](../README.md)
- C library: [`../c/README.md`](../c/README.md)
- Public API header: [`../c/include/hexpr.h`](../c/include/hexpr.h)
- Python binding: [`../python/README.md`](../python/README.md)
- Index encoding: [`../docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
- Methodology paper: [DOI placeholder until preprint is up]
