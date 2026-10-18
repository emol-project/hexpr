# OH X 2Pi spin-orbit splitting: PySCF + hexpr vs GAMESS HSO1 and experiment

OH at R = 0.9697 A (given in Bohr), cc-pVTZ (spherical), O 1s core,
CAS(7,5) on O 2s 2p and H 1s. Orbitals from a state-averaged CASSCF over the
two components of X 2Pi, so that they are equivalent and the components stay
exactly degenerate. The one-electron Breit-Pauli spin-orbit operator is built
between the doublet CSFs with `hexpr.rank1_pair` (shared code in
`../so_tools.py`), the state-interaction matrix over the four lowest doublet
roots is diagonalized, each level is named by Omega = |J_z| (J_z = L_z + S_z,
diagonalized inside each Kramers pair), and

    A = E(2Pi_3/2) - E(2Pi_1/2).

## Files

| file | what |
|---|---|
| `run_pyscf_oh.py` | SA-CASSCF (PySCF), CAS and spin-orbit matrix (hexpr), determinant check, levels with Omega, A; `--zeff`, `--gamess` |
| `run_all.sh` | the two hexpr runs below (after `gamess/run_all.sh`) |
| `gamess/make_inp.sh`, `gamess/run_all.sh` | GAMESS ROHF -> SA-MCSCF (CISTEP=GUGA) -> SO-CI (RUNTYP=TRANSITN, OPERAT=HSO1) |
| `results/*.txt` | outputs of 2026-10-16 |
| `gamess/*.log` | GAMESS 2024 R2 Patch 1 outputs (not distributed; written by `gamess/run_all.sh`) |

## Results (2026-10-16)

| | hexpr | GAMESS HSO1 | difference |
|---|---|---|---|
| ROHF energy (Eh) | -75.4144656124 | -75.4144656124 | 0 |
| SA-CASSCF energy, both Pi components (Eh) | -75.4376161432 | -75.4376161432 | 0 |
| A, true nuclear charges (cm-1) | -213.7549 | -213.7544 | 5e-4 |
| A, Z_eff O 6.0, H 1.0 (`ZEFTYP=3-21G`) (cm-1) | -160.6357 | -160.6353 | 4e-4 |

- The determinant check (every doublet CSF pair, both M_S) agrees with hexpr
  to 6e-19 Eh.
- Both programs put 2Pi_3/2 below 2Pi_1/2 (inverted), as observed.
- GAMESS values are divided by its hartree-to-cm-1 factor, 1.00001703 times
  CODATA 2018 (measured in `../ch2_so/README.md`); it prints the levels to
  three decimals, so the differences are within its rounding.
- GAMESS's effective charges for `ZEFTYP=3-21G` are printed as O 6.0, H 1.0;
  hexpr with exactly 6.0 reproduces its A, so the value is 6.0.
- The higher roots (the Sigma states, not averaged in the MCSCF) differ by
  1.8e-5 Eh between the two programs: GAMESS ends its MCSCF with "the
  Lagrangian matrix is not entirely symmetric", i.e. slightly different
  orbitals in directions to which the averaged energy is stationary but the
  other states are not. This does not reach A beyond its printed digits.

Experiment: A = -139.2 cm-1 (X 2Pi, v = 0). The one-electron operator with
true charges overestimates it by 54 %, as expected without the two-electron
(shielding) terms; the effective charges recover most of that (-160.6 cm-1,
15 % high) in this small CAS and basis.

### Two-electron terms (`--somf`)

hexpr takes one-electron operators, so the two-electron Breit-Pauli terms
enter as a mean field built from PySCF's `int2e_p1vxp1`
(`so_tools.somf_2e`): sum_kl D_kl [g(pq|kl) - 3/2 g(pk|lq) - 3/2 g(lq|pk)]
with g = (alpha^2/2) 1j `int2e_p1vxp1` (Hess et al.).

| operator | hexpr | GAMESS | difference |
|---|---|---|---|
| 1e + mean field of the O 1s core (exact core-active part) | -144.9735 | HSO2P -144.9735 | 7e-5 |
| 1e + mean field of the SA-CASSCF density (SOMF) | -135.6014 | HSO2 (exact 2e) -136.7607 | 1.16 |
| experiment | | | -139.2 |

- With the closed-shell core density the mean field is exact for the
  core-active two-electron terms, which is what GAMESS `OPERAT=HSO2P`
  computes; the two agree to its printed digits. This fixes the sign of
  `int2e_p1vxp1` and the -3/2 exchange factors.
- With the full density the mean field approximates the active-active
  two-electron terms that GAMESS `OPERAT=HSO2` treats exactly; the
  difference, 1.2 cm-1 (0.9 %), is the mean-field error of SOMF here.

## Run

```
bash gamess/run_all.sh          # GAMESS 2024 R2 in ~/gamess
PYTHON=python bash run_all.sh   # one-electron operator, true and effective charges
PYTHON=python bash run_somf.sh  # two-electron mean field vs HSO2P / HSO2
```
