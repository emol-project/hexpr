# CH2 spin-orbit coupling: PySCF + hexpr vs determinants and GAMESS HSO1

CH2 / cc-pVDZ (spherical), RHF orbitals, CAS(6,6) with the C 1s frozen.
The one-electron Breit-Pauli spin-orbit operator is evaluated between the
singlet and triplet CAS states with `hexpr.rank1_pair` (reduced matrix
elements of (a+_p a~_q)^1) and one 3j-symbol per magnetic component.

## Files

| file | what |
|---|---|
| `run_pyscf_so.py` | PySCF integrals, spin-free CAS from hexpr, spin-orbit matrix from hexpr; determinant check; optional `--gamess LOG`, `--zeff C=..,H=..` |
| `gamess/ch2_rhf.inp` | GAMESS RHF (same Bohr geometry and basis) |
| `gamess/make_so_inp.sh` | writes the GAMESS SO-CI input (`RUNTYP=TRANSITN`, `OPERAT=HSO1`) with the RHF `$VEC` |
| `gamess/*.log` | GAMESS 2024 R2 Patch 1 output (not distributed; written by the commands below) |
| `results/*.txt` | output of `run_pyscf_so.py` for the two comparisons |

## How to run

```
cd gamess
~/gamess/rungms ch2_rhf 00 1 > ch2_rhf.log
bash make_so_inp.sh TRUE   && ~/gamess/rungms ch2_so_true  00 1 > ch2_so_true.log
bash make_so_inp.sh 3-21G  && ~/gamess/rungms ch2_so_3_21g 00 1 > ch2_so_3_21g.log
cd ..
python run_pyscf_so.py --gamess gamess/ch2_so_true.log
python run_pyscf_so.py --zeff C=3.9,H=1.0 --gamess gamess/ch2_so_3_21g.log
```

## Results (2026-10-16)

1. **Determinant check.** Every singlet-triplet and triplet-triplet CSF pair
   (175 and 189 CSFs), every magnetic component: hexpr and the determinant
   construction agree to 1.5e-19 Eh (elements up to 9e-5 Eh).
2. **Spin-free part.** RHF (-38.8670698382 Eh) and the CAS energies agree
   with GAMESS and with PySCF CASCI to all printed digits.
3. **GAMESS HSO1.** All 444 printed matrix elements between the 16 GAMESS
   states agree in magnitude to the print rounding (5e-5 cm-1), for true
   nuclear charges (`ZEFTYP=TRUE`) and for GAMESS's effective charges
   (`ZEFTYP=3-21G`, C 3.9, H 1.0), once GAMESS's hartree-to-cm-1 factor is
   taken into account: it is 1.00001703 times the CODATA 2018 value used
   here, measured independently from GAMESS's own list of spin-orbit levels
   (spread 2e-8). The same factor appears in GAMESS's spin-free energy
   differences, so it is a unit, not a property of the spin-orbit elements.

| S0-T1 coupling, sqrt(sum_M \|<S0\|H_SO\|T1,M>\|^2) | hexpr | GAMESS |
|---|---|---|
| true charges | 19.8487 cm-1 | 19.849 cm-1 |
| Z_eff (C 3.9, H 1.0) | 13.0244 cm-1 | 13.024 cm-1 |

Notes on the GAMESS side:

- The C1 Davidson solver converges only to states whose irreducible
  representation (of C2v) is present in its start vectors, so some singlet
  roots (exact roots 2, 4, 8, ...) never appear. The script matches GAMESS
  states to the exact roots by energy.
- `IROOTS = NSTATE`: with `NSTATE` larger than `IROOTS` this GAMESS version
  stopped in the spin-orbit step ("I/O past end of record").
- `SYMTOL=1e-10`: the default 1e-4 leaves up to 2e-3 cm-1 error in the
  printed elements.
- Integral convention: PySCF's `int1e_pnucxp` is (p* | nuc x p |), so the
  Hermitian operator is h_pq = (alpha^2/2) * 1j * `int1e_pnucxp`; for
  effective charges `int1e_pnucxp` is replaced by
  -sum_A Z_eff(A) `int1e_prinvxp`(A).
