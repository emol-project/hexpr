# Rank-1 one-electron spin operators (`hexpr_rank1_pair`, since 1.2)

This note describes what `hexpr_rank1_pair` returns, the conventions behind
it, how it is evaluated with the per-orbital factors of the spin-free path,
and how it was checked. It is the reference for the API comment in
`c/include/hexpr.h`; the spin-free formulation it builds on is that of the
hexpr paper (preprint, doi:10.5281/zenodo.22977329).

## 1. Operators

With the spin-adapted annihilator (Sasaki's convention, as for the spin-free
operators of hexpr)

    a~_{q,+1/2} = a_{q,beta},    a~_{q,-1/2} = -a_{q,alpha},

the one-electron triplet operators `T^1(pq) = (a+_p a~_q)^1` have the
spherical components

    T^1_{+1}(pq) =  a+_{p alpha} a_{q beta}
    T^1_{0}(pq)  = -(a+_{p alpha} a_{q alpha} - a+_{p beta} a_{q beta}) / sqrt(2)
    T^1_{-1}(pq) = -a+_{p beta} a_{q alpha}

The one-electron spin tensor
`s_m(pq) = sum_{sigma tau} (s_m)_{sigma tau} a+_{p sigma} a_{q tau}` is then

    s_m(pq) = -T^1_m(pq) / sqrt(2)          (m = +1, 0, -1),

the same form as the number operator `n_p = -sqrt(2) [a+_p a~_p]^0` of the
spin-free path. Summed over orbitals, `sum_p T^1(pp) = -sqrt(2) S`.

## 2. What is returned: reduced matrix elements (Edmonds)

For a bra CSF `Psi` of total spin `S'` and a ket CSF `Phi` of total spin `S`,

    <S'M'| T^1_m(pq) |SM> = (-1)^(S'-M') ( S' 1 S ; -M' m M ) R_pq(Psi, Phi),

with the Wigner 3j symbol `( S' 1 S ; -M' m M )`. `hexpr_rank1_pair` returns
the list of non-zero `(p, q, R_pq)` (0-based orbitals, `|R| > 1e-12`). `R`
contains neither integrals nor `M`: it depends only on the two step vectors,
so it can be stored once per CSF pair and combined with any operator and any
magnetic component afterwards.

The two CSFs may have different `S`; no DRT is built (a DRT has one head,
i.e. one total spin). Step codes are those of `hexpr_csf_steps`:
`0..3 = e, u, d, f`, `step[0]` the lowest orbital.

### Spin-orbit operator

For a one-electron spin-orbit operator with spatial vector integrals
`h_pq = (h^x, h^y, h^z)_pq`, spherical components
`h_{pq,+-1} = -+(h^x +- i h^y)/sqrt(2)`, `h_{pq,0} = h^z`,

    H_SO = sum_pq sum_m (-1)^m h_{pq,-m} s_m(pq)

and

    <S'M'| H_SO |SM> = -(1/sqrt 2) (-1)^(S'-M')
                       sum_m (-1)^m ( S' 1 S ; -M' m M ) sum_pq h_{pq,-m} R_pq.

`example_ci/so_tools.py` (`so_blocks`) does exactly this with PySCF
Breit-Pauli integrals; `example_ci/ch2_so` and `example_ci/oh_so` use it.

### Identities

- Total spin: `sum_p R_pp(I, I) = -sqrt(2) sqrt(S(S+1)(2S+1))` for every CSF
  `I` with `S > 0`.
- Adjoint: `R_qp(Phi, Psi) = (-1)^(S'-S) R_pq(Psi, Phi)`. Both `p > q` and
  `p < q` are returned, so the user never has to apply it.

## 3. Evaluation

### 3.1 The operator chain

Orbitals are taken from the bottom up, as in the spin-free path. The chain
carries an intermediate rank; `theta` below is the multiplicity `2s+1` of
the coupled operator after that orbital.

| case | operator | lower orbital | upper orbital | sign `sigma` |
|---|---|---|---|---|
| `p > q` | `(a+_p a~_q)^1` | `q`: `a~`, theta 2 | `p`: `a+`, theta 3 | +1 |
| `p < q` | `-(a~_q a+_p)^1` | `p`: `a+`, theta 2 | `q`: `a~`, theta 3 | -1 |
| `p = q` | `[a+_p a~_p]^1` | `p`: `[a+ a~]^1`, theta 3 | -- | +1 |

The sign for `p < q`: exchanging two spin-1/2 factors gives
`(AB)^k = (-1)^(1-k) (BA)^k`, and the two fermion operators on different
orbitals anticommute; for `k = 0` this is the spin-free +1, for `k = 1` it is
-1.

### 3.2 Reduced matrix element

    R_pq(Psi, Phi) = (-1)^(1 + S' - S) sqrt(2S + 1) sigma  prod_{j=1}^{N_orb} L_j

where `L_j` are the per-orbital factors of the spin-free path, unchanged
(the general recoupling coefficient of hexpr, which is already written for
any operator rank). On orbitals without an operator the factor is the unit
operator and `theta` is carried over. Three things differ from the spin-free
(rank 0) evaluation:

1. the rank of the chain after the last operator is 1 (theta 3), not 0;
2. the product runs to the top orbital: the walks do not close into a loop;
3. the final normalisation is that of a reduced matrix element,
   `(-1)^(1+S'-S) sqrt(2S+1)`, instead of the factor that turns a scalar
   into an ordinary matrix element.

The same formula with rank 0 gives `sqrt(2S+1)` times the spin-free value,
i.e. the Edmonds reduced matrix element of a scalar; this was used during
development as a consistency check of the convention.

### 3.3 Which orbitals contribute

- **Below the lowest operator** the rank is 0: `L_j = 1` where the two walks
  agree, and a difference there makes the element zero. Evaluation starts at
  `min(p, q)`.
- **Between the two operators**: as in the spin-free path (theta 2).
- **Above the last operator** (the tail) the rank stays 1 up to the top.
  Unlike the spin-free case, `L_j` is not 1 on singly occupied orbitals even
  when the walks agree, and a crossing (`u` against `d`) can be non-zero; this
  is the triplet-operator loop of GUGA that reaches the head. But an orbital
  that is empty in both CSFs or doubly occupied in both has `L_j = 1` for any
  partial spins, so those orbitals are skipped and evaluation over occupied
  orbitals only (as used by the symbolic-matrix method of hexpr2) remains
  valid in the tail.

### 3.4 Selection rules and size

- Same number of electrons; occupations equal, or differing by one electron
  moved from `q` (in the ket) to `p` (in the bra).
- Equal occupations leave only `p = q` on orbitals that are open shells in
  both CSFs.
- `|S' - S| <= 1 <= S' + S` (singlet-singlet elements vanish).
- Hence at most one term for different occupations and at most the number
  of open shells for equal ones; `hexpr_rank1_max_terms(norb) = norb` bounds
  both.
- Cost: `O(N_orb)` per `(p, q)`, no DRT.
- Limits: `norb <= 256`, partial `2S <= 14` along both walks
  (`HEXPR_ERR_INVALID_ARG` / `HEXPR_ERR_RANGE` otherwise). Thread-safe.

## 4. Worked examples

Steps are written bottom orbital first (`"udf"` = orbital 1 `u`, orbital 2
`d`, orbital 3 `f`); orbitals in `R` are 1-based here, 0-based in the API.
These values are asserted by the tests of every binding.

| bra (S') | ket (S) | element | value |
|---|---|---|---|
| `udf` (0) | `fuu` (1) | `R_31` | `+sqrt(3/2)` |
| `uuu` (3/2) | `udu` (1/2) | `R_11` | `+sqrt(2)` |
| `ud` (0) | `uu` (1) | `R_11` | `-sqrt(3/2)` |
| `ud` (0) | `uu` (1) | `R_22` | `+sqrt(3/2)` |

For `ud`/`uu` the singlet-triplet spin density is `R_11 = -R_22`, as it
must be, since the total spin operator `sum_p T^1(pp)` has no
singlet-triplet element.

Python:

```python
import hexpr
p, q, R = hexpr.rank1_pair([1, 2], [1, 1])     # bra "ud", ket "uu"
# p = [0, 1], q = [0, 1], R = [-1.2247..., 1.2247...]
```

## 5. How it was checked

- During development, the formulation was compared with determinants built
  by `hexpr_csf_to_dets`: every `M`, every `m`, every `(p, q)` and every CSF
  pair up to 4 orbitals, all electron counts, `2S <= 4`. Results: `R`
  independent of `M` (Wigner-Eckart) to 9e-16; the rank-0 and rank-1
  normalisations as in Sec. 3.2; the adjoint sign; the total-spin identity;
  tail factors of empty and doubly occupied orbitals equal to 1 for
  `2S <= 7`. Chains of up to four operators with final rank 0, 1 and 2 on
  4 orbitals follow `(-1)^(k+S'-S) sqrt(2S+1) prod_j L_j`.
- The C evaluator against determinants: every pair and `(p, q)` up to 4
  orbitals and random pairs at 6 and 8 orbitals, largest difference 4.4e-15.
- `c/tests/test_rank1.c` (run by `make check`): the total-spin identity and
  the adjoint relation over every CSF pair up to 5 orbitals, the worked
  examples, the selection rules and argument checks. Python, Fortran, Julia
  and Ruby have the same checks in their test suites.
- `example_ci/ch2_so`, `example_ci/oh_so`: spin-orbit matrices from PySCF
  integrals agree with the determinant construction to 1e-18 Eh and with
  GAMESS HSO1 to its printed precision. With a mean field for the
  two-electron terms, the closed-shell core part reproduces GAMESS HSO2P
  (7e-5 cm-1 in the OH splitting); against the exact HSO2 the difference,
  1.2 cm-1, is the error of the mean-field approximation. See their READMEs.
- `validation/safety/`: sanitizers, valgrind, static analysis and fuzzing of
  `hexpr_rank1_pair`.

## 6. Not covered

- Rank 2 (spin-spin) and the two-electron spin-orbit operator: the chains
  are the two-electron ones of the spin-free path with the final rank 1 or
  2. The normalisation of Sec. 5 holds for them, but there is no C
  evaluator yet; the examples use the spin-orbit mean field for the
  two-electron part.
- Whole-space expressions (the loop-driven `hexpr_expr_build`) for rank 1:
  only the pair evaluator exists.
