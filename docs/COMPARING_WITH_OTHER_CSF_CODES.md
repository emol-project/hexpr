# Comparing hexpr with other CSF codes

> This document is informative. The conventions it relies on are
> normative and are stated in [`c/include/hexpr.h`](../c/include/hexpr.h);
> where the two differ, the header is right.

This is for anyone who has another CSF-based code and wants to know
whether it and hexpr agree. In general the two will not use the same CSF
basis, so their matrix elements cannot be compared one for one even when
both codes are correct. What follows is a comparison that does not
depend on the coupling scheme, the data each side has to produce for it,
and the ways it goes wrong.

No tool is provided. The comparison needs the other code's output, there
is no common format for it, and the method is short enough to write
against whatever format the other code has.
[`validation/csf2det_vs_pyscf.py`](../validation/csf2det_vs_pyscf.py) is
a complete, working instance of one part of it (Level 2a below), with
PySCF as the other side.

---

## Why CSFs cannot be compared directly

A CSF is fixed by four things:

1. its spatial configuration (the occupation of each orbital: 0, 1 or 2),
2. its total spin S,
3. the order and scheme in which the individual spins are coupled,
4. sign conventions.

Only the first two are physical. hexpr uses the genealogical coupling
defined by its DRT walk, taken in orbital order 0, 1, ..., norb-1, with
Condon-Shortley Clebsch-Gordan phases. Other codes use the Serber or
Rumer schemes, or genealogical coupling in a different orbital order, or
the same coupling with different phases.

For a given configuration and S, the CSFs of every scheme span the same
space, and two orthonormal schemes differ by an orthogonal transformation
inside that space. (Rumer functions are not orthogonal; see
"Normalization and orthogonality" below.) Two things follow:

- Individual CSFs, coupling coefficients and matrix elements `H[I,J]`
  differ between two codes that are both correct.
- The transformation between two codes' bases is **block diagonal over
  configurations**. It never mixes CSFs of different spatial
  occupation.

The second point is what makes Levels 1 and 3 below possible.

## Determinants as the common basis

A Slater determinant, once an operator-order convention is fixed, is
specified by its alpha and beta occupations alone. It has no coupling
and no phase freedom. Every CSF code can expand its CSFs into
determinants, and in the determinant basis two codes' results can be
compared element by element.

hexpr's expansion is `hexpr_csf_to_dets` in C and
`Drt.csf_to_dets(index, ms_x2)` in Python, which returns the alpha and
beta masks as `(ndets, nwords)` arrays together with the coefficients;
the other bindings expose the same call. Its conventions, stated in full
in `hexpr.h`, are:

- operators ordered all alpha first, then all beta, each in ascending
  spatial orbital;
- bit `p` of word `w` is spatial orbital `64*w + p`;
- Ms passed as `ms_x2 = 2*Ms`, any value allowed by S;
- every CSF normalized at every Ms;
- zero coefficients emitted, not dropped.

Two conventions are left that both sides must share before any
determinant-level comparison means anything: the numbering of the
orbitals and the operator order. Both are covered under "Pitfalls".

---

## Level 1: invariants

**Needs:** the Hamiltonian in each code's CSF basis, over the same CSF
space and with the same integrals, and the configuration of each CSF.
No determinants.

For hexpr, the configuration of CSF `i` comes from `hexpr_csf_steps`:
steps E, U, D, F give orbital occupations 0, 1, 1, 2.

Group the CSFs of each code by configuration. Because the transformation
between the two bases is orthogonal and block diagonal over
configurations, the following agree between the two codes:

- the eigenvalues of the whole matrix;
- the eigenvalues of each diagonal block `H[c,c]`;
- the singular values of each off-diagonal block `H[c,c']`, and so its
  Frobenius norm.

None of these depends on how either code orders its CSFs within a
configuration. The block-level checks are much stronger than the
eigenvalues alone, and they localize a disagreement to a pair of
configurations.

What Level 1 cannot see is anything that is itself a change of basis
within a configuration. That is exactly the difference between two
coupling schemes, so Level 1 says nothing about either code's CSF
expansion or phase conventions. It is the cheapest level and the one to
run first.

## Level 2: matrices in the determinant basis

**Needs:** for every CSF, its expansion in determinants at one agreed
Ms.

Write each code's expansion as a matrix `U` (`ncsfs x ndets`), one row
per CSF, the columns indexed by the union of the `(alpha, beta)` keys
that occur. If the CSFs are orthonormal, so are the rows of `U`.

### 2a. Against a determinant code

Take `H_det` from a determinant CI program over the same determinants,
and check

```
U H_det U^T  ==  H_csf
```

element by element. This tests one CSF code at a time: a pass means the
code's coupling coefficients and its own determinant expansion are
consistent with each other and with the determinant program.
`validation/csf2det_vs_pyscf.py` does this for hexpr against PySCF, in 14
cases, to a tolerance of `1e-10`.

Eigenvalues are no substitute here. They are invariant under any
orthogonal change of basis inside the CSF space and so are blind to
exactly the signs this level exists to check.

### 2b. Between two CSF codes, with no determinant program

Form, for each code,

```
P = U^T H_csf U          (ndets x ndets)
```

If `H_csf = U H_det U^T`, then `P = Q H_det Q` with `Q = U^T U`, the
projector onto the CSF space. `P` depends on the space and not on the
basis chosen in it, so

```
P_hexpr  ==  P_other
```

element by element, whatever the two coupling schemes. A difference
means one of: a code whose coupling coefficients are inconsistent with
its own determinant expansion, CSF spaces that differ, or a convention
(orbital numbering, operator order, Ms) that the two sides do not share.
Level 2a, run on each side separately, tells the first apart from the
others.

## Level 3: the basis change

**Needs:** the same as Level 2, plus `H_csf` from each code.

```
T = U_hexpr U_other^T    (ncsfs x ncsfs),   T[i,j] = <hexpr CSF i | other CSF j>
```

Check, in this order:

- `T T^T == I`. If it fails, the two CSF spaces differ, or the orbital
  numbering, operator order or Ms does. See "Pitfalls".
- `T` is block diagonal over configurations. An element connecting two
  different configurations means the orbital numbering does not match.
- `H_hexpr == T H_other T^T`.

Then read `T`:

- For a configuration with a single CSF, the block is `+1` or `-1`.
  The pattern of signs is the difference between the two codes' phase
  conventions. It is not an error; record it.
- If the other code is genealogical in the same orbital order, `T` is a
  permutation with signs.
- A dense block **is** the recoupling matrix between the two schemes for
  that open-shell pattern. It comes out numerically, with no 6j or 9j
  symbols.
- `T` does not depend on Ms, provided both codes phase the Ms
  components of a spin state as the ladder operators do under
  Condon-Shortley: `S_-` is the adjoint of `S_+` and lowers both sides
  by the same real factor, so
  `<Phi(S,M-1)|Psi(S,M-1)> = <Phi(S,M)|Psi(S,M)>`. hexpr's expansion,
  built from Clebsch-Gordan coefficients with Condon-Shortley phases,
  is of this kind at every Ms. A sign that changes with Ms therefore
  means the other side phases its Ms components differently.

`T` is cheap. The rows of `U` belonging to one configuration touch only
that configuration's determinants, so `T` can be built block by block,
never forming anything of size `ndets`.

---

## What to ask the other side for

- **The orbitals:** which orbitals are active, in which order along the
  coupling chain, and how they map onto hexpr's `0 .. norb-1`. Which
  orbitals are frozen.
- **The integrals**, or agreement to use one common set. hexpr takes
  `h1eff` and `eri` in chemists' notation `(pq|rs)`, with the constant
  core energy outside; see
  [`docs/PQRS_INDEX_SPEC.md`](PQRS_INDEX_SPEC.md).
- **For every CSF:** its spatial configuration, and its determinant
  expansion at a stated Ms as a list of (alpha bit string, beta bit
  string, coefficient), with the bit-to-orbital mapping stated.
- **The operator-order convention** of their determinants.
- **Ms**, if the expansion is not at an agreed value.
- **The Hamiltonian in their CSF basis**, for Levels 1, 2b and 3.
- **Whether their CSFs are orthonormal.**

---

## Pitfalls

### Orbital numbering

In hexpr, step position `k`, determinant bit `k` and integral index `k`
are all orbital `k`, and the genealogical coupling runs in that order.
Coupling depends on the order: a genealogical code that couples the
same orbitals in a different order -- the reverse order, say --
produces a different basis, and `T` is dense where you would expect
signs only. Singlets are the exception for the reverse order: when the
total S is 0, coupling from the other end gives the same functions up
to order and sign, so `T` is a signed permutation and a reversed
orbital order goes unnoticed until a case with S > 0 is compared. Map
the orbitals explicitly; do not assume.

Symptoms: `T` not block diagonal means the numbering is wrong. `T` block
diagonal but dense where a sign permutation was expected means the
coupling order differs.

### Operator order

A determinant written with the operators interleaved by orbital differs
from hexpr's alpha-then-beta block by a sign that depends on the
determinant. [`c/src/csf/csf2det.c`](../c/src/csf/csf2det.c) gives the
rule. Convert the other side's coefficients to one convention before
building `U`.

Symptom: Level 1 passes, Level 2b fails, and the failing elements follow
the determinants rather than the CSFs.

### Ms

hexpr accepts any Ms allowed by S; many codes expand only at `Ms = S`.
Determinants at different Ms are disjoint, so a mismatch shows up as a
key union with no overlap. Agree on Ms first. `Ms = 0` (or `1/2`) is the
most searching choice: it has the largest expansion and the most room
for an alpha/beta ordering error.

### Normalization and orthogonality

hexpr's CSFs are orthonormal. Rumer and other valence-bond functions are
not. With `S = U_other U_other^T` their overlap matrix:

- Level 1 compares the eigenvalues of the generalized problem
  `H c = E S c`, over the whole space and block by block.
- Level 2b uses `P_other = U_other^T S^-1 H_other S^-1 U_other`.
- Level 3 checks `H_hexpr == T S^-1 H_other S^-1 T^T`.

Or orthonormalize the other side first (for example Loewdin
`S^(-1/2)`), after which everything above applies unchanged.

### Different spaces

Excitation-level restrictions, a different frozen core, a different
reference, or spatial symmetry in the other code (hexpr uses none) all
give CSF spaces that differ. Symptoms: `T T^T != I` and `P_hexpr !=
P_other`.

Restrict both sides to the configurations they have in common. Because
`T` is block diagonal, the block-level checks of Levels 1 and 3 remain
valid configuration by configuration; the whole-matrix eigenvalues and
Level 2b do not. On the hexpr side, take every CSF of the common
configurations -- all of them, or the span changes -- and build a CSFSET
sub-DRT from them with `hexpr_drt_from_csfset`.

### Determinant keys

hexpr emits zero coefficients; most codes drop them. hexpr orders the
determinants of one CSF lexicographically by which open shells carry
alpha; no other code is obliged to do the same. Index the columns of
`U` by the `(alpha, beta)` key, never by position.

### Integrals with accidental zeros

Spatial symmetry makes many integrals vanish, and a coupling coefficient
multiplied by a zero integral is not tested. None of the three levels
needs physical integrals. Real random values with the 8-fold
permutational symmetry of `(pq|rs)`, and a symmetric `h1eff`, exercise
every coefficient.

### The overall sign of a CSF

It is a convention. Never compare an individual coefficient's sign
across codes; compare through `T`.

---

## Cost

A CSF with `n_open` open shells has `C(n_open, n_alpha_open)`
determinants at a given Ms, with `n_alpha_open = (n_open + 2*Ms) / 2`.
At `Ms = 0`, with the number of singlet CSFs per configuration beside
it:

| n_open | determinants per CSF | singlet CSFs per configuration |
|-------:|---------------------:|-------------------------------:|
| 2  | 2      | 1     |
| 4  | 6      | 2     |
| 6  | 20     | 5     |
| 8  | 70     | 14    |
| 10 | 252    | 42    |
| 12 | 924    | 132   |
| 16 | 12870  | 1430  |
| 20 | 184756 | 16796 |

Both columns grow exponentially.

- **Level 1** costs nothing beyond the two Hamiltonians and can be run
  on the full space.
- **Level 3**, built block by block, costs per configuration the number
  of CSFs times the number of determinants, and is practical to twelve
  or more open shells.
- **Level 2b** as written forms dense `ndets x ndets` matrices over the
  whole space and is the expensive one. Evaluate it one pair of
  configurations at a time, or on small spaces only.
- **Level 2a** is limited by the determinant program.

In practice: Level 1 on the space of interest, and Levels 2 and 3 on
small spaces that contain every open-shell pattern you care about. That
is how the validation script is built -- H4 and H6, up to six open
shells, closed shells between open ones -- rather than on a large
molecule.

---

## See also

- [`c/include/hexpr.h`](../c/include/hexpr.h) -- the normative statement
  of the expansion conventions
- [`c/src/csf/csf2det.c`](../c/src/csf/csf2det.c) -- their derivation,
  including the operator-order sign rule
- [`c/README.md`](../c/README.md), "Note: CSF -> determinant expansion"
- [`validation/csf2det_vs_pyscf.py`](../validation/csf2det_vs_pyscf.py)
  -- Level 2a against PySCF
- [`docs/PQRS_INDEX_SPEC.md`](PQRS_INDEX_SPEC.md) -- the integral
  conventions
