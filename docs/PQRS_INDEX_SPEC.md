# PQRS index encoding specification

> This specification is normative for hexpr v1.0.0 and later.

This document describes how hexpr encodes Hamiltonian matrix element indices
in the `ij` and `pqrs` integer arrays returned by `hexpr_expr_read` (C API)
or the language-binding equivalents: `terms_full` / `terms_one` in
Python, `read_full` / `read_one` in Julia and Ruby, and
`hexpr_expr_read_full` / `hexpr_expr_read_one` in Fortran.

---

## Source

The C implementation in `c/src/loopgen/expr.c` is authoritative.
This encoding was originally developed under JSPS KAKENHI
Grant Number JP16K05660 (PI: H. Honda, Kyushu University,
FY2016-2017). See the accompanying methodology paper (DOI
placeholder until preprint is up) for the formal description.

Two helper functions appear throughout:

```
triangle(a, b) = max(a,b)*(max(a,b)+1)//2 + min(a,b)   # used for ij
tri(p, q)      = max(p,q)*(max(p,q)-1)//2 + min(p,q)   # used for pqrs
```

Note: `triangle` and `tri` both symmetrize the pair (order-independent), but
they use different formulas.

---

## `ij` -- CSF pair encoding

Each term contributes to one Hamiltonian matrix element `H[I, J]`.

```
ij = triangle(wI, wJ)
```

where `wI` and `wJ` are the Shavitt path weights of CSF `I` and CSF `J`
respectively.  In hexpr's DRT the path weight of CSF `i` (0-based) equals `i`
(the ordering is by path weight).  Therefore `ij = triangle(I, J)`.

### Properties

- `ij` is 0-based (the Python API does NOT add 1; the text dump does).
- `triangle` is symmetric: `triangle(I,J) = triangle(J,I)`.  Each
  off-diagonal pair (I != J) appears **once** in the expression; the caller
  must symmetrize: `H[I,J] += val` and `H[J,I] += val`.
- Diagonal elements (I = J) appear once; no extra symmetrization needed.
- `I >= J` after decoding (because `triangle` puts the larger index first).

### `ij` in the per-pair, block and subspace APIs

`hexpr_pair_eval`, `hexpr_block_eval` and `hexpr_subspace_eval` take
INPUT-ORDER CSF indices, but the `ij` they emit follows the definition
above: it is `triangle` of the two Shavitt path weights, which is the
CANONICAL index. For SOCI and FOCI the two coincide and decoding `ij`
gives back the indices you passed. For a CSFSET sub-DRT whose canonical
walk reorders, extends or dedups the input, they differ, and decoding
`ij` gives canonical positions that are not valid arguments to those
calls.

Use the `bra_pos` / `ket_pos` outputs to attribute a term to a pair:
they are positions within the index arrays you supplied. `ij` remains
the row address of the matrix element itself, which is what
`hexpr_expr_read` over a full DRT needs and what the H-matrix
construction below uses.

### Decode

```python
import math

def decode_ij(ij):
    """Return (I, J) with I >= J >= 0."""
    I = int((math.sqrt(8*ij + 1) - 1) / 2)   # floor
    J = ij - I*(I+1)//2
    return I, J
```

Verification: `ij=0 -> (0,0)`, `ij=1 -> (1,0)`, `ij=2 -> (1,1)`, `ij=3 -> (2,0)`.

---

## `pqrs` -- orbital index encoding

`pqrs` encodes either a 1-electron or a 2-electron orbital quadruple.

Let `norb` = number of active orbitals in the expression,
and `n_oel = norb*(norb+1)//2`.

### Threshold

| Condition | Type |
|---|---|
| `pqrs < n_oel` | 1-electron: orbital pair `(p, q)` |
| `pqrs >= n_oel` | 2-electron: orbital quadruple `(p, q, r, s)` |

### Helper: `tri` inverse

Both 1-electron and 2-electron decodes use the inverse of `tri`:

```python
def tri_inv(t):
    """Decode tri(p,q) value t to 1-based (p, q) with p >= q >= 1."""
    p = math.ceil((-1 + math.sqrt(1 + 8*t)) / 2)
    q = t - p*(p-1)//2
    return p, q
```

### 1-electron decode

```
pqrs = tri(p, q) - 1          (p, q are 1-based orbital levels)
```

Decode:

```python
def decode_1e(pqrs):
    """Return (p, q) 0-based."""
    p1, q1 = tri_inv(pqrs + 1)   # 1-based
    return p1 - 1, q1 - 1        # 0-based
```

The integral is `h1eff[p, q]` (symmetric; `h[p,q] = h[q,p]` for real MOs).

### 2-electron decode

```
pqrs = tri( tri(p,q), tri(r,s) ) - 1 + n_oel
```

Decode:

```python
def decode_2e(pqrs, norb):
    """Return (p, q, r, s) 0-based."""
    n_oel = norb*(norb+1)//2
    outer = pqrs - n_oel + 1          # = tri(tri(p,q), tri(r,s))
    a1, b1 = tri_inv(outer)            # a1=tri(p,q) >= b1=tri(r,s)
    p1, q1 = tri_inv(a1)               # 1-based
    r1, s1 = tri_inv(b1)               # 1-based
    return p1-1, q1-1, r1-1, s1-1     # 0-based
```

The integral is `eri[p, q, r, s]` in **chemists' notation** `(pq|rs)`:

```
(pq|rs) = ∫∫ φ_p*(r1) φ_q(r1) (1/r12) φ_r*(r2) φ_s(r2) dr1 dr2
```

PySCF's `ao2mo.kernel(..., compact=False)` returns this convention directly.

### 8-fold symmetry

Both `tri` and the outer `tri` are symmetric, so the 4-index quadruple
satisfies `(pq|rs) = (qp|rs) = (pq|sr) = (qp|sr) = (rs|pq) = ...`
The decode is unique up to these equivalences (which are all the same value).

---

## H matrix construction

```python
H = np.zeros((ncsfs, ncsfs))
for k in range(len(ij_arr)):
    ij_val   = int(ij_arr[k])
    pqrs_val = int(pqrs_arr[k])
    c        = float(coef_arr[k])

    I, J = decode_ij(ij_val)

    if pqrs_val < n_oel:
        p, q = decode_1e(pqrs_val)
        val = h1eff[p, q]
    else:
        p, q, r, s = decode_2e(pqrs_val, norb)
        val = eri_act[p, q, r, s]

    H[I, J] += c * val
    if I != J:
        H[J, I] += c * val          # symmetrize

E_active = np.linalg.eigh(H)[0]
E_total  = ecore + E_active[0]      # ecore from PySCF mc.get_h1eff()
```

---

## Verification on `s_full_4o4e` (norb=4, ncsfs=20)

CSF 0 is `[e,e,f,f]` -- orbitals 2 and 3 doubly occupied (0-based).
Diagonal element `ij=0` -> `(I,J)=(0,0)`.  From the reference file:

| pqrs | Type | Decode | coef | Physical meaning |
|---|---|---|---|---|
| 5 | 1e | (2,2) | +2 | 2·h[2,2] |
| 9 | 1e | (3,3) | +2 | 2·h[3,3] |
| 30 | 2e | (2,2,2,2) | +1 | J₂₂ = (22\|22) |
| 60 | 2e | (3,3,2,2) | +4 | 4·J₂₃ = 4·(22\|33) |
| 54 | 2e | (3,2,3,2) | −2 | −2·K₂₃ = −2·(23\|32) |
| 64 | 2e | (3,3,3,3) | +1 | J₃₃ = (33\|33) |

This matches the standard RHF energy formula for two doubly-occupied orbitals
`m=2, n=3`:

```
E = 2·h_mm + 2·h_nn + J_mm + 4·J_mn − 2·K_mn + J_nn
```

All six entries decoded correctly.  The convention is confirmed: chemists'
`(pq|rs)` = `eri[p, q, r, s]`.

---

## Orbital indexing in PySCF

For `s_full_6o6e` (H2O / 3-21G, CASCI(6,6)):

- `norb = 6` (active orbitals)
- `n_oel = 21`
- Decoded indices (0-based) run `0..5`, mapping directly into
  `h1eff[0..5, 0..5]` and `eri_act[0..5, 0..5, 0..5, 0..5]`
  from `mc.get_h1eff()` and `ao2mo.kernel(mol, mo_act, compact=False)`.
- `ecore` from `mc.get_h1eff()` = `E_NN + E_frozen` (the constant part
  outside the active-space CI).

For H2O / 6-31G(d) SOCI: same formulas, `norb = 18`, `n_oel = 171`.
The one frozen orbital (O 1s) is handled by PySCF's `mc.ncore=1`;
hexpr's expression then spans the 18-orbital active space.
