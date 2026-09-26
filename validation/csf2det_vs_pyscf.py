"""
csf2det_vs_pyscf.py -- settle the csf2det operator-ordering convention
against an independent determinant-basis code.

WHY THIS EXISTS
---------------
The C tests (c/tests/test_csf2det.c) pin everything about the expansion that
can be checked without an oracle: determinant count, coefficient magnitudes,
unit norm, Sz, <I|J> = delta_IJ across every Ms, and S^2. None of those
constrains the SIGN of a coefficient against anything outside hexpr: the first
five are quadratic forms or bit counts, and S^2 is a product of S_- and S_+ in
which a sign common to a whole CSF cancels.

So this script settles the rest: whether hexpr_csf_to_dets puts its
coefficients on the determinants a determinant-basis program would recognise,
with the signs that program would use.

It reaches csf2det through the public Python binding rather than a private
ctypes layer, so the binding is verified here too, at no extra cost.

WHAT IS COMPARED
----------------
Not eigenvalues. Eigenvalues are invariant under any orthogonal change of
basis inside the CSF space, so they are blind to exactly the sign in
question. What is compared is the matrix:

    U            (ncsfs x ndets)   from hexpr_csf_to_dets
    H_det        (ndets x ndets)   from PySCF, determinant basis
    H_csf        (ncsfs x ncsfs)   from hexpr, CSF basis

    U H_det U^T  must equal  H_csf  element by element.

U is not square, so H_det is projected into the CSF space. U's rows are
already known to be orthonormal from the C tests, so the untested content of
U is which determinant each column refers to and with what sign -- which is
what this product exercises.

ON THE MOLECULES
----------------
The geometries are chosen for numerical convenience, not chemistry. What the
comparison needs from a system is CSF structure -- how many open shells, how
closed shells interleave with them, how long the coupling chain gets -- and
integrals that are not accidentally degenerate. It does not need the system
to be interesting, or even bound.

H6 is therefore three H2 molecules at their own bond length, set far apart
from each other. That makes the RHF converge without fuss while still giving
six orbitals and, through SOCI/FOCI, CSFs with up to six open shells.

ONE SET OF INTEGRALS PER GEOMETRY. The RHF is run once, closed-shell, and its
MO integrals are used for the singlet and triplet DRTs alike. The SCF is only
a way to obtain a reasonable h1eff and eri; hexpr's expression path and
csf2det are what is under test, and both sides of every comparison are fed the
same numbers, so the reference determinant used to generate the orbitals is
irrelevant. PySCF's side is told the spin sector through nelec=(na,nb), not
through the orbitals.

RUN
---
    source $HOME/env/bin/activate
    cd <worktree root>
    HEXPR_LIBRARY=$PWD/c/libhexpr.so \
    LD_LIBRARY_PATH=$PWD/c \
    PYTHONPATH=$PWD/python \
        python validation/csf2det_vs_pyscf.py
"""

import os
import sys
import time

import numpy as np

TOL = 1e-10
STEPCH = ".ud2"          # 0=E 1=U 2=D 3=F


# ===========================================================================
# Which library is under test.
#
# csf2det is now reached through the public Python binding, so this script
# verifies the binding as well as the C code. Everything below used to be a
# private ctypes layer; it was written that way on purpose while the
# convention was still unsettled, because a binding written before the
# convention was right would have spread a wrong convention across five
# languages -- and the convention was in fact wrong at first.
#
# The binding's own search order is HEXPR_LIBRARY, then the in-tree build,
# then the system linker path. The installed copy reports the same version
# string as the in-tree one, so picking up the wrong library is invisible at
# run time. Hence the check below: say out loud which file is being tested.
# ===========================================================================
def report_library():
    import hexpr._ctypes as _c
    path = os.path.abspath(_c.lib._name)
    print("library under test: %s" % path)
    if not os.path.exists(path):
        sys.exit("the resolved library does not exist: %s" % path)
    return path


def steps_str(step):
    return "".join(STEPCH[s] for s in step)


def has_closed_between_open(step):
    """True if an F sits between two open shells -- the configuration whose
    treatment in csf2det.c rests on a textbook result, not on measurement."""
    pos = [k for k, s in enumerate(step) if s in (1, 2)]
    if len(pos) < 2:
        return False
    return any(step[k] == 3 for k in range(pos[0], pos[-1]))


def build_u(drt, ms_x2):
    """Return (U, det_keys, steps): U[i, d] is the coefficient of determinant
    det_keys[d] in CSF i; det_keys are (alpha_mask, beta_mask)."""
    norb = drt.norb
    ncsfs = drt.ncsfs
    if (norb + 63) // 64 != 1:
        sys.exit("this script assumes norb <= 64")

    steps = [drt.csf_steps(i).tolist() for i in range(ncsfs)]

    cols = {}
    rows = []
    for i in range(ncsfs):
        a, b, c = drt.csf_to_dets(i, ms_x2)
        n = len(c)
        row = {}
        for k in range(n):
            key = (int(a[k, 0]), int(b[k, 0]))
            if key not in cols:
                cols[key] = len(cols)
            row[cols[key]] = row.get(cols[key], 0.0) + c[k]
        rows.append(row)

    U = np.zeros((ncsfs, len(cols)))
    for i, row in enumerate(rows):
        for d, v in row.items():
            U[i, d] = v
    det_keys = [None] * len(cols)
    for key, d in cols.items():
        det_keys[d] = key
    return U, det_keys, steps


# ===========================================================================
# H in the CSF basis, through the public Python binding.
# ===========================================================================
def tri_inv(t):
    p = int(np.ceil((-1.0 + np.sqrt(1.0 + 8.0 * t)) / 2.0))
    if p * (p - 1) // 2 >= t:
        p -= 1
    return p, t - p * (p - 1) // 2


def assemble_h(ij_list, pqrs_list, coef_list, bra, ket,
               ncsfs, norb, h1eff, eri_act):
    """Shared assembly loop. When bra/ket are None, ij is the packed triangle
    index and is unpacked here (the Drt/Expression path); otherwise bra/ket
    are used directly (the csf_subspace_eval path)."""
    n_oel = norb * (norb + 1) // 2
    H = np.zeros((ncsfs, ncsfs))
    for k in range(len(pqrs_list)):
        if bra is not None:
            I, J = int(bra[k]), int(ket[k])
        else:
            ij = int(ij_list[k])
            I = int((np.sqrt(8 * ij + 1) - 1) // 2)
            J = ij - I * (I + 1) // 2
        c = coef_list[k]
        t = int(pqrs_list[k])
        if t < n_oel:
            p, q = tri_inv(t + 1)
            v = h1eff[p - 1, q - 1]
        else:
            a_, b_ = tri_inv(t - n_oel + 1)
            p, q = tri_inv(a_)
            r, s = tri_inv(b_)
            v = eri_act[p - 1, q - 1, r - 1, s - 1]
        H[I, J] += c * v
        if I != J:
            H[J, I] += c * v
    return H


# ===========================================================================
# H in the determinant basis, from PySCF.
# ===========================================================================
def build_h_det(det_keys, norb, nelec, h1eff, eri_act):
    from pyscf import fci
    from pyscf.fci import cistring

    na, nb = nelec
    sa = cistring.make_strings(range(norb), na)
    sb = cistring.make_strings(range(norb), nb)
    index = {}
    for ia, x in enumerate(sa):
        for ib, y in enumerate(sb):
            index[(int(x), int(y))] = (ia, ib)

    missing = [k for k in det_keys if k not in index]
    if missing:
        sys.exit("csf2det produced determinants PySCF does not enumerate "
                 "at nelec=%s: %s" % ((nelec,), missing[:5]))

    dim = (len(sa), len(sb))
    h2e = fci.direct_spin1.absorb_h1e(h1eff, eri_act, norb, nelec, .5)

    n = len(det_keys)
    H = np.zeros((n, n))
    for j, kj in enumerate(det_keys):
        ia, ib = index[kj]
        civec = np.zeros(dim)
        civec[ia, ib] = 1.0
        hc = fci.direct_spin1.contract_2e(h2e, civec, norb, nelec)
        for i, ki in enumerate(det_keys):
            ja, jb = index[ki]
            H[i, j] = hc[ja, jb]
    return H


# ===========================================================================
# Integrals. One closed-shell RHF per geometry; the same MO integrals serve
# every DRT built on it, whatever its spin.
# ===========================================================================
def integrals_h2():
    """Hard-coded, from example_ci/h2_3csfs/run.py, so the two sides cannot
    disagree about the integrals rather than about the basis."""
    h1 = np.array([[-1.252463573564898, 0.0],
                   [0.0, -0.475948715220964]])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0] = 0.674488766356838
    eri[1, 1, 1, 1] = 0.697393767423027
    eri[0, 0, 1, 1] = eri[1, 1, 0, 0] = 0.663468096423568
    eri[0, 1, 0, 1] = eri[0, 1, 1, 0] = \
        eri[1, 0, 0, 1] = eri[1, 0, 1, 0] = 0.181288808211496
    return h1, eri


def cas_integrals(atom, ncas, nelecas):
    """Closed-shell RHF, then a CASCI harness for the MO integrals."""
    from pyscf import gto, scf, mcscf, ao2mo
    mol = gto.Mole()
    mol.atom = atom
    mol.basis = "sto-3g"
    mol.charge = 0
    mol.spin = 0
    mol.unit = "Angstrom"
    mol.symmetry = False
    mol.verbose = 0
    mol.build()

    mf = scf.RHF(mol)
    mf.kernel()
    if not mf.converged:
        sys.exit("RHF did not converge for:\n%s" % atom)

    mc = mcscf.CASCI(mf, ncas, nelecas)
    h1, _ecore = mc.get_h1eff()
    mo = mf.mo_coeff[:, mc.ncore:mc.ncore + mc.ncas]
    eri = ao2mo.kernel(mol, mo, compact=False).reshape(
        ncas, ncas, ncas, ncas)
    return h1, eri


# H4, linear: four orbitals is the smallest size that admits a closed shell
# between two open shells.
ATOM_H4 = """
H 0.0 0.0 0.00
H 0.0 0.0 1.00
H 0.0 0.0 2.00
H 0.0 0.0 3.00
"""

# H6 as three well-separated H2 molecules at their own bond length. The RHF
# converges without fuss; the point of the system is six orbitals and long
# coupling chains, not chemistry.
ATOM_H6 = """
H 0.0 0.0 0.000
H 0.0 0.0 0.741
H 0.0 4.0 0.000
H 0.0 4.0 0.741
H 0.0 8.0 0.000
H 0.0 8.0 0.741
"""


# ===========================================================================
# One case
# ===========================================================================
def run_case(label, drt, ncsfs, norb, ms_x2, nelec,
             h1eff, eri_act, H_csf):
    print("\n" + "-" * 70)
    print("%s   ms_x2=%+d  nelec=%s" % (label, ms_x2, (nelec,)))
    print("-" * 70)

    t0 = time.time()
    U, det_keys, steps = build_u(drt, ms_x2)
    print("  ncsfs=%d  determinants spanned=%d  norb=%d"
          % (ncsfs, len(det_keys), norb))

    nf = sum(1 for s in steps if has_closed_between_open(s))
    nopen_max = max((sum(1 for x in s if x in (1, 2)) for s in steps),
                    default=0)
    print("  CSFs with a closed shell BETWEEN open shells: %d" % nf)
    print("  largest number of open shells in any CSF:     %d" % nopen_max)
    if nf == 0:
        print("  NOTE: this case does NOT exercise interleaved closed shells")

    err_g = np.abs(U @ U.T - np.eye(ncsfs)).max()
    print("  max |U U^T - I|           = %.3e" % err_g)

    H_det = build_h_det(det_keys, norb, nelec, h1eff, eri_act)
    proj = U @ H_det @ U.T
    diff = proj - H_csf
    err = np.abs(diff).max()
    scale = max(np.abs(H_csf).max(), 1e-30)
    print("  max |U H_det U^T - H_csf| = %.3e   (relative %.3e)"
          % (err, err / scale))
    print("  [%.1fs]" % (time.time() - t0))

    if err > TOL:
        i, j = np.unravel_index(np.abs(diff).argmax(), diff.shape)
        print("  WORST element [%d, %d]: hexpr %+.12f  pyscf %+.12f"
              % (i, j, H_csf[i, j], proj[i, j]))
        print("    csf %d = %s  (f between open: %s)"
              % (i, steps_str(steps[i]), has_closed_between_open(steps[i])))
        print("    csf %d = %s  (f between open: %s)"
              % (j, steps_str(steps[j]), has_closed_between_open(steps[j])))
        bad = np.abs(diff) > TOL
        involved = sorted(set(np.where(bad.any(axis=1))[0].tolist()))
        withf = [k for k in involved if has_closed_between_open(steps[k])]
        print("    %d of %d CSFs appear in a bad element; %d of those have an "
              "f between open shells" % (len(involved), ncsfs, len(withf)))
        if involved and len(withf) == len(involved):
            print("    EVERY involved CSF has that pattern -- look at the "
                  "closed-shell handling in csf2det.c, not det_block_sign()")
        elif involved and not withf:
            print("    NONE of them has that pattern -- the closed-shell "
                  "handling is not the cause; look at det_block_sign()")

    return err, err_g, (nf > 0), nopen_max


def run_drt_case(hexpr, label, kind, nve, spin_x2,
                 ncore, nval, next_, ms_list, h1, eri, results, state):
    """Build the DRT, assemble H_csf once, compare at each Ms.

    One DRT now serves both sides: the expression path and csf2det both go
    through the binding, so there is no second library handle that could
    resolve to a different libhexpr.so.
    """
    try:
        drt = hexpr.Drt(kind, nve=nve, spin_x2=spin_x2,
                        ncore=ncore, nval=nval, next_=next_)
    except Exception as e:                      # noqa: BLE001
        print("\n%s: SKIP, Drt() failed: %s" % (label, e))
        return
    expr = hexpr.Expression(drt)
    ij_a, pqrs_a, coef_a = expr.terms_full()
    ncsfs, norb = drt.ncsfs, drt.norb
    H_csf = assemble_h(ij_a, pqrs_a, coef_a, None, None,
                       ncsfs, norb, h1, eri)

    for ms in ms_list:
        na, nb = (nve + ms) // 2, (nve - ms) // 2
        r = run_case(label, drt, ncsfs, norb, ms, (na, nb),
                     h1, eri, H_csf)
        results.append(("%s ms=%+d" % (label, ms), r[0], r[1]))
        state["f"] = state["f"] or r[2]
        state["open"] = max(state["open"], r[3])


# ===========================================================================
def main():
    import hexpr

    hexpr.init()

    print("=" * 70)
    print("csf2det vs PySCF -- matrix comparison in the determinant basis")
    print("=" * 70)
    report_library()

    results = []
    state = {"f": False, "open": 0}

    # ---------------- (a) H2 CAS(2,2), explicit CSF set ----------------
    h1, eri = integrals_h2()
    csfs = np.array([[3, 0], [0, 3], [1, 2]], dtype=np.uint8)
    drt = hexpr.Drt.from_csfset(csfs, norb=2)

    for i in range(3):
        got = drt.csf_steps(i).tolist()
        if got != csfs[i].tolist():
            sys.exit("index convention mismatch at %d: %s vs %s"
                     % (i, got, csfs[i].tolist()))
    print("\nindex convention: csf_steps(i) == input row i, confirmed")

    ij, pqrs, coef, bra, ket = hexpr.csf_subspace_eval(csfs, which='full')
    H_csf = assemble_h(ij, pqrs, coef, bra, ket, 3, 2, h1, eri)
    r = run_case("(a) H2 CAS(2,2) singlet, explicit CSF set",
                 drt, 3, 2, 0, (1, 1), h1, eri, H_csf)
    results.append(("(a) H2 singlet", r[0], r[1]))
    state["f"] = state["f"] or r[2]
    state["open"] = max(state["open"], r[3])
    del drt

    # ---------------- (b) H4 CAS(4,4) ----------------
    print("\n[H4 CAS(4,4): one closed-shell RHF, integrals shared by every DRT]")
    h1, eri = cas_integrals(ATOM_H4, 4, 4)
    for label, kind, spin_x2, ms_list in [
        ("(b) H4 SOCI singlet", "soci", 0, [0]),
        ("(b) H4 FOCI singlet", "foci", 0, [0]),
        ("(b) H4 SOCI triplet", "soci", 2, [2, 0, -2]),
    ]:
        run_drt_case(hexpr, label, kind, 4, spin_x2,
                     1, 2, 1, ms_list, h1, eri, results, state)

    # ---------------- (c) H6 CAS(6,6) ----------------
    # Three separated H2. Six orbitals, up to six open shells, and closed
    # shells interleaved among them in more ways than H4 allows.
    print("\n[H6 CAS(6,6), three separated H2: one closed-shell RHF, "
          "integrals shared by every DRT]")
    h1, eri = cas_integrals(ATOM_H6, 6, 6)
    for label, kind, spin_x2, ms_list in [
        ("(c) H6 SOCI singlet", "soci", 0, [0]),
        ("(c) H6 FOCI singlet", "foci", 0, [0]),
        ("(c) H6 SOCI triplet", "soci", 2, [2, 0, -2]),
        ("(c) H6 FOCI triplet", "foci", 2, [2, 0, -2]),
    ]:
        run_drt_case(hexpr, label, kind, 6, spin_x2,
                     1, 4, 1, ms_list, h1, eri, results, state)

    # ---------------- summary ----------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    worst = 0.0
    for name, err, err_g in results:
        flag = "ok  " if (err < TOL and err_g < TOL) else "FAIL"
        print("  %-40s %s  max|diff| = %.3e" % (name, flag, err))
        worst = max(worst, err)

    print("\ncoverage actually reached by this run:")
    print("  closed shell between open shells : %s"
          % ("yes" if state["f"] else "NO -- still unverified"))
    print("  largest open-shell count         : %d" % state["open"])
    print("  negative Ms                      : %s"
          % ("yes" if any("ms=-" in n for n, _, _ in results) else "no"))

    if worst >= TOL:
        print("\nRESULT: MISMATCH -- see the worst-element diagnostics above.")
        return 1
    if state["f"] and state["open"] >= 5:
        print("\nRESULT: the shipped block ordering reproduces the "
              "determinant-basis matrix in every case above, including")
        print("interleaved closed shells and coupling chains through %d open "
              "shells." % state["open"])
        return 0
    print("\nRESULT: every case matched, but the coverage lines above show "
          "what was not reached.")
    print("Treat the untouched configurations as still unverified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
