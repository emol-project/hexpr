"""
OH X 2Pi spin-orbit splitting from PySCF integrals + hexpr.

OH at R = 0.9697 A (given in Bohr) along z, O 1s core, CAS(7,5) on the
valence 2s 2p of O and 1s of H. Orbitals from a state-averaged CASSCF over
the two components of X 2Pi, so that they are equivalent. The one-electron
Breit-Pauli spin-orbit operator (true or effective nuclear charges) is
built between the doublet CSFs with hexpr.rank1_pair and checked against
the determinant construction; the state-interaction matrix over the lowest
doublet roots is diagonalized, each level is named by Omega = |<J_z>| with
J_z = L_z + S_z, and A = E(2Pi_3/2) - E(2Pi_1/2) is reported (experiment:
A = -139.2 cm-1, inverted).

    python run_pyscf_oh.py [--basis cc-pvtz] [--zeff O=5.6,H=1.0] [--nroot 4]
"""
import argparse, os, sys
import numpy as np
from pyscf import gto, scf, mcscf, ao2mo, fci

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from so_tools import (HARTREE_CM, CasSpace, DetSpace, so_blocks, so_integrals,
                      spherical, hexpr)

ap = argparse.ArgumentParser()
ap.add_argument("--basis", default="cc-pvtz")
ap.add_argument("--zeff", help="effective charges, e.g. O=5.6,H=1.0 (default: true)")
ap.add_argument("--nroot", type=int, default=4, help="doublet roots in the SI matrix")
ap.add_argument("--gamess", help="GAMESS SO-CI log (gamess/oh_so_*.log) to compare A with")
ap.add_argument("--somf", choices=["core", "full"],
                help="add the two-electron mean field: 'core' = O 1s density only "
                     "(exact core-active part, GAMESS HSO2P), 'full' = state-averaged "
                     "CASSCF density (SOMF)")
args = ap.parse_args()
zeff = None
if args.zeff:
    zeff = {k: float(v) for k, v in (kv.split("=") for kv in args.zeff.split(","))}

R_BOHR = 1.8324674229907405          # 0.9697 A with PySCF's Bohr radius
mol = gto.M(atom=f"O 0 0 0; H 0 0 {R_BOHR}", basis=args.basis, unit="Bohr",
            spin=1, verbose=0)
mf = scf.ROHF(mol)
mf.conv_tol = 1e-12
mf.kernel()

ncas, nelecas = 5, 7
mc = mcscf.CASSCF(mf, ncas, (4, 3))
mc.fcisolver = fci.addons.fix_spin_(fci.direct_spin1.FCI(mol), ss=0.75)
mc = mc.state_average_([0.5, 0.5])
mc.conv_tol = 1e-12
mc.conv_tol_grad = 1e-7
mc.kernel()
print(f"basis {args.basis}, nao {mol.nao}; ROHF {mf.e_tot:.10f}; "
      f"SA-CASSCF(7,5) Pi components {mc.e_states[0]:.10f} {mc.e_states[1]:.10f}")

ci = mcscf.CASCI(mf, ncas, nelecas)
ci.mo_coeff = mc.mo_coeff
h1, ecore = ci.get_h1eff()
C_act = mc.mo_coeff[:, ci.ncore:ci.ncore + ncas]
eri = ao2mo.kernel(mol, C_act, compact=False).reshape((ncas,) * 4)

hexpr.init()
dbl = CasSpace(ncas, nelecas, 1, h1, eri)
print(f"doublet CAS: {dbl.drt.ncsfs} CSFs, lowest roots (Eh) "
      + " ".join(f"{ecore + e:.10f}" for e in dbl.e[:args.nroot]))

D_somf = None
if args.somf == "core":
    C_core = mc.mo_coeff[:, :ci.ncore]
    D_somf = 2.0 * C_core @ C_core.T
elif args.somf == "full":
    D_somf = mc.make_rdm1()                    # state-averaged, spin-summed, AO
    print(f"SOMF density: {np.einsum('ij,ji->', D_somf, mf.get_ovlp()):.6f} electrons")
h_cart = so_integrals(mol, C_act, zeff, D_somf)
X = so_blocks(dbl, dbl, spherical(h_cart), ncas)

# independent check in the determinant basis
dets = DetSpace(ncas, nelecas)
Hso_det = dets.h_so(h_cart)
U = {ms: dets.csf_matrix(dbl, ms) for ms in (-1, 1)}
worst = max(np.abs(X[mb, mk] - U[mb].T @ Hso_det @ U[mk]).max()
            for mb in (-1, 1) for mk in (-1, 1))
print(f"doublet-doublet CSF spin-orbit matrix, hexpr vs determinants: max |diff| = {worst:.1e} Eh")

# state interaction over the lowest doublet roots
basis = [(n, ms) for n in range(args.nroot) for ms in (-1, 1)]
Hsi = np.zeros((len(basis), len(basis)), dtype=complex)
for a, (nb, mb) in enumerate(basis):
    Hsi[a, a] = ecore + dbl.e[nb]
    for b, (nk, mk) in enumerate(basis):
        Hsi[a, b] += dbl.c[:, nb] @ X[mb, mk] @ dbl.c[:, nk]
assert np.abs(Hsi - Hsi.conj().T).max() < 1e-12
lev, vec = np.linalg.eigh(Hsi)

# J_z = L_z + S_z in the same basis; L_z = -1j * (x d/dy - y d/dx), origin on the axis
mol.set_common_orig((0, 0, 0))
Lz_act = -1j * np.einsum("ip,ij,jq->pq", C_act, mol.intor("int1e_cg_irxp", comp=3)[2], C_act)
Lz_det = dets.spin_free(Lz_act)
Jz = np.zeros_like(Hsi)
for a, (nb, mb) in enumerate(basis):
    for b, (nk, mk) in enumerate(basis):
        if mb == mk:
            Jz[a, b] = (U[mb] @ dbl.c[:, nb]) @ Lz_det @ (U[mk] @ dbl.c[:, nk])
    Jz[a, a] += mb / 2
# Kramers pairs are degenerate, so eigh mixes them arbitrarily: diagonalize
# J_z inside each degenerate group of levels (J_z commutes with H for a
# linear molecule), and name the group by |J_z|.
omega = np.zeros(len(lev))
i = 0
while i < len(lev):
    j = i
    while j + 1 < len(lev) and lev[j + 1] - lev[i] < 1e-9:
        j += 1
    V = vec[:, i:j + 1]
    omega[i:j + 1] = np.abs(np.linalg.eigvalsh(V.conj().T @ Jz @ V))
    i = j + 1
print("\nspin-orbit levels (cm-1 above the lowest), Omega = |J_z|:")
for e, w in zip(lev, omega):
    print(f"  {(e - lev[0]) * HARTREE_CM:12.4f}   Omega = {w:.4f}")
x = [(e, w) for e, w in zip(lev[:4], omega[:4])]
e32 = np.mean([e for e, w in x if abs(w - 1.5) < 0.05])
e12 = np.mean([e for e, w in x if abs(w - 0.5) < 0.05])
label = "true charges" if zeff is None else "Z_eff " + args.zeff
if args.somf:
    label += f", two-electron mean field ({args.somf})"
A_mine = (e32 - e12) * HARTREE_CM
print(f"\nA = E(2Pi_3/2) - E(2Pi_1/2) = {A_mine:.4f} cm-1  "
      f"({label}; experiment -139.2 cm-1)")

if args.gamess:
    # GAMESS SO-CI log: the four lowest spin-orbit levels and their <JZ> labels.
    # Its hartree -> cm-1 factor is 1.00001703 x CODATA 2018 (ch2_so/README.md).
    import re
    txt = open(args.gamess).read()
    rel = [float(x) for x in re.findall(r"RELATIVE E=\s+(-?\d+\.\d+) CM-1", txt)][:4]
    jz = [max(((float(j), float(w)) for j, w in re.findall(r"(-?\d\.\d)\(\s*(\d+\.\d)%\)", line)),
              key=lambda t: t[1])[0]
          for line in re.findall(r"STATE\s+\d+: <JZ>=(.*)", txt)][:4]
    g32 = np.mean([e for e, j in zip(rel, jz) if abs(abs(j) - 1.5) < 0.01])
    g12 = np.mean([e for e, j in zip(rel, jz) if abs(abs(j) - 0.5) < 0.01])
    A_g = (g32 - g12) / 1.00001703
    print(f"GAMESS ({os.path.basename(args.gamess)}): A = {A_g:.4f} cm-1 in CODATA 2018 units "
          f"({g32 - g12:.3f} as printed); difference {A_mine - A_g:+.1e} cm-1")
