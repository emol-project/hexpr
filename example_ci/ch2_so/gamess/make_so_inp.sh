#!/bin/bash
# Build ch2_so_<zeftyp>.inp for GAMESS RUNTYP=TRANSITN (SO-CI, OPERAT=HSO1):
# singlet and triplet CAS(6,6) on the RHF orbitals of ch2_rhf (the $VEC of its
# .dat in GAMESS's USERSCR). Usage: ./make_so_inp.sh [TRUE|3-21G]
# Roots: the C1 Davidson solver converges only to states of the irreps present
# in its start vectors, so it skips some C2v singlets (exact roots 2, 4, ...);
# run_pyscf_so.py matches GAMESS states to the exact roots by energy. 8 roots
# are kept per spin. IROOTS = NSTATE: with NSTATE > IROOTS (8 vs 3) this
# GAMESS (2024 R2 P1) stopped in the SOC step with "I/O past end of record".
# SYMTOL: the default 1e-4 leaves up to 2e-3 cm-1 error in HSO elements.
set -e
Z=${1:-TRUE}
DAT=${DAT:-$HOME/gamess/restart/ch2_rhf.dat}
tag=$(echo "$Z" | tr 'A-Z-' 'a-z_')
out=ch2_so_${tag}.inp
cat > "$out" <<EOF
 \$CONTRL SCFTYP=NONE CITYP=GUGA RUNTYP=TRANSITN ISPHER=1 UNITS=BOHR MULT=1 \$END
 \$SYSTEM MWORDS=50 \$END
 \$BASIS GBASIS=CCD \$END
 \$GUGDIA CVGTOL=1.0E-10 \$END
 \$TRANST OPERAT=HSO1 ZEFTYP=$Z NUMCI=2 NFZC=1 NOCC=7 SYMTOL=1.0E-10
         IROOTS(1)=8,8 NSTATE(1)=8,8 \$END
 \$DRT1 GROUP=C1 FORS=.TRUE. NFZC=1 NDOC=3 NVAL=3 \$END
 \$DRT2 GROUP=C1 FORS=.TRUE. NFZC=1 NDOC=2 NALP=2 NVAL=2 \$END
 \$DATA
CH2 cc-pVDZ CAS(6,6) singlet/triplet, HSO1 ZEFTYP=$Z
C1
C   6.0   0.0000000000   0.0000000000   0.1899174755
H   1.0   0.0000000000   1.8721516716  -1.1400717710
H   1.0   0.0000000000  -1.8721516716  -1.1400717710
 \$END
EOF
# first \$VEC ... \$END block of the RHF .dat, renamed \$VEC1
awk 'BEGIN{p=0} /^ \$VEC/ && !done {p=1; print " $VEC1"; next}
     p && /^ \$END/ {print; p=0; done=1; next} p {print}' "$DAT" >> "$out"
echo "wrote $out"
