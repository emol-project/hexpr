#!/bin/bash
# GAMESS inputs for OH X 2Pi, same geometry (Bohr) and basis as run_pyscf_oh.py.
#   ./make_inp.sh rohf            -> oh_rohf.inp
#   ./make_inp.sh mcscf           -> oh_mcscf.inp   ($VEC from oh_rohf.dat)
#   ./make_inp.sh so TRUE|3-21G   -> oh_so_<z>.inp  ($VEC from oh_mcscf.dat)
# The .dat files are read from GAMESS's USERSCR (~/gamess/restart).
set -e
R=${DATDIR:-$HOME/gamess/restart}
geom() {
cat <<EOF
 \$DATA
OH X2Pi cc-pVTZ, R = 0.9697 A
C1
O   8.0   0.0   0.0   0.0
H   1.0   0.0   0.0   1.8324674229907405
 \$END
EOF
}
vec() {   # first \$VEC ... \$END of a .dat; for MCSCF the last one (converged)
  awk -v which="$2" 'BEGIN{n=0} /^ \$VEC/{n++} {print n" "$0}' "$1" > /tmp/vec.$$
  last=$(awk '{print $1}' /tmp/vec.$$ | sort -n | tail -1)
  [ "$which" = first ] && last=1
  awk -v k="$last" 'BEGIN{p=0} $1==k && /\$VEC/{p=1; print " $VEC1"; next}
       p && /\$END/{sub(/^[0-9]+ /,""); print; p=0; next}
       p{sub(/^[0-9]+ /,""); print}' /tmp/vec.$$
  rm -f /tmp/vec.$$
}
case "$1" in
rohf)
  { echo ' $CONTRL SCFTYP=ROHF MULT=2 RUNTYP=ENERGY ISPHER=1 UNITS=BOHR $END'
    echo ' $SYSTEM MWORDS=100 $END'
    echo ' $BASIS GBASIS=CCT $END'
    echo ' $SCF CONV=1.0E-10 $END'
    geom; } > oh_rohf.inp ;;
mcscf)
  # state average over the two Pi components (equal weights)
  { echo ' $CONTRL SCFTYP=MCSCF MULT=2 RUNTYP=ENERGY ISPHER=1 UNITS=BOHR $END'
    echo ' $SYSTEM MWORDS=100 $END'
    echo ' $BASIS GBASIS=CCT $END'
    echo ' $GUESS GUESS=MOREAD NORB=44 $END'
    echo ' $DRT GROUP=C1 FORS=.TRUE. NMCC=1 NDOC=3 NALP=1 NVAL=1 $END'
    echo ' $GUGDIA NSTATE=4 CVGTOL=1.0E-10 $END'
    echo ' $GUGDM2 WSTATE(1)=1.0,1.0 $END'
    echo ' $MCSCF CISTEP=GUGA FULLNR=.TRUE. ACURCY=1.0E-8 MAXIT=100 $END'
    geom; vec "$R/oh_rohf.dat" first | sed 's/\$VEC1/$VEC/'; } > oh_mcscf.inp ;;
so)
  # ./make_inp.sh so ZEFTYP [OPERAT]: HSO1 runs are tagged by ZEFTYP
  # (oh_so_true, oh_so_3_21g), HSO2P/HSO2 by the operator (oh_so_hso2p, ...)
  Z=${2:-TRUE}; OP=${3:-HSO1}
  if [ "$OP" = HSO1 ]; then tag=$(echo "$Z" | tr 'A-Z-' 'a-z_')
  else tag=$(echo "$OP" | tr 'A-Z' 'a-z'); fi
  { echo ' $CONTRL SCFTYP=NONE CITYP=GUGA RUNTYP=TRANSITN MULT=2 ISPHER=1 UNITS=BOHR $END'
    echo ' $SYSTEM MWORDS=100 $END'
    echo ' $BASIS GBASIS=CCT $END'
    echo ' $GUGDIA CVGTOL=1.0E-10 $END'
    echo " \$TRANST OPERAT=$OP ZEFTYP=$Z NUMCI=1 NFZC=1 NOCC=6 SYMTOL=1.0E-10"
    echo '         IROOTS(1)=4 NSTATE(1)=4 $END'
    echo ' $DRT1 GROUP=C1 FORS=.TRUE. NFZC=1 NDOC=3 NALP=1 NVAL=1 $END'
    geom; vec "$R/oh_mcscf.dat" last; } > "oh_so_${tag}.inp" ;;
*) echo "usage: $0 rohf|mcscf|so [TRUE|3-21G] [HSO1|HSO2P|HSO2]"; exit 1 ;;
esac
