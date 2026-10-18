#!/bin/bash
# Run every GAMESS step of the OH comparison, in order.
set -e
cd "$(dirname "$0")"
G=${RUNGMS:-$HOME/gamess/rungms}
bash make_inp.sh rohf  && $G oh_rohf  00 1 > oh_rohf.log  2>&1
bash make_inp.sh mcscf && $G oh_mcscf 00 1 > oh_mcscf.log 2>&1
for z in TRUE 3-21G; do
  t=$(echo "$z" | tr 'A-Z-' 'a-z_')
  bash make_inp.sh so "$z" && $G "oh_so_$t" 00 1 > "oh_so_$t.log" 2>&1
done
for op in HSO2P HSO2; do
  t=$(echo "$op" | tr 'A-Z' 'a-z')
  bash make_inp.sh so TRUE "$op" && $G "oh_so_$t" 00 1 > "oh_so_$t.log" 2>&1
done
for f in oh_rohf oh_mcscf oh_so_true oh_so_3_21g oh_so_hso2p oh_so_hso2; do
  printf '%-12s ' "$f"; grep -c 'TERMINATED NORMALLY' "$f.log"
done
