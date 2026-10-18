#!/bin/bash
# hexpr side of the two-electron comparison (after gamess/run_all.sh).
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python}
mkdir -p results
$PY run_pyscf_oh.py --somf core --gamess gamess/oh_so_hso2p.log > results/hexpr_somf_core.txt
$PY run_pyscf_oh.py --somf full --gamess gamess/oh_so_hso2.log  > results/hexpr_somf_full.txt
grep -h -E "SOMF density|vs determinants|^A =|^GAMESS" results/hexpr_somf_core.txt results/hexpr_somf_full.txt
