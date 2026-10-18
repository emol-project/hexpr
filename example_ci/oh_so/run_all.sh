#!/bin/bash
# hexpr side of the OH comparison; run gamess/run_all.sh first.
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python}
mkdir -p results
$PY run_pyscf_oh.py --gamess gamess/oh_so_true.log > results/hexpr_true.txt
$PY run_pyscf_oh.py --zeff O=6.0,H=1.0 --gamess gamess/oh_so_3_21g.log > results/hexpr_zeff_3_21g.txt
grep -h -E "^A =|^GAMESS|vs determinants" results/*.txt
