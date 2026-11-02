#!/bin/bash
# 16 parallel dynesty fits; logs out/q42/runs.log
cd /workspace; export OMP_NUM_THREADS=1
{ for m in base drop1 inf30; do for s in cmb nominal; do for sd in 1 2; do echo "python3 scripts/q42_dyn.py $m $s $sd 2000"; done; done; done
  echo "python3 scripts/q42_dyn.py base cmb 3 4000"; echo "python3 scripts/q42_dyn.py drop1 cmb 3 4000"
  echo "python3 scripts/q42_dyn.py base cmb 1 2000 sep"; echo "python3 scripts/q42_dyn.py drop1 cmb 1 2000 sep"; } | xargs -P 16 -I{} bash -c '{} 2>&1 | tail -1' > out/q42/runs.log
