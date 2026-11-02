#!/bin/bash
# Q42 round 2: dynesty rslice, -inf for invalid points, 4 parallel jobs (cgroup limit 4 cpus). logs out/q42/runs2.log
cd /workspace; export OMP_NUM_THREADS=1 Q42_SL=8 Q42_OUT=/workspace/out/q42/runs2
{ for sd in 1 2 3; do for s in cmb nominal; do for m in base drop1 inf30; do echo "python3 scripts/q42_dyn2.py $m $s $sd 1500"; done; done; done
  echo "python3 scripts/q42_dyn2.py base cmb 1 1500 sep"; echo "python3 scripts/q42_dyn2.py drop1 cmb 1 1500 sep"; } | xargs -P 4 -I{} bash -c '{} 2>&1 | tail -1' >> out/q42/runs2.log
