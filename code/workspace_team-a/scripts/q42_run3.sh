#!/bin/bash
# Q42 round 3: high-nlive (5000) repeats of base/drop1 (cmb), 2 seeds each -> out/q42/runs3 ; plus profile-likelihood checks
cd /workspace; export OMP_NUM_THREADS=1 Q42_SL=8 Q42_OUT=/workspace/out/q42/runs3
for m in base drop1; do for sd in 11 12; do nohup python3 scripts/q42_dyn2.py $m cmb $sd 5000 > out/q42/runs3_${m}_$sd.log 2>&1 & done; done
nohup python3 scripts/q42_prof.py "" cmb > out/q42/prof_base.log 2>&1 &
nohup python3 scripts/q42_prof.py drop1 cmb > out/q42/prof_drop1.log 2>&1 &
wait
