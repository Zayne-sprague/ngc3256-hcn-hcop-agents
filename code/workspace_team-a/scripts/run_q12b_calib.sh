#!/bin/bash
# Bayes-factor / chi2 calibration: simulate the 8 ratio data from the best-fit one-zone (or two-comp) N model with the observed errors, refit with 1z and 2cT; 40 sims each.
cd /workspace; export Q12_NP=4000
python3 -I - <<'PY'
import numpy as np
for t in ('1z','2cT'):
    z=np.load(f'/workspace/out/q12b/runs/base_{t}_N_r8_tbg2.73_s1.npz'); np.savez(f'/workspace/out/q12b/truth_{t}.npz',qbest=z['qbest'])
PY
{ for t in 1z 2cT; do for s in $(seq 1 40); do for m in 1z 2cT; do echo "$m N r8 g273 1 /workspace/out/q12b/truth_$t.npz:$s"; done; done; done; } | \
 xargs -P 15 -L 1 bash -c 'Q12_TAG=cal Q12_SIM=$5 python3 -I scripts/q12b_twocomp.py $0 $1 $2 $3 $4 >/dev/null 2>>out/q12b/errors_cal.log'
