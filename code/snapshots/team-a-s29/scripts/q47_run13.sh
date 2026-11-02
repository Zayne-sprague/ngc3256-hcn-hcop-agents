#!/bin/bash
# Q47: N-only NS with/without the 1-0 13C/12C ratios in the likelihood; base: data A,T,none x seeds 1-3; drop1(1-0 12C dropped): A,none x seeds 1-2
cd /workspace; export OMP_NUM_THREADS=1; mkdir -p out/q47/runs
for s in 1 2 3; do for d in A T none; do python3 scripts/q47_dyn13.py N "" $s $d 1500 > out/q47/runs/log_N_base_${d}_$s.txt 2>&1 & done; done
for s in 1 2; do for d in A none; do python3 scripts/q47_dyn13.py N drop1 $s $d 1500 > out/q47/runs/log_N_drop1_${d}_$s.txt 2>&1 & done; done
wait
