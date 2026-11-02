#!/bin/bash
# Q50: 8 single-nucleus NS fits (N,S x base,drop1 x seeds 1,2); ~10-15 min each in parallel
cd /workspace; export OMP_NUM_THREADS=1
for s in 1 2; do for r in N S; do for m in "" drop1; do python3 scripts/q50_dyn.py $r "$m" $s 1500 > out/q50/log_${r}_${m:-base}_$s.txt 2>&1 & done; done; done; wait
