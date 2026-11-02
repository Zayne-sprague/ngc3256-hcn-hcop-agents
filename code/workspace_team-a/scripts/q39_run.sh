#!/bin/bash
# runs all q39 fits (16 seeds) in parallel; results out/q39/runs/*.csv|npz ; logs out/q39/runs.log
cd /workspace; export OMP_NUM_THREADS=1; xargs -P 14 -I{} bash -c '{} 2>&1 | tail -1' < scripts/q39_jobs.txt > out/q39/runs.log
