#!/bin/bash
# runs the Q12 model grid in parallel (14 procs); outputs out/q12/runs/*.csv ; then run scripts/q12_summary.py
cd /workspace; export Q12_NP=8000; mkdir -p out/q12
{
for suf in g273 _tbg5; do for reg in N S; do for seed in 1 2; do
 for ds in r8 r8a; do for m in 1z 2cT 2cI 2cTX; do echo "$m $reg $ds $suf $seed"; done; done
 for m in 2cT 2cI 2cTX; do echo "$m $reg full $suf $seed"; done
done; done; done
} | xargs -P 14 -L 1 bash -c 'python3 -I scripts/q12_twocomp.py $0 $1 $2 $3 $4 > /dev/null 2>>out/q12/errors.log'
