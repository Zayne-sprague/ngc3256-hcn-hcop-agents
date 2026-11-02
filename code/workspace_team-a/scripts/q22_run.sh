#!/bin/bash
# Q22 battery: models x regions x data x Tbg x seeds (baseline), then sensitivity variants (core size, free envelope size at 10", 10" error scale)
cd /workspace; mkdir -p out/q22/runs
run(){ python3 scripts/q22_joint.py "$@" 2>&1 | tail -1; }
for reg in S N; do for ds in a ab; do for suf in '' _tbg5; do for m in M1 M1E M2T M2 M2H; do for s in 1 2; do
  Q12_TAG=base run $m $reg $ds "$suf" $s &
  while [ $(jobs -r | wc -l) -ge 6 ]; do sleep 1; done
done; done; done; done; done; wait
for reg in S N; do for m in M1E M2 M2H; do for s in 1 2; do
  Q12_TAG=tsc1.5 Q20_TSC=1.5 run $m $reg ab "" $s &
  Q12_TAG=lr Q22_LR=1 run $m $reg ab "" $s &
  Q12_TAG=sg15 Q22_SG=0.15 run $m $reg ab "" $s &
  Q12_TAG=aps2 Q22_APS=2 run $m $reg ab "" $s &
  while [ $(jobs -r | wc -l) -ge 6 ]; do sleep 1; done
done; done; done; wait
