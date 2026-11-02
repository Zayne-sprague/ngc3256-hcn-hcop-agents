#!/bin/bash
# Q22 CO-bound variants: Tk floor for both zones from the Q32/p0040 CO(2-1) ff=1 limits (S: 30 window-mean, 40 peak-native; N: 20) -> Q20_TMIN
cd /workspace
run(){ python3 scripts/q22_joint.py "$@" 2>&1 | tail -1; }
for cfg in "S 30" "S 40" "N 20"; do set -- $cfg; reg=$1; tm=$2
 for m in M1E M2T M2; do for s in 1 2; do Q12_TAG=tmin$tm Q20_TMIN=$tm run $m $reg ab "" $s & done; done
 wait
done
