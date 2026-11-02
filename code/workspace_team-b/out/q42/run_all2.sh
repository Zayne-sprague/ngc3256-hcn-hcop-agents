#!/bin/bash
cd /workspace
while pgrep -f run_all.sh >/dev/null; do sleep 30; done
for spec in "N Tp 0,1" "S Tp 0,1" "C Tp 0,1" "N T150 0,1" "S T150 0,1" "C T150 0,1"; do set -- $spec; python3 -I scripts/s250b_q42_ratio_fit.py $1 $2 $3; done
