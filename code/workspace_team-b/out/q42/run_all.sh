#!/bin/bash
cd /workspace
for spec in "N H4c 0,1,2" "N H4g 0,1" "N H4d 0,1" "N T 0,1" "N T32 0,1" "N T432 0,1" "S T 0,1" "C T 0,1" "S T432 0,1" "C T432 0,1"; do
  set -- $spec
  python3 -I scripts/s250b_q42_ratio_fit.py $1 $2 $3
done
