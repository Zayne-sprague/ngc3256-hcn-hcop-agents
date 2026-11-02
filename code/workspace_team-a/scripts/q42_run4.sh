#!/bin/bash
# Q42 round 4: 10" band-scale prior sigma variants (Q22_SG=0.15,0.25; baseline 0.10), cmb, nlive 1500 seed 1 -> out/q42/runs4_sg<sg>
cd /workspace; export OMP_NUM_THREADS=1 Q42_SL=8
for sg in 0.15 0.25; do for m in base drop1; do Q22_SG=$sg Q42_OUT=/workspace/out/q42/runs4_sg$sg nohup python3 scripts/q42_dyn2.py $m cmb 1 1500 > out/q42/runs4_${m}_$sg.log 2>&1 & done; done; wait
