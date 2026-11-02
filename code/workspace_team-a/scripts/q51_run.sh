#!/bin/bash
# Q51 run matrix: dset{r7a,r7} x model{2cT,1z} x reg{N,S} x geom{geo,base} x h13{0,1} x data{with,drop1} x seeds 1-3, nlive 1000; 4 CPUs in the cgroup -> 8 parallel
cd /workspace; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; mkdir -p out/q51/runs out/q51/logs
for s in 1 2 3; do for ds in r7a r7; do for m in 2cT 1z; do for r in N S; do for g in geo base; do for h in 0 1; do for d in with drop1; do
 [ -f out/q51/runs/q51_${m}_${r}_${g}_h${h}_${d}_${ds}_s${s}.json ] || echo "$m $r $g $h ${d}_${ds} $s $ds"
done; done; done; done; done; done; done | xargs -P 8 -L 1 bash -c 'python3 scripts/q51_dyn.py $0 $1 $2 $3 $4 $5 1000 > out/q51/logs/$0_$1_$2_$3_$4_$5_$6.txt 2>&1'
