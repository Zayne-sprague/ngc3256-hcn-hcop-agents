#!/bin/bash
# Q20: eta variants. LQ=0: same eta for HCN & HCO+; LQ=1: eta_HCO+/eta_HCN free within x10.
cd /workspace
for reg in N S; do for ds in r8 r8a; do for lq in 0 1; do for suf in '' _tbg5; do for m in 1z 2cT; do for s in 1 2; do
 Q20_LQ=$lq Q12_TAG=lq$lq python3 scripts/q20b_eta.py $m $reg $ds "$suf" $s 2>&1 | tail -1 &
done; done; wait; done; done; done; done
