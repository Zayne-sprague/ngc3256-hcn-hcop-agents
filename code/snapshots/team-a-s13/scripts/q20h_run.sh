#!/bin/bash
# Q20h: (i) size/Tint scaling TSC (1.5 = Q5 Gaussian core, 2.5 ~ s0.14", 4 ~ s0.10" floor), (ii) T floor 50 K (sub-thermal single phase with warm Tk): n range needed.
cd /workspace
for reg in S N; do
 for tsc in 1.5 2.5 4; do for s in 1 2; do Q20_LQ=0 Q20_TSC=$tsc Q12_TAG=tsc$tsc python3 scripts/q20b_eta.py 1z $reg r8a "" $s 2>&1 | tail -1 & done; done
 for s in 1 2; do Q20_LQ=0 Q20_TMIN=50 Q12_TAG=tmin50 python3 scripts/q20b_eta.py 1z $reg r8a "" $s 2>&1 | tail -1 &
                  Q20_LQ=0 Q20_TMIN=50 Q12_TAG=tmin50r8 python3 scripts/q20b_eta.py 1z $reg r8 "" $s 2>&1 | tail -1 & done
 wait
done
