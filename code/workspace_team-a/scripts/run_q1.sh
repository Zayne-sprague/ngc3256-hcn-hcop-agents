#!/bin/sh
# Q1 pipeline (run from anywhere); outputs in /workspace/out/q1
cd /workspace && for s in q1_ratios q1_maps q1_wings q1_disk q1_fluxrec; do python3 -I scripts/$s.py > out/q1/$s.log 2>&1 || echo FAIL $s; done
python3 -I scripts/build_report_q1.py
# revision 2 (deblending, velocity ratios, disk recovery model, RADEX, report)
for s in q1_deblend q1_velratio q1_redwing_check q1_disk_model q1_radex; do python3 -I scripts/$s.py > out/q1/$s.log 2>&1 || echo FAIL $s; done
python3 -I scripts/build_report_q1b.py
