#!/bin/bash
# full re-runnable pipeline for Q12c (skips finished runs)
cd /workspace
python3 -I scripts/q12c_geometry.py
python3 -I scripts/q12c_variants.py
python3 -I scripts/run_q12c.py
python3 -I scripts/run_q12c_var.py
python3 -I scripts/run_q12c_calib.py 200
python3 -I scripts/q12c_summary.py > out/q12c/summary.log
python3 -I scripts/q12c_calib_summary.py > out/q12c/calib_summary.log
python3 -I scripts/run_q12c_calib2.py 200 > out/q12c/calib2.log 2>&1
python3 -I scripts/q12c_branches.py > out/q12c/branches.log
python3 -I scripts/q12c_calib_summary.py > out/q12c/calib_summary.log
