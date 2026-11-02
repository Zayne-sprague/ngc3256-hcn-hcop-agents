#!/bin/sh
# Q12 pipeline (about 6 min on 16 cores). python3 -I; scripts add /workspace/scripts to sys.path themselves.
cd /workspace/scripts
for o in Xb4f Xb4d Xb4b; do for l in HCO+ HCN; do python3 -I s56a_maps.py $o $l & done; done; wait
python3 -I s56a_maps.py Xb4f HCO+ 31; python3 -I s56a_maps.py Xb4f HCN 31     # wide crop (bowl check)
python3 -I s56b_measure.py
python3 -I s56c_transfer_fit.py      # negative result (chi2_red 25-72), ~2 min
python3 -I s56d_recovery_model.py
python3 -I s56e_extrapolate.py
python3 -I s57_update_report_q12.py
# ---- round 2 (reviewer fixes): empirical ACA filtering, clump decontamination, sensitivity, forward model
python3 -I s58a_fullmaps.py; python3 -I s58b_bowl.py; python3 -I s58c_contrecov.py; python3 -I s58d_clumpfit.py; python3 -I s58e_sensitivity.py; python3 -I s58f_simulate.py
python3 -I s59_update_report_q12r2.py
