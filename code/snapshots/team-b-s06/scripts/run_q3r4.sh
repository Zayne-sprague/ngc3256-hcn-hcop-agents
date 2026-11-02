#!/bin/bash
# Q3 round 4 pipeline (python3 -I; ~6 min). Order matters: grid -> beam dilution -> centre fit -> disk contrast -> 1-0/derived -> disk positions -> report.
cd /workspace/scripts
python3 -I s62_radex_grid.py && python3 -I s63_beam_dilution.py && python3 -I s64_centre_fit.py && python3 -I s65_disk_contrast.py && python3 -I s66_j10_derived.py && python3 -I s67_disk_positions.py && python3 -I s68_update_report_q3r4.py
