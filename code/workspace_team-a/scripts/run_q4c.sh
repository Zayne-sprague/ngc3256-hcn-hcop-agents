#!/bin/sh
# re-run Q4c (needs /workspace/lit/harada/harada.txt = text of arxiv.org/html/1801.05941, lit/other/mich.txt)
set -e; cd /workspace
for s in q4c_lit q4c_lit_ratios q4c_grid q4c_joint q4c_sens q4c_ns q4c_phys q4c_cs q4c_h13cn q4c_so_scope; do python3 -I scripts/$s.py > out/q4c/$s.log 2>&1; done
python3 scripts/build_report_q4c.py
