#!/bin/bash
# all Q2 sensitivity configurations (floor test at Tbg 2.73 K; Tbg test with floors 0 and 30 K)
cd /workspace
for t in 0 15 20 30; do python3 scripts/q2_sens.py none $t; done
for s in _tbg5 _tbg10 _tbg20; do for t in 0 30; do python3 scripts/q2_sens.py $s $t; done; done
