#!/bin/bash
# usage: q36_run.sh <env> <grid: ''|q36cmb|q36nominal|q36dusty|q36high> <seed> [envsuf]   (12000 particles, as q22b_runone.sh base run)
cd /workspace; Q36_ENVSUF=$4 Q22_TAG=${5:-q36} Q22_NP=12000 Q22_CROSS=0.5 Q22_SG=0.10 python3 scripts/q36_q22.py $1 "$2" $3 2>&1 | tail -1
