#!/bin/bash
# reduced secondary-variant set (time): XSTEP grid + XEXT, 10/15 and 20/25 priors
cd /workspace/scripts
for job in "N,S,C P10" "N,S P10C" "N,S PCS" "N Pw" "N Pp"; do
  set -- $job
  XSTEP=1 python3 -I s260b_q44_fit.py $1 $2 0,1
  XEXT=1 python3 -I s260b_q44_fit.py $1 $2 0,1
done
