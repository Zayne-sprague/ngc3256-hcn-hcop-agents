#!/bin/bash
cd /workspace/scripts
python3 -I s260b_q44_fit.py N,S,C P 0,1,2
XEXT=1 python3 -I s260b_q44_fit.py N,S,C P 0,1,2
for v in P10 PCS P10C Pw Pp; do
  XSTEP=1 python3 -I s260b_q44_fit.py N,S,C $v 0,1
  XEXT=1 python3 -I s260b_q44_fit.py N,S,C $v 0,1
done
