#!/bin/bash
cd /workspace/scripts
for v in M; do for ext in 0 1; do
  if [ $ext = 1 ]; then export XEXT=1; else unset XEXT; fi
  python3 -I s261_q44_matched_fit.py N,S,C $v 0,1; done; done
unset XEXT; export XSTEP=1
for v in M10 MC M12 M10C M32; do python3 -I s261_q44_matched_fit.py N,S,C $v 0,1; XEXT=1 python3 -I s261_q44_matched_fit.py N,S,C $v 0,1; done
python3 -I s261_q44_matched_fit.py N M4 0,1; XEXT=1 python3 -I s261_q44_matched_fit.py N M4 0,1
