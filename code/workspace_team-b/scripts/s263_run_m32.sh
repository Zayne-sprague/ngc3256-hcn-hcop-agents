#!/bin/bash
# Q44 round 2: M3 (3-2 12C only; X-independent -> one X) and M32 (3-2 12C + 13C 1-0) for N,S,C; band priors k=0,1 ; then k=2 (flat 40/40%) for N
cd /workspace/scripts
LXONLY=1 python3 -I s261_q44_matched_fit.py N,S,C M3 0,1,2
export XSTEP=1
python3 -I s261_q44_matched_fit.py N,S,C M32 0,1
XEXT=1 python3 -I s261_q44_matched_fit.py N,S,C M32 0,1
python3 -I s261_q44_matched_fit.py N M32 2
XEXT=1 python3 -I s261_q44_matched_fit.py N M32 2
