#!/bin/bash
# time-limited remainder: CS variants (3-species fits are ~4x slower) for N only on the 0.1-dex X grid; Pw, Pp for N with XSTEP+XEXT
cd /workspace/scripts
XSTEP=1 python3 -I s260b_q44_fit.py N PCS 0,1
XSTEP=1 python3 -I s260b_q44_fit.py N P10C 0,1
XSTEP=1 python3 -I s260b_q44_fit.py N Pw 0,1
XEXT=1 python3 -I s260b_q44_fit.py N Pw 0,1
XSTEP=1 python3 -I s260b_q44_fit.py N Pp 0,1
XEXT=1 python3 -I s260b_q44_fit.py N Pp 0,1
