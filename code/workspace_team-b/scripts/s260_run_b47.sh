#!/bin/bash
cd /workspace/scripts
for v in B6 B6n B22 B22n; do python3 -I s260d_q44_b47_fit.py N,S,C $v 0 1 2; done
