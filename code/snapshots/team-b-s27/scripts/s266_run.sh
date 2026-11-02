#!/bin/bash
cd /workspace/scripts; export XSTEP=1
for v in M M32; do python3 -I s266_q44_core_ring_fit.py Nc,Nr,Sc,Sr,Cc,Cr $v 0,1; XEXT=1 python3 -I s266_q44_core_ring_fit.py Nc,Nr,Sc,Sr,Cc,Cr $v 0,1; done
