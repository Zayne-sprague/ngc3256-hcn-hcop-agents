#!/bin/bash
# Q12b run matrix (variants via env vars, see header of q12b_twocomp.py). Outputs out/q12b/runs/. Then run q12b_summary.py
cd /workspace; export Q12_NP=8000; mkdir -p out/q12b
export V_base="" V_u="Q12_U=0.5" V_us="Q12_U=0.5 Q12_SK=0.15" V_us30="Q12_U=0.5 Q12_SK=0.30"
#V_u="Q12_U=0.5"; V_us="Q12_U=0.5 Q12_SK=0.15"; V_us30="Q12_U=0.5 Q12_SK=0.30"
export V_nc6="Q12_NCMAX=6" V_nc7="Q12_NCMAX=7" V_flo2="Q12_FLO=-2" V_floor25="Q12_FLOOR=0.25" V_hsig50="Q12_HSIG=0.5"
export V_us_nc7="Q12_U=0.5 Q12_SK=0.15 Q12_NCMAX=7" V_us_nc6="Q12_U=0.5 Q12_SK=0.15 Q12_NCMAX=6"
{
for reg in N S; do for seed in 1 2; do
 for v in base u us us30 nc6 nc7 flo2 floor25 hsig50 us_nc7 us_nc6; do for m in 1z 2cT 2cI; do echo "$v|$m $reg r8 g273 $seed"; done; done
 for v in base us; do for m in 2cTX; do echo "$v|$m $reg r8 g273 $seed"; done
   for ds in r8a full; do for m in 1z 2cT 2cI; do echo "$v|$m $reg $ds g273 $seed"; done; done
   for m in 1z 2cT 2cI; do echo "$v|$m $reg r8 _tbg5 $seed"; done; done
done; done
} | grep -v '^u|1z' | xargs -P 15 -I{} bash -c 'a="{}"; v=${a%%|*}; rest=${a#*|}; eval vv=\$V_$v; env $vv Q12_TAG=$v python3 -I scripts/q12b_twocomp.py $rest > /dev/null 2>>out/q12b/errors.log'
