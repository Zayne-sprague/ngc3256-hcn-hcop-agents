# Q12c calibration of dlnZ(2cT-1z) and chi2 under truths: (T1) one-zone best fit; (T2) two-comp (2cT) best fit; (T3) MIS-SPECIFIED one-zone truth with different geometry (2-1 size x1.35, 4-3 size x0.74);
# (T4) MIS-SPECIFIED two-comp truth with u (extended filling at 1.7") x2; (T3b) one-zone truth with opposite geometry shift (2-1 x0.74, 4-3 x1.35).
# Sim data = truth predictions + Gaussian noise with the observed errors (7 independent data r7); refit with 1z and 2cT (setting 'base' = no geometry nuisance; 'geo' = measured-geometry priors u,3-2 size + 15% 2-1/4-3 size).
# NS sims per truth/setting (default 200). Output out/q12c/calib/*.csv ; summarised by q12c_calib_summary.py
import subprocess, os, sys
from concurrent.futures import ThreadPoolExecutor
os.chdir('/workspace'); os.makedirs('out/q12c/calib',exist_ok=True); NS=int(sys.argv[1]) if len(sys.argv)>1 else 200
R='out/q12c/runs/base_%s_N_r7_tbg2.73_s1.csv'
base=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
def dump(name,model,**ov):
    env=dict(base,Q12_DUMPQ=f'out/q12c/calib/truth_{name}.npz',Q12_BESTFROM=R%model,**{k:str(v) for k,v in ov.items()})
    r=subprocess.run(['python3','-I','scripts/q12c_twocomp.py',model,'N','r7','g273','1'],env=env,capture_output=True,text=True); print(name,r.stdout.strip(),r.stderr[-300:])
import numpy as np
LN=np.log
dump('T1','1z'); dump('T2','2cT'); dump('T3','1z',Q12_OVK2=LN(1.35),Q12_OVK4=LN(0.74)); dump('T3b','1z',Q12_OVK2=LN(0.74),Q12_OVK4=LN(1.35)); dump('T4','2cT',Q12_OVU=LN(2.0))
SET=dict(base={},geo=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15))
jobs=[]
for t in ('T1','T3','T3b','T2','T4'):
    for st in (('base','geo') if t in ('T1','T3','T3b') else ('base',)):
        for s in range(1,NS+1):
            for m in ('1z','2cT'): jobs.append((t,st,s,m))
def run(j):
    t,st,s,m=j
    if os.path.exists(f'out/q12c/runs/cal_{t}_{st}_{m}_N_r7_tbg2.73_s1_sim{s}_truth_{t}.csv'): return
    env=dict(base,Q12_NP='3000',Q12_TAG=f'cal_{t}_{st}',Q12_SIM=f'out/q12c/calib/truth_{t}.npz:{s}',**{k:str(v) for k,v in SET[st].items()})
    r=subprocess.run(['python3','-I','scripts/q12c_twocomp.py',m,'N','r7','g273','1'],env=env,capture_output=True,text=True)
    if r.returncode: open('out/q12c/errors_cal.log','a').write(str(j)+r.stderr[-300:]+'\n')
with ThreadPoolExecutor(4) as ex: list(ex.map(run,jobs))
print(len(jobs),'calibration jobs done')
