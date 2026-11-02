# Q12c calibration 2: mis-specified truths that STILL FIT THE OBSERVED N DATA. The truth model is first refit to the real N data with the geometry overridden
# (M1: 1z with 2-1 compact size x1.35 and 4-3 x0.74; M2: 1z with 2-1 x0.74, 4-3 x1.35; M3: 2cT with u=2 [extended filling at 1.7" x2]; M4: 2cT with u=0.5), its best-fit predictions
# are the truth; simulated data (obs. errors) are then fitted with the geometry-naive models ('base': sizes common, u=1) and with the measured-geometry-prior setting ('geo').
# Output out/q12c/runs/cal_<truth>_<fit>_*; summarised by q12c_calib_summary.py
import subprocess, os, sys, numpy as np
from concurrent.futures import ThreadPoolExecutor
os.chdir('/workspace'); NS=int(sys.argv[1]) if len(sys.argv)>1 else 200; LN=np.log
base=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
TR=dict(M1=('1z',dict(Q12_OVK2=LN(1.35),Q12_OVK4=LN(.74))),M2=('1z',dict(Q12_OVK2=LN(.74),Q12_OVK4=LN(1.35))),M3=('2cT',dict(Q12_OVU=LN(2.))),M4=('2cT',dict(Q12_OVU=LN(.5))))
def sh(args,env): return subprocess.run(['python3','-I','scripts/q12c_twocomp.py']+args,env=env,capture_output=True,text=True)
for t,(m,ov) in TR.items():
    e=dict(base,Q12_NP='8000',Q12_TAG='ov'+t,**{k:str(v) for k,v in ov.items()}); r=sh([m,'N','r7','g273','1'],e); print(t,'truth fit to real N data:',r.stdout.strip()[:140])
    e=dict(base,Q12_DUMPQ=f'out/q12c/calib/truth_{t}.npz',Q12_BESTFROM=f'out/q12c/runs/ov{t}_{m}_N_r7_tbg2.73_s1.csv',**{k:str(v) for k,v in ov.items()}); print(sh([m,'N','r7','g273','1'],e).stdout.strip())
SET=dict(base={},geo=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15))
jobs=[(t,st,s,m) for t in TR for st in ('base','geo') for s in range(1,NS+1) for m in ('1z','2cT')]
def run(j):
    t,st,s,m=j
    if os.path.exists(f'out/q12c/runs/cal_{t}_{st}_{m}_N_r7_tbg2.73_s1_sim{s}_truth_{t}.csv'): return
    env=dict(base,Q12_NP='3000',Q12_TAG=f'cal_{t}_{st}',Q12_SIM=f'out/q12c/calib/truth_{t}.npz:{s}',**{k:str(v) for k,v in SET[st].items()})
    r=sh([m,'N','r7','g273','1'],env)
    if r.returncode: open('out/q12c/errors_cal.log','a').write(str(j)+r.stderr[-300:]+'\n')
with ThreadPoolExecutor(4) as ex: list(ex.map(run,jobs))
print(len(jobs),'calibration-2 jobs done')
