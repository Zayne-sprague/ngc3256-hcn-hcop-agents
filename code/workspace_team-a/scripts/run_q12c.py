# Q12c run matrix driver (7 independent data r7*, measured geometry priors 'geo*', H13CN term 'h_*'). Outputs out/q12c/runs/*.csv/npz; then run q12c_summary.py
import subprocess, os, itertools, sys
from concurrent.futures import ThreadPoolExecutor
os.chdir('/workspace'); os.makedirs('out/q12c',exist_ok=True)
V=dict(base={},us=dict(Q12_U=.5,Q12_SK=.15),us30=dict(Q12_U=.5,Q12_SK=.30),floor25=dict(Q12_FLOOR=.25),
 geo0=dict(Q12_GEOM=1,Q12_U=1),geo=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15),geo30=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.30),geo_floor25=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15,Q12_FLOOR=.25),
 nc6=dict(Q12_NCMAX=6),nc7=dict(Q12_NCMAX=7),geo_nc7=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15,Q12_NCMAX=7),geo_nc6=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15,Q12_NCMAX=6))
jobs=[]
for reg,seed,h in itertools.product('NS',(1,2,3),(0,1)):
    for m in ('1z','2cT','2cI'):
        for t in V: jobs.append((t,m,reg,'r7','g273',seed,h))
        for t in ('base','geo'):
            for ds in ('r7a','r7full'): jobs.append((t,m,reg,ds,'g273',seed,h))
            if h==0: jobs.append((t,m,reg,'r7','_tbg5',seed,0))
    for t in ('base','geo'): jobs.append((t,'2cTX',reg,'r7','g273',seed,h))
def run(j):
    t,m,reg,ds,suf,seed,h=j
    tagn=('h_' if h else '')+t; fn=f"out/q12c/runs/{tagn}_{m}_{reg}_{ds}_{'tbg2.73' if suf=='g273' else suf}_s{seed}.csv"
    if os.path.exists(fn): return
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',Q12_NP='8000',Q12_TAG=('h_' if h else '')+t,**{k:str(v) for k,v in V[t].items()})
    if h: env['Q12_H13']='0.033' if reg=='N' else '0.05'
    r=subprocess.run(['python3','-I','scripts/q12c_twocomp.py',m,reg,ds,suf,str(seed)],env=env,capture_output=True,text=True)
    if r.returncode: open('out/q12c/errors.log','a').write(str(j)+r.stderr[-500:]+'\n')
with ThreadPoolExecutor(4) as ex: list(ex.map(run,jobs))
print(len(jobs),'jobs done')
