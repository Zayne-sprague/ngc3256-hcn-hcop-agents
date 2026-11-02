# Q12c: joint N-S propagation of the deblend ensemble: for each of the 16 template variants (same variant for N and S = draw-by-draw) fit N and S (1z, 2cT; settings base/geo; with and without H13CN).
# Output out/q12c/runs/var<k>_<tag>_... ; combined by q12c_summary.py into the joint Delta T posterior (mixture over variants).
import subprocess, os, itertools
from concurrent.futures import ThreadPoolExecutor
os.chdir('/workspace')
base=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',Q12_NP='5000')
SET=dict(base_={},geo=dict(Q12_GEOM=1,Q12_U=1,Q12_SK=.15))
jobs=list(itertools.product(range(16),('base_','geo'),('N','S'),('1z','2cT'),(0,1)))
def run(j):
    k,st,reg,m,h=j; ov=SET[st]
    if os.path.exists(f'out/q12c/runs/var{k}_{"h_" if h else ""}{st}_{m}_{reg}_r7_tbg2.73_s1.csv'): return
    env=dict(base,Q12_VAR=str(k),Q12_TAG=f'var{k}_{"h_" if h else ""}{st}',**{a:str(b) for a,b in ov.items()})
    if h: env['Q12_H13']='0.033' if reg=='N' else '0.05'
    r=subprocess.run(['python3','-I','scripts/q12c_twocomp.py',m,reg,'r7','g273','1'],env=env,capture_output=True,text=True)
    if r.returncode: open('out/q12c/errors_var.log','a').write(str(j)+r.stderr[-300:]+'\n')
with ThreadPoolExecutor(4) as ex: list(ex.map(run,jobs))
print(len(jobs),'variant jobs done')
