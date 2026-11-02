# Q50: N-only / S-only nested-sampling fits (dynesty, same likelihood as q39_q22.py/q42_dyn2.py) with a free own envelope and no cross-nucleus leakage.
# usage: q50_dyn.py <reg N|S> <mode ''|drop1> <seed> [nlive]  -> out/q50/runs/q50_<reg>_<mode>_s<seed>.npz/.json  (scenario cmb, env sh = single envelope of that nucleus, sigma_band 0.10)
import sys,os,json,time
reg,mode,seed=sys.argv[1],sys.argv[2],int(sys.argv[3]); nlive=int(sys.argv[4]) if len(sys.argv)>4 else 1500
os.environ['Q39_MODE']=mode; os.environ['Q22_CROSS']='0'
src=open('/workspace/scripts/q39_q22.py').read().split('def smc(')[0]
R=("REGS=('N','S')","REGS=(os.environ['Q50_REG'],)")
E1=("co,_=zones(p,o); fc=10**p['lfc_'+r]; fco=10**p['lfc_'+o]","co,_=zones(p,r if len(REGS)==1 else o); fc=10**p['lfc_'+r]; fco=10**p['lfc_'+(r if len(REGS)==1 else o)]")
E2=("(p['ne']>=p['nc_N'])|(p['ne']>=p['nc_S'])","(p['ne']>=p['nc_'+REGS[0]])")
for a,b in (R,E1,E2): assert a in src,a; src=src.replace(a,b)
os.environ['Q50_REG']=reg
sys.argv=['q39_q22.py','sh','q36cmb',str(seed)]; exec(compile(src,'q39_q22.py','exec'))
import dynesty
from dynesty import NestedSampler
def logl(x):
    c=chi2(x[None,:])[0]
    return -0.5*c if c<1e8 else -np.inf
def ptf(u): return LO+(HI-LO)*u
OUT='/workspace/out/q50/runs'; os.makedirs(OUT,exist_ok=True); t0=time.time()
ns=NestedSampler(logl,ptf,ND,nlive=nlive,sample='rslice',slices=int(os.environ.get('Q42_SL',8)),bound='multi',rstate=np.random.default_rng(seed))
ns.run_nested(dlogz=0.1,print_progress=False)
r=ns.results; w=np.exp(r.logwt-r.logz[-1]); w/=w.sum(); S=r.samples; tag=f'{reg}_{mode or "base"}_s{seed}'
np.savez(f'{OUT}/q50_{tag}.npz',S=S,w=w,names=names,logl=r.logl)
out=dict(reg=reg,mode=mode or 'base',seed=seed,nlive=nlive,logZ=float(r.logz[-1]),logZerr=float(r.logzerr[-1]),ess=float(w.sum()**2/(w**2).sum()),sec=time.time()-t0,chi2min=float(-2*r.logl.max()),names=list(names))
json.dump(out,open(f'{OUT}/q50_{tag}.json','w')); print(out)
