# Q42: convergent nested-sampling (dynesty) rerun of the Q39 joint N+S fit. Reuses the q39_q22.py likelihood verbatim (source truncated before the SMC sampler).
# usage: q42_dyn.py <mode ''|drop1|inf30> <scen cmb|nominal> <seed> [nlive] [env sh|sep]   -> out/q42/runs/q42_<mode>_<scen>_<env>_s<seed>.npz/.json
import sys,os,json,time
mode,scen,seed=sys.argv[1],sys.argv[2],int(sys.argv[3]); nlive=int(sys.argv[4]) if len(sys.argv)>4 else 1500; envn=sys.argv[5] if len(sys.argv)>5 else 'sh'
os.environ['Q39_MODE']=mode; src=open('/workspace/scripts/q39_q22.py').read().split('def smc(')[0]
sys.argv=['q39_q22.py',envn,'q36'+scen,str(seed)]; exec(compile(src,'q39_q22.py','exec'))
import dynesty
from dynesty import NestedSampler
def logl(x):
    c=chi2(x[None,:])[0]
    return -0.5*c if c<1e8 else -np.inf
def ptf(u): return LO+(HI-LO)*u
SAMP=os.environ.get('Q42_SAMP','rslice'); SL=int(os.environ.get('Q42_SL',20)); WK=int(os.environ.get('Q42_WK',100)); OUT=os.environ.get('Q42_OUT','/workspace/out/q42/runs2'); os.makedirs(OUT,exist_ok=True)
t0=time.time()
ns=NestedSampler(logl,ptf,ND,nlive=nlive,sample=SAMP,slices=SL,walks=WK,bound='multi',rstate=np.random.default_rng(seed))
ns.run_nested(dlogz=0.1,print_progress=False,maxcall=int(os.environ.get("Q42_MAXCALL",10**9)),maxiter=None)
r=ns.results; w=np.exp(r.logwt-r.logz[-1]); w/=w.sum(); S=r.samples
ess=w.sum()**2/(w**2).sum()
np.savez(f'{OUT}/q42_{mode or "base"}_{scen}_{envn}_s{seed}.npz',S=S,w=w,names=names,logl=r.logl)
out=dict(mode=mode or 'base',scen=scen,env=envn,seed=seed,nlive=nlive,logZ=float(r.logz[-1]),logZerr=float(r.logzerr[-1]),ess=float(ess),ncall=int(r.ncall.sum()),sec=time.time()-t0,chi2min=float(-2*r.logl.max()))
json.dump(out,open(f'{OUT}/q42_{mode or "base"}_{scen}_{envn}_s{seed}.json','w')); print(out)
