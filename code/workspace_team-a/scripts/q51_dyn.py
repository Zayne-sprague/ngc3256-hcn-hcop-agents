# Q51: nested-sampling (dynesty) refit of the Q12c two-zone (and one-zone) model; likelihood = scripts/q12c_twocomp.py verbatim (source cut before the SMC call).
# usage: q51_dyn.py <model 1z|2cT> <reg N|S> <geom base|geo> <h13 0|1> <data with|drop1>[_r7|_r7a] <seed> [nlive]   (r7 = 7 ratio data [default]; r7a = + compact intrinsic Tb of HCN/HCO+ 3-2, which fixes fc and hence the core mass)
#   geom 'base' = fixed geometry (Q12_GEOM=0), 'geo' = measured-geometry priors (Q12_GEOM=1,Q12_U=1,Q12_SK=.15 as run_q12c.py);  h13 1 = H13CN 4-3 limit term (alpha<0.033 N, <0.05 S)
#   drop1: masks the three 1.7" data that involve J=1-0 (R10, HCN 3-2/1-0, HCO+ 3-2/1-0; R32 already redundant/masked) -> 4 deblended 10"-ladder ratios only
# output: out/q51/runs/q51_<model>_<reg>_<geom>_h<h13>_<data>_s<seed>.npz/.json (weighted dynesty samples + names + logZ)
import sys,os,json,time
model,reg,geom,h13,dat,seed=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),sys.argv[5],int(sys.argv[6]); dat,dset=(dat.split('_')+['r7'])[:2]; nlive=int(sys.argv[7]) if len(sys.argv)>7 else 1000
if geom=='geo': os.environ.update(Q12_GEOM='1',Q12_U='1',Q12_SK='.15')
if h13: os.environ['Q12_H13']='0.033' if reg=='N' else '0.05'
src=open('/workspace/scripts/q12c_twocomp.py').read().split('NPART=int(')[0]
sys.argv=['q12c_twocomp.py',model,reg,dset,'g273',str(seed)]; exec(compile(src,'q12c_twocomp.py','exec'))
if dat=='drop1': MASK[5:8]=0
import dynesty
from dynesty import NestedSampler
def logl(x):
    c=chi2(x[None,:])[0]
    return -0.5*c if c<1e8 else -np.inf
def ptf(u): return LO+(HI-LO)*u
OUT='/workspace/out/q51/runs'; os.makedirs(OUT,exist_ok=True); t0=time.time()
ns=NestedSampler(logl,ptf,ND,nlive=nlive,sample='rslice',slices=int(os.environ.get('Q51_SL',6)),bound='multi',rstate=np.random.default_rng(seed))
ns.run_nested(dlogz=0.1,print_progress=False)
r=ns.results; w=np.exp(r.logwt-r.logz[-1]); w/=w.sum(); tag=f'{model}_{reg}_{geom}_h{h13}_{dat}_{dset}_s{seed}'
np.savez(f'{OUT}/q51_{tag}.npz',S=r.samples,w=w,names=names,logl=r.logl)
out=dict(model=model,reg=reg,geom=geom,h13=h13,dat=dat,dset=dset,seed=seed,nlive=nlive,logZ=float(r.logz[-1]+LOGZ_CORR),logZerr=float(r.logzerr[-1]),ess=float(w.sum()**2/(w**2).sum()),ncall=int(r.ncall.sum()),sec=time.time()-t0,chi2min=float(-2*r.logl.max()),names=list(names),ND=ND)
json.dump(out,open(f'{OUT}/q51_{tag}.json','w')); print(out)
