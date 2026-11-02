# Q47 (fix 1): q50_dyn.py with the 2.2" 1-0 13C/12C ratios (H13CN/HCN, H13CO+/HCO+, out/q47/ratios10.json) INCLUDED IN THE NESTED-SAMPLING LIKELIHOOD (X=12C/13C marginalised on a 12-point
# quantile grid of the prior LN(45,0.2 dex)<=100, 10% model error); replaces importance reweighting (Zfac 0.2-1.5%).  usage: q47_dyn13.py <reg> <mode> <seed> <A|T|none> [nlive]
# Q50: N-only / S-only nested-sampling fits (dynesty, same likelihood as q39_q22.py/q42_dyn2.py) with a free own envelope and no cross-nucleus leakage.
# usage: q50_dyn.py <reg N|S> <mode ''|drop1> <seed> [nlive]  -> out/q50/runs/q50_<reg>_<mode>_s<seed>.npz/.json  (scenario cmb, env sh = single envelope of that nucleus, sigma_band 0.10)
import sys,os,json,time
reg,mode,seed=sys.argv[1],sys.argv[2],int(sys.argv[3]); DAT=sys.argv[4]; nlive=int(sys.argv[5]) if len(sys.argv)>5 else 1500
os.environ['Q39_MODE']=mode; os.environ['Q22_CROSS']='0'
src=open('/workspace/scripts/q39_q22.py').read().split('def smc(')[0]
R=("REGS=('N','S')","REGS=(os.environ['Q50_REG'],)")
E1=("co,_=zones(p,o); fc=10**p['lfc_'+r]; fco=10**p['lfc_'+o]","co,_=zones(p,r if len(REGS)==1 else o); fc=10**p['lfc_'+r]; fco=10**p['lfc_'+(r if len(REGS)==1 else o)]")
E2=("(p['ne']>=p['nc_N'])|(p['ne']>=p['nc_S'])","(p['ne']>=p['nc_'+REGS[0]])")
for a,b in (R,E1,E2): assert a in src,a; src=src.replace(a,b)
os.environ['Q50_REG']=reg
SAVEARGV=sys.argv; sys.argv=['q39_q22.py','sh','q36cmb',str(seed)]; exec(compile(src,'q39_q22.py','exec'))
import dynesty
from dynesty import NestedSampler

from scipy.ndimage import map_coordinates
from scipy.stats import norm
G12={m:np.load(f'/workspace/out/q36/grid_{m}_{reg}_cmb.npz') for m in('hcn','hco+')}; G13=np.load('/workspace/out/q47/grid_h13cn_10.npz')
LTK=np.log10(G12['hcn']['TK']); LNg=G12['hcn']['LN']; LNC=G12['hcn']['LNC']; LNC13=G13['LNC']
def lg(a):
    a=np.array(a,float); a[~((a>0)&(a<1e3))]=np.nan; return np.log10(a)
L12={'hcn':lg(G12['hcn']['TR'][...,0]),'hco+':lg(G12['hco+']['TR'][...,0])}; L13=lg(G13['TR'])
def ip(arr,lT,n,lN,axN):
    c=np.array([np.interp(lT,LTK,np.arange(len(LTK)))*np.ones_like(lN),np.interp(n,LNg,np.arange(len(LNg)))*np.ones_like(lN),np.interp(lN,axN,np.arange(len(axN)))])
    return 10**map_coordinates(arr,c,order=1,mode='nearest')
XG=np.exp(np.log(45)+0.2*np.log(10)*norm.ppf((np.arange(12)+0.5)/12)); XG=XG[XG<=100]; LX=np.log10(XG)
RJ=json.load(open('/workspace/out/q47/ratios10.json')); DD=None if DAT=='none' else {'hcn':RJ[f'HCN|{reg}|{DAT}'],'hcop':RJ[f'HCO+|{reg}|{DAT}']}
def lnL13(x):
    p=unpack(x[None,:]); g=lambda k: float(np.asarray(p[k]).ravel()[0])
    lT,n,NH,NO=g('lTc_'+reg),g('nc_'+reg),g('NHc_'+reg),g('NOc_'+reg)
    h12=ip(L12['hcn'],lT,n,np.array([NH]),LNC)[0]; o12=ip(L12['hco+'],lT,n,np.array([NO]),LNC)[0]
    rh=ip(L13,lT,n,NH-LX,LNC13)/h12; ro=ip(L12['hco+'],lT,n,NO-LX,LNC)/o12
    ll=-0.5*((rh-DD['hcn']['R'])/np.hypot(DD['hcn']['err'],0.10*rh))**2-0.5*((ro-DD['hcop']['R'])/np.hypot(DD['hcop']['err'],0.10*ro))**2
    ll=np.nan_to_num(ll,nan=-1e3); m=ll.max(); return m+np.log(np.mean(np.exp(ll-m)))
def logl(x):
    c=chi2(x[None,:])[0]
    if not c<1e8: return -np.inf
    return -0.5*c+(lnL13(x) if DD else 0.)
def ptf(u): return LO+(HI-LO)*u
OUT='/workspace/out/q47/runs'; os.makedirs(OUT,exist_ok=True); t0=time.time()
ns=NestedSampler(logl,ptf,ND,nlive=nlive,sample='rslice',slices=int(os.environ.get('Q42_SL',8)),bound='multi',rstate=np.random.default_rng(seed))
ns.run_nested(dlogz=0.1,print_progress=False)
r=ns.results; w=np.exp(r.logwt-r.logz[-1]); w/=w.sum(); S=r.samples; tag=f'{reg}_{mode or "base"}_{DAT}_s{seed}'
np.savez(f'{OUT}/q47_{tag}.npz',S=S,w=w,names=names,logl=r.logl)
out=dict(data=DAT,reg=reg,mode=mode or 'base',seed=seed,nlive=nlive,logZ=float(r.logz[-1]),logZerr=float(r.logzerr[-1]),ess=float(w.sum()**2/(w**2).sum()),sec=time.time()-t0,chi2min=float(-2*r.logl.max()),names=list(names))
json.dump(out,open(f'{OUT}/q47_{tag}.json','w')); print(out)
