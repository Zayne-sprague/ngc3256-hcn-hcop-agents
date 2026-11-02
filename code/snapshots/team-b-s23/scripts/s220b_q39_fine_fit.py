# Q39: Q35 one-gas / two-gas RADEX column-ratio fit on a FINE continuous grid (no 0.2 dex cells).
# Grid: logTk 30 pts (10-300 K, 0.051 dex) x log n 3-8 (0.1) x log N each species 11.5-19 (0.05 dex, cubic interpolation of ln I from the 0.1-dex RADEX grid out/q39/grid_fine.npz).
# Data: Harada+2018 Table 13, 1.7" peak pixel, K km/s: HCN 1-0, 3-2; HCO+ 1-0, 3-2; CS 2-1, 5-4 (as s210c). sigma_i^2 = (eI/I)^2 + 0.05^2.
# Nuisance: ln eta (filling, flat, <=0, common to all lines of one gas) and band scales z3 (1-0/2-1 lines, sigma 0.10), z6 (3-2/5-4 lines, sigma 0.15); both marginalised
#   by working with the per-band effective offsets o3=ln eta+z3, o6=ln eta+z6 on a 0.05 grid, with the analytic prior pi(o3,o6) = (o3-o6)^2/(s3^2+s6^2) - 2 ln Phi(-mu/s)  (eta integrated analytically).
# Model species factorise given (Tk,n,o3,o6): L = L_HCN(N_H) L_HCO(N_O) L_CS(N_C) exp(-pi/2) -> exact sums over a 0.05 dex N grid, no sampling.
# Outputs: marginal posteriors and profile likelihoods (chi2) of log N ratios (A=HCN/HCO+, CS/HCO+, CS/HCN) as histograms on a 0.05 dex lag grid; two-gas model (HCN+HCO+ in gas 1, CS in gas 2) variants.
# usage: python3 -I s220b_q39_fine_fit.py  -> out/q39/fit_<pos>_<variant>.npz and summary printed by s220c
import sys,numpy as np,pandas as pd,json,time
from multiprocessing import Pool
from scipy.interpolate import CubicSpline
from scipy.special import log_ndtr
import os
O='/workspace/out/q39/'
NCAP=float(os.environ.get('NCAP','99')); TAG=os.environ.get('TAG','')   # optional prior cap on log N (sensitivity run)
d=np.load(O+'grid_fine.npz'); TK=d['TK']; LN=d['LN']; LNC=d['LNC']
NF=np.round(np.arange(11.5,19.0001,0.05),2); nN=len(NF)
def lnI(a):   # a[T,n,N,line,(TR,Tex,tau)] -> ln(1.0645*150*TR) on NF, lines 0,1. TR<=0 (Tex<Tbg absorption at low n) and rare RADEX failures floored at 1e-3 K (never fit the bright data anyway).
    TR=np.nan_to_num(a[...,:2,0],nan=1e-3); l=np.log(1.0645*150*np.maximum(TR,1e-3))
    return CubicSpline(LNC,l,axis=2)(NF)
MOD={'H':lnI(d['hcn']),'O':lnI(d['hcop']),'C':lnI(d['cs'])}
H=pd.read_csv('/workspace/lit/harada_table13_parsed.csv'); H['sp']=H.species_transition.str.replace(' ','')
def g(p,s): r=H[(H.pos==p)&(H.sp.str.startswith(s))].iloc[0]; return r.I_Kkms,r.eI_Kkms
LINES=['HCN1','HCN3','HCO+1','HCO+3','CS2','CS5']
SPL={'H':(0,1),'O':(2,3),'C':(4,5)}
def nuis(s3,s6,step=0.05):
    if s3==0:   # all band scales fixed to 1: o3=o6=ln eta (flat <=0)
        o=np.round(np.arange(-3.6,0.0001,step),3); return o,o.copy(),np.zeros_like(o)
    o=np.round(np.arange(-3.6,0.6001,step),3); a,b=np.meshgrid(o,o,indexing='ij'); a=a.ravel(); b=b.ravel(); v=s3**2+s6**2
    mu=(a*s6**2+b*s3**2)/v; s=s3*s6/np.sqrt(v); pi=(a-b)**2/v-2*log_ndtr(-mu/s)
    k=pi<45; return a[k],b[k],pi[k]
class Ctx: pass
X=Ctx()
def setup(pos,s3,s6,species):
    X.y=np.array([np.log(g(pos,s)[0]) for s in LINES]); X.v=np.array([(g(pos,s)[1]/g(pos,s)[0])**2+0.05**2 for s in LINES])
    X.o3,X.o6,X.pi=nuis(s3,s6); X.sp=species
def chi_s(it,jn,s,o3=None,o6=None):
    M=MOD[s][it,jn]                       # (nN,2)
    if np.isnan(M).any(): return None
    i0,i1=SPL[s]; o3=X.o3 if o3 is None else o3; o6=X.o6 if o6 is None else o6
    f=(X.y[i0]-M[None,:,0]-o3[:,None])**2/X.v[i0]+(X.y[i1]-M[None,:,1]-o6[:,None])**2/X.v[i1]     # (C,nN)
    return np.where(NF[None,:]>NCAP,1e4,f)
def cell_min(c):
    it,jn=c; m=X.pi.copy()
    for s in X.sp:
        f=chi_s(it,jn,s)
        if f is None: return np.inf
        m=m+f.min(1)
    return m.min()
LAGS={'A':np.arange(-20,29),'CO':np.arange(-40,41),'CH':np.arange(-40,41)}   # in 0.05 dex steps; A: log N_HCN-log N_HCO+ -1..1.4; CS: N_CS/N_other -2..2
def pair(fX,fY,wX,wY,rest_w,rest_f,lags):
    """returns hist[L], prof[L] for stat N_X-N_Y = L steps."""
    h=np.zeros(len(lags)); p=np.full(len(lags),np.inf)
    for k,L in enumerate(lags):
        if L>=0: sx=slice(L,nN); sy=slice(0,nN-L)
        else: sx=slice(0,nN+L); sy=slice(-L,nN)
        h[k]=np.einsum('kn,kn,k->',wX[:,sx],wY[:,sy],rest_w); p[k]=(fX[:,sx]+fY[:,sy]+rest_f[:,None]).min()
    return h,p
def cell_pass2(args):
    c,chibest,thr=args; it,jn=c
    f={s:chi_s(it,jn,s) for s in X.sp}
    if any(v is None for v in f.values()): return None
    mins={s:f[s].min(1) for s in X.sp}; tot=X.pi+sum(mins.values())
    # exact cell weight
    w={s:np.exp(-0.5*(f[s]-chibest)) for s in X.sp}   # (C,nN)
    S={s:w[s].sum(1) for s in X.sp}; pw=np.exp(-0.5*X.pi)
    Wc=(pw*np.prod([S[s] for s in X.sp],0)).sum()
    out=dict(Wc=Wc,chimin=tot.min(),hist={},prof={})
    k=tot<chibest+thr
    out['best']=None
    if k.any():
        ib=np.argmin(tot); s_best={s:float(NF[np.argmin(f[s][ib])]) for s in X.sp}; out['best']=(tot[ib],s_best,X.o3[ib],X.o6[ib])
        specs=[('A','H','O')]+([('CO','C','O'),('CH','C','H')] if 'C' in X.sp else [])
        for key,a,b in specs:
            rest=[s for s in X.sp if s not in (a,b)]
            rw=pw[k]*np.prod([S[s][k] for s in rest],0) if rest else pw[k]
            rf=X.pi[k]+sum(mins[s][k] for s in rest)
            out['hist'][key],out['prof'][key]=pair(f[a][k],f[b][k],w[a][k],w[b][k],rw,rf,LAGS[key])
    # marginal over Tk, n done by caller (Wc); also per-species N marginals not needed
    return out
def run(pos,s3,s6,species):
    setup(pos,s3,s6,species); Pool_=Pool(16)
    cells=[(i,j) for i in range(len(TK)) for j in range(len(LN))]
    m=np.array(Pool_.map(cell_min,cells,chunksize=20)); cb=m.min(); print(pos,species,s3,'chi2min(data+prior)',cb,'at cell',cells[int(m.argmin())],flush=True)
    res=Pool_.map(cell_pass2,[(c,cb,45.) for c in cells],chunksize=20)
    Wc=np.array([r['Wc'] if r else 0 for r in res]).reshape(len(TK),len(LN))
    out=dict(chibest=cb,Wc=Wc)
    for key in LAGS:
        hs=[r['hist'][key] for r in res if r and key in r['hist']]
        if hs: out['h_'+key]=np.sum(hs,0); out['p_'+key]=np.min([r['prof'][key] for r in res if r and key in r['prof']],0)
    bs=[r['best'] for r in res if r and r['best']]; b=min(bs,key=lambda t:t[0]); ib=int(np.argmin(m))
    Pool_.close(); out['best']=dict(chi2=float(b[0]),Nbest=b[1],o3=float(b[2]),o6=float(b[3]),logT=float(np.log10(TK[cells[ib][0]])),logn=float(LN[cells[ib][1]]))
    return out
if __name__=='__main__':
    VARS={'one':['H','O','C'],'hh':['H','O'],'cs':['C']}
    jobs=sys.argv[1:] or ['one','hh']
    for sig in ((0.10,0.15),(0,0)):
      for pos in ('N','S','C'):
        for v in jobs:
            if v=='cs': continue
            t=time.time(); o=run(pos,sig[0],sig[1],VARS[v])
            tag=f'{pos}_{v}_{"sb" if sig[0] else "nosb"}{TAG}'
            np.savez(O+f'fit_{tag}.npz',**{k:v_ for k,v_ in o.items() if k!='best'}); json.dump(o['best'],open(O+f'fit_{tag}.json','w'))
            print(tag,o['best'],f'{time.time()-t:.0f}s',flush=True)
