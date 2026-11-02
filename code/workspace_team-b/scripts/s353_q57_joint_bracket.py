# Q57: joint N+S one-zone non-LTE (RADEX LVG) bracket: how much of the N/S HCN/HCO+ contrast can excitation/column/opacity differences carry vs [HCN]/[HCO+] abundance.
# Data (one 2.2" circular beam, +-150 km/s, K km/s; Q44 s261 data vector = Q29 s183/s185 1-0 and Q34 s192 12m+7m 3-2; 13C 1-0 from Q59-consistent Q29 s185_box):
#   HCN 1-0, HCN 3-2, HCO+ 1-0, HCO+ 3-2, H13CN 1-0, H13CO+ 1-0 at N and S.
# Model: RADEX LVG slab (dv=150, Tbg 2.73, H2), grids out/q39/grid_fine.npz (HCN, HCO+), out/q42/grid_h13cn.npz (H13CN own LAMDA file); H13CO+ = HCO+ file at N/X.
#   Free per nucleus: log Tk (10-300, 30 pts), log n (3-8 step 0.2), log N(HCO+) (13-18 step 0.1), log A = N(HCN)/N(HCO+) (-0.5..1.0 step 0.1), log X=12C/13C (1.3-2.3), beam-filling ln eta<=0 (flat, marginalised on a 0.1 grid).
#   Shared nuisance: 3-2 vs 1-0 band scale g (ln g = -0.1,0,+0.1 with weights .25,.5,.25), identical for N and S.
#   sigma^2 = stat^2 + (0.05 I_obs)^2 (cross-EB / 12m+7m scale) + (0.10 I_mod)^2 (LAMDA/geometry).  Linear-Gaussian likelihood (13C S/N is low).
# Hypotheses (flat priors in log of each parameter; Bayesian evidence by exact grid sums, shared parameters summed jointly for N and S):
#   Hexc : A, N(HCO+), X shared (only Tk, n, eta differ)         -> excitation only
#   Hcol : A, X shared;  N(HCO+) free per nucleus                -> excitation + column/opacity
#   HcolX: A shared; N(HCO+) and X free per nucleus               -> + nuclear-wide 13C/12C difference
#   Habun: Tk, n, N(HCO+), X shared; A free per nucleus           -> pure abundance
#   Hfull: everything free per nucleus (A free, X free)
#   HfullXs: A free, X shared
# usage: python3 -I s353_q57_joint_bracket.py [XLO XHI]  (prior range of log X, default 1.6 2.0 = X 40-100) -> out/q57/s353_<tag>.txt/json ; likelihood cube cached in out/q57/cube_*.npy
import sys,os,json,numpy as np
from multiprocessing import Pool
O='/workspace/out/q57/'; K=1.00927
FSYS=float(os.environ.get('Q57_SYS','0.05')); FMOD=float(os.environ.get('Q57_MOD','0.10')); SUF='' if (FSYS,FMOD)==(0.05,0.10) else f'_s{FSYS}_m{FMOD}'
g39=np.load('/workspace/out/q39/grid_fine.npz'); g42=np.load('/workspace/out/q42/grid_h13cn.npz')
TK=g39['TK']; LN=g39['LN']; LNC=g39['LNC']
assert np.allclose(TK,g42['TK']) and np.allclose(LN,g42['LN']) and np.allclose(LNC,g42['LNC'])
HCN=g39['hcn'][...,0]; HCO=g39['hcop'][...,0]; H13=g42['h13cn'][...,0]   # TR [T,n,N,line]
TAU={'hcn':g39['hcn'][...,2],'hco+':g39['hcop'][...,2]}
sel_n=np.arange(0,len(LN),2); LNs=LN[sel_n]
NOg=np.round(np.arange(13.0,18.0001,0.1),2); LAg=np.round(np.arange(-0.5,1.0001,0.1),2); LXg=np.round(np.arange(1.3,2.3001,0.1),2)
def ci(x): return np.round((x-LNC[0])/0.1).astype(int)
def I_of(tr): return 1.0645*150.0*tr   # K km/s (RJ Tb, 150 km/s slab; units as s62)
def model_cube(X):
    """returns I (6 lines) arrays [T,n,NO,A] for given log X"""
    T=len(TK); n=len(LNs); iO=ci(NOg); iA=ci(NOg[:,None]+LAg[None,:]); iHX=np.clip(ci(NOg[:,None]+LAg[None,:]-X),0,len(LNC)-1); iOX=np.clip(ci(NOg-X),0,len(LNC)-1)
    hn=HCN[:,sel_n]; ho=HCO[:,sel_n]; h13=H13[:,sel_n]
    out=np.empty((6,T,n,len(NOg),len(LAg)),np.float32)
    out[0]=I_of(hn[:,:,iA,0]); out[1]=I_of(hn[:,:,iA,1])
    out[2]=I_of(ho[:,:,iO,0])[...,None]; out[3]=I_of(ho[:,:,iO,1])[...,None]
    out[4]=I_of(h13[:,:,iHX,0]); out[5]=I_of(ho[:,:,iOX,0])[...,None]
    return np.nan_to_num(out,nan=1e-6).clip(1e-6)
def data(pos):
    import json
    r183=json.load(open('/workspace/out/q29/s183.json')); r192=json.load(open('/workspace/out/q34/s192.json'))['res']; box=json.load(open('/workspace/out/q29/s185_box.json'))
    g={'N':'N','S':'S'}[pos]; a=r183['2015.1.00993 X961|2.2|'+g]; c=r192['2p2|comb|'+g]
    I=[a['HCN']['h150_ap'][0],c['HCN32|ap|h150'][0],a['HCOp']['h150_ap'][0],c['HCOp32|ap|h150'][0]]
    e=[a['HCN']['h150_ap'][1],c['HCN32|ap|h150'][1],a['HCOp']['h150_ap'][1],c['HCOp32|ap|h150'][1]]
    hn=box[f'pm150|{g}|H13CN']; hc=box[f'pm150|{g}|H13COp']
    I+= [hn['amp'],hc['amp']]; e+=[hn['err'],hc['err']]
    return np.array(I)*K,np.array(e)*K
GS=[-0.1,0.0,0.1]; GW=np.array([.25,.5,.25]); LETA=np.arange(-3.0,0.0001,0.1)
def like(args):
    pos,X=args; I,e=data(pos); M=model_cube(X); out=np.zeros((3,)+M.shape[1:],np.float64)
    for gi,lg in enumerate(GS):
        sc=np.array([1,np.exp(lg),1,np.exp(lg),1,1])
        acc=np.zeros(M.shape[1:])
        for le in LETA:
            chi=0; lns=0
            for k in range(6):
                mod=M[k]*np.exp(le)*sc[k]
                s2=e[k]**2+(FSYS*I[k])**2+(FMOD*mod)**2
                chi=chi+(I[k]-mod)**2/s2; lns=lns+0.5*np.log(s2)
            acc+=np.exp(-0.5*chi-lns)
        out[gi]=acc/len(LETA)
    return out.astype(np.float64)
def build(pos):
    fn=O+f'cube_{pos}{SUF}.npy'
    if os.path.exists(fn): return np.load(fn)
    with Pool(11) as P: R=P.map(like,[(pos,X) for X in LXg])
    C=np.stack(R,1)   # [g,X,T,n,NO,A]
    np.save(fn,C.astype(np.float32)); return C
def main():
    xlo,xhi=(float(sys.argv[1]),float(sys.argv[2])) if len(sys.argv)>2 else (1.6,2.0)
    tag=f'X{int(round(10**xlo))}-{int(round(10**xhi))}'; out=[]
    def P(t): print(t,flush=True); out.append(t)
    CN=build('N').astype(np.float64); CS=build('S').astype(np.float64)
    xs=(LXg>=xlo-1e-6)&(LXg<=xhi+1e-6); CN=CN[:,xs]; CS=CS[:,xs]; nX=xs.sum()
    # axes: g,X,T,n,NO,A.  Unit prior weights; evidence = mean over grid (flat) -> ratios of evidences are what matter.
    nT,nn,nO,nA=CN.shape[2:]
    def ev(f):  # weighted over g
        return sum(GW[i]*f(i) for i in range(3))
    # clean explicit computation
    res={}
    def Zcalc(i):
        N_,S_=CN[i],CS[i]            # [X,T,n,NO,A]
        z={}
        z['Hexc']=(N_.sum((1,2))*S_.sum((1,2))).sum()/((nT*nn)**2*nO*nA*nX)
        z['Hcol']=(N_.sum((1,2,3))*S_.sum((1,2,3))).sum()/((nT*nn*nO)**2*nA*nX)
        z['HcolX']=(N_.sum((0,1,2,3))*S_.sum((0,1,2,3))).sum()/((nT*nn*nO*nX)**2*nA)
        z['Habun']=(N_.sum(4)*S_.sum(4)).sum()/((nT*nn*nO*nX)*nA**2)
        z['HfullXs']=(N_.sum((1,2,3,4))*S_.sum((1,2,3,4))).sum()/((nT*nn*nO*nA)**2*nX)
        z['Hfull']=N_.sum()*S_.sum()/((nT*nn*nO*nA*nX)**2)
        # pure abundance with X free per nucleus? skip
        return z
    zz=[Zcalc(i) for i in range(3)]
    Zf={k:sum(GW[i]*zz[i][k] for i in range(3)) for k in zz[0]}
    P(f'# prior log X {xlo}-{xhi} (X {10**xlo:.0f}-{10**xhi:.0f}); evidences relative to Hfull (Bayes factor = Z/Zfull); Occam penalties apply')
    for k,v in Zf.items(): P(f'{k:8s} lnZ-lnZfull = {np.log(v/Zf["Hfull"]):7.2f}   BF vs full = {v/Zf["Hfull"]:.3g}')
    # best-fit chi2 equivalents: (max likelihood = max of cube joint) per hypothesis (g=0 slice)
    i=1; N_,S_=CN[i],CS[i]
    mx={}
    mx['Hfull']=N_.max()*S_.max()
    mx['HfullXs']=(N_.max((1,2,3,4))*S_.max((1,2,3,4))).max()
    mx['Hexc']=(N_.max((1,2))*S_.max((1,2))).max()
    mx['Hcol']=(N_.max((1,2,3))*S_.max((1,2,3))).max()
    mx['HcolX']=(N_.max((1,2,3)).max(0)*S_.max((1,2,3)).max(0)).max()
    mx['Habun']=(N_.max(4)*S_.max(4)).max()
    P('max-likelihood ratio vs full (g=0): '+', '.join(f'{k}: {2*np.log(mx["Hfull"]/v):.1f}' for k,v in mx.items())+'  (= delta chi2 incl. sigma-normalisation terms)')
    # posteriors under Hfull: A_N, A_S, X, contrast
    pAN=CN.sum((0,1,2,3,4)); pAN/=pAN.sum(); pAS=CS.sum((0,1,2,3,4)); pAS/=pAS.sum()
    LA=LAg; q=lambda p,x:x[np.searchsorted(np.cumsum(p),[0.16,0.5,0.84])]
    P(f'Hfull: log A=N(HCN)/N(HCO+): N {10**q(pAN,LA)}, S {10**q(pAS,LA)}')
    J=pAN[:,None]*pAS[None,:]; ratio=(LA[:,None]-LA[None,:]).ravel(); w=J.ravel(); o=np.argsort(ratio); cw=np.cumsum(w[o])/w.sum()
    r=10**ratio[o][np.searchsorted(cw,[0.05,0.16,0.5,0.84,0.95])]
    P(f'Hfull: A_N/A_S (5,16,50,84,95%) = {np.round(r,2)}; P(A_N/A_S>1)={w[ratio>0.001].sum()/w.sum():.3f}; P(>1.3)={w[ratio>np.log10(1.3)].sum()/w.sum():.3f}')
    # same under HcolX/Hcol (A shared) not meaningful; give posterior of shared A and best-fit N/S contrast
    # abundance ratio required when X_N<X_S etc: posterior of X_N,X_S under Hfull
    for nm,C in (('N',CN),('S',CS)):
        pX=C.sum((0,2,3,4,5)); pX/=pX.sum(); xl=LXg[xs]; P(f'Hfull posterior X_{nm}: '+' '.join(f'{10**x:.0f}:{p:.2f}' for x,p in zip(xl,pX)))
    # model-predicted observables for the best fit of each hypothesis and for the excitation-only family
    return out,CN,CS,xs,mx,Zf
if __name__=='__main__':
    out,*_=main(); open(O+f's353_{sys.argv[1]+"-"+sys.argv[2] if len(sys.argv)>2 else "40-100"}.txt','w').write('\n'.join(out))
