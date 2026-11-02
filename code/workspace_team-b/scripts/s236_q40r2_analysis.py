# Q40 round 2 analysis of s235 outputs (N, optionally C):
#  A. 13C 3-2/1-0 ratio ALONE (H13CO+, Harada 1.7"): chi2(Tk,n) with N(H13CO+) free; dchi2 cold vs rest.
#  B. per-datum decomposition (12C lines, 13C lines, band-scale/eta prior pi) of the best cell in each Tk branch (cold<=30, mid 30-100, hot>=100, hot capped at 150 K), for X and prior widths.
#  C. X marginalisation: lognormal priors on 12C/13C, posterior probability of Tk branches, Bayes factors; Tk-profile; Tk-cap sensitivity.
# usage: python3 -I s236_q40r2_analysis.py  -> out/q40r2/analysis.txt
import sys,os,numpy as np
sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']
import s235_q40r2_Xmarg as S
s=S.s; m=s.m; O=S.O; TK=m.TK; LN=m.LN
L=[]
def P(x=''): print(x); L.append(x)
# ---------- A. ratio alone
d=s.d; hc=d['hcop']      # [T,n,N,line,(TR,Tex,tau)]
LNC=m.LNC
tr=np.nan_to_num(hc[...,:2,0],nan=np.nan)
Iobs=(9.4,0.95,2.6,0.49); r_obs=Iobs[2]/Iobs[0]
ferr_obs=np.hypot(Iobs[1]/Iobs[0],Iobs[3]/Iobs[2])
for lab,fb in (('default band scales 10%/15%',np.hypot(.10,.15)),('wide 20%/25%',np.hypot(.2,.25))):
    ferr=np.hypot(np.hypot(ferr_obs,fb),0.10)       # obs + band-scale + 10% model/LAMDA
    ok=(tr[...,0]>1e-3)&(tr[...,1]>1e-3)
    R=np.where(ok,tr[...,1]/np.where(ok,tr[...,0],1),np.nan)           # [T,n,N]
    chi=(np.log(R/r_obs)/ferr)**2
    # N free
    c=np.nanmin(chi,axis=2); cm=np.nanmin(c); 
    P(f'A. ratio-only, {lab}: obs R32/10={r_obs:.3f} (frac err obs {ferr_obs:.2f}, total {ferr:.2f}); chi2min {cm:.2f}; min chi2 by branch (N13 free): cold(Tk<=30) {np.nanmin(c[TK<=30]):.2f}, mid {np.nanmin(c[(TK>30)&(TK<100)]):.2f}, hot {np.nanmin(c[TK>=100]):.2f}')
    # N13 tied to physically allowed: N13>=10^14.5 (thick 12C requires N13>~1e15/X) -> restrict to log N13 in [14.0,17.0] (HCO+ N(12)~1e15-1e17, X 25-100)
    sel=(LNC>=13.5)&(LNC<=17.0); c2=np.nanmin(chi[...,sel],axis=2)
    P(f'   restricted to logN13 13.5-17: cold {np.nanmin(c2[TK<=30]):.2f}, mid {np.nanmin(c2[(TK>30)&(TK<100)]):.2f}, hot {np.nanmin(c2[TK>=100]):.2f}')
    for t in (10,20,30,50,100,150,300):
        k=np.argmin(abs(TK-t)); P(f'   Tk={TK[k]:.0f}: min chi2 over n,N13 = {np.nanmin(c[k]):.2f}')
# ---------- B,C helpers
pri=lambda k: S.SIGS[k]
def load(pos,var,k,lx):
    f=O+f'fit_{pos}_{var}_s{k}_lX{int(round(lx*100))}.npz'
    return np.load(f) if os.path.exists(f) else None
def decomp(pos,var,k,lx,mask):
    S.BEAMSYS.clear()
    if var=='A2': S.BEAMSYS['H']=0.30
    q=S.iso2(var,pos); s.X_iso={'H':q['H'],'O':q['O'],'C':[]}
    s.MOD13['H']=s.lnI13(d['hcn'],lx); s.MOD13['O']=s.lnI13(d['hcop'],lx)
    m.setup(pos,*pri(k),['H','O','C'])
    r=load(pos,var,k,lx); cm=r['chimin'].copy(); cm[~mask,:]=np.inf
    it,jn=np.unravel_index(cm.argmin(),cm.shape)
    tot=m.X.pi.copy(); F={}
    for sp in 'HOC': F[sp]=m.chi_s(it,jn,sp); tot=tot+F[sp].min(1)
    ib=int(np.argmin(tot)); out=dict(Tk=TK[it],logn=LN[jn],chi2=tot[ib],pi=m.X.pi[ib],o3=m.X.o3[ib],o6=m.X.o6[ib],terms={})
    for sp in 'HOC':
        iN=int(np.argmin(F[sp][ib])); i0,i1=m.SPL[sp]; Mo=m.MOD[sp][it,jn][iN]
        off=(m.X.o3[ib],m.X.o6[ib]); y=m.X.y; v=m.X.v
        t12=[(y[i0]-Mo[0]-off[0])**2/v[i0],(y[i1]-Mo[1]-off[1])**2/v[i1]]
        out['terms'][sp]=dict(logN=m.NF[iN],t12=t12,pull12=[ (y[i0]-Mo[0]-off[0])/np.sqrt(v[i0]),(y[i1]-Mo[1]-off[1])/np.sqrt(v[i1])])
        t13=[]
        if sp in s.MOD13:
            M3=s.MOD13[sp][it,jn][iN]
            for (li,I,e) in s.X_iso[sp]:
                mod=np.exp(M3[li]+off[li]); s2=e**2+(0.10*I)**2+(0.10*mod)**2+(BEAMSYS_f(sp)*I)**2+(BEAMSYS_f(sp)*mod)**2
                t13.append((li,I,mod,(mod-I)**2/s2))
        out['terms'][sp]['t13']=t13
    return out
def BEAMSYS_f(sp): return S.BEAMSYS.get(sp,0.0)
def show(o,tag):
    h=o['terms']['H'];op=o['terms']['O'];c=o['terms']['C']
    s13=' '.join(f"[{('1-0' if li==0 else '3-2')} obs {I:.2f} mod {mod:.2f} chi2 {t:.2f}]" for li,I,mod,t in op['t13'])
    P(f"  {tag}: Tk={o['Tk']:.0f} logn={o['logn']:.1f} chi2={o['chi2']:.2f} = prior(band/eta) {o['pi']:.2f} + 12C HCN {sum(h['t12']):.2f} HCO+ {sum(op['t12']):.2f} CS {sum(c['t12']):.2f} + 13C {sum(t for *_,t in op['t13']):.2f}  (band offsets ln: {o['o3']:.2f},{o['o6']:.2f}); H13CO+ {s13}")
    return o
branches=[('cold Tk<=30',TK<=30),('mid 30<Tk<100',(TK>30)&(TK<100)),('hot Tk>=100',TK>=100),('Tk<=150',TK<=150)]
pos=os.environ.get('POS','N'); var=os.environ.get('VAR','H')
P(f'\nB. per-datum decomposition, position {pos}, variant {var}')
for k in range(3 if var=='H' else 2):
    for lx in (1.4,1.6,1.8,1.95):
        P(f' band-scale prior {pri(k)}, X={10**lx:.0f}')
        for lab,mk in branches:
            if load(pos,var,k,lx) is None: continue
            show(decomp(pos,var,k,lx,mk),lab)
# ---------- C. marginalisation
P(f'\nC. X marginalisation, position {pos}, variant {var}')
def tkprof(r): return r['chimin'].min(1)
NSP=3   # Wc = sum exp(-0.5 pi) prod_s sum_N exp(-0.5 (f_s-chibest)) -> Z = Wc exp(-0.5*NSP*chibest) (Q40r1 summary.txt used -0.5*chibest: wrong across X, harmless at fixed X)
def lnZ(r): return -0.5*NSP*r['chibest']+np.log(r['Wc'].sum())
priors={'LN(X0=45,s=0.2dex)':(np.log10(45),0.2),'LN(X0=60,s=0.2dex)':(np.log10(60),0.2),'LN(X0=40,s=0.15dex)':(np.log10(40),0.15),'flat 22-100 (logX)':None}
for k in range(3 if var=='H' else 2):
    R={lx:load(pos,var,k,lx) for lx in S.LOGX}; R={a:b for a,b in R.items() if b is not None}
    if len(R)<len(S.LOGX): P(f'  (incomplete set for prior {pri(k)}: {len(R)} X values)')
    lxs=np.array(sorted(R)); lz=np.array([lnZ(R[a]) for a in lxs])
    P(f' band-scale prior {pri(k)}: lnZ(X) = '+' '.join(f'{10**a:.0f}:{z:.1f}' for a,z in zip(lxs,lz)))
    for nm,pr in priors.items():
        lp=np.zeros(len(lxs)) if pr is None else -0.5*((lxs-pr[0])/pr[1])**2
        w=np.exp(lz+lp-(lz+lp).max()); 
        def Pbr(mask_T):
            num=0
            for a,wa in zip(lxs,w):
                r=R[a]; W=r['Wc']*np.exp(-0.5*r['chibest']+lz.max()) ; num+=wa*0   # placeholder
            return None
        # branch posterior: sum_X prior(X) * Z_branch(X)
        Zb=np.zeros((len(lxs),len(TK)))
        for i,a in enumerate(lxs):
            r=R[a]; Zb[i]=r['Wc'].sum(1)*np.exp(-0.5*NSP*r['chibest']-lz.max())
        pw=np.exp(lp)[:,None]*Zb; tot=pw.sum()
        br={lab:pw[:,mk].sum()/tot for lab,mk in (('cold<=30',TK<=30),('mid',(TK>30)&(TK<100)),('hot>=100',TK>=100),('<=150',TK<=150))}
        # prior mass of branches (log-uniform Tk grid)
        pxs=w/w.sum(); mean_lx=(pxs*lxs).sum(); 
        P(f'   X prior {nm}: P(Tk<=30)={br["cold<=30"]:.3f} P(30-100)={br["mid"]:.3f} P(>=100)={br["hot>=100"]:.3f} P(<=150)={br["<=150"]:.3f}; posterior X median {10**np.interp(.5,np.cumsum(pxs),lxs):.0f}; BF(hot>=100 : cold<=30)={br["hot>=100"]/br["cold<=30"]*(np.sum(TK<=30)/np.sum(TK>=100)):.1f} (prior-mass-corrected)')
    # Tk profile at X=45 and 60 (chi2 min over n,N), and cap sensitivity
    for lx in (1.65,1.8):
        if lx in R:
            pr_=tkprof(R[lx]); P(f'   profile chi2(Tk)-min at X={10**lx:.0f} (cell chi2): '+' '.join(f'{t:.0f}K:{c-pr_.min():.1f}' for t,c in zip(TK[::3],pr_[::3])))
    # best-X cap
    for cap in (30,50,100,150,300):
        b=min(R[a]['chimin'][TK<=cap].min() for a in lxs); b0=min(R[a]['chimin'].min() for a in lxs); P(f'   cap Tk<={cap}: best chi2 over X = {b:.2f} (dchi2 vs uncapped {b-b0:.2f})')
open(O+f'analysis_{pos}_{var}.txt','w').write('\n'.join(L))
