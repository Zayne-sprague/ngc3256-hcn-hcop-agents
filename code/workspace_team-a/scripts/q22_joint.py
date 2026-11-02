# Q22: joint core + envelope non-LTE model (RADEX LVG grids out/q4c) fitted to the 1.7" ladders (Q12b/Q20b data: 8 ratios), the 0.24" core T_b (Tint) and the team-b p0028 10" ladders
# (6 intensities, Q20f error model: per-line stat+aperture terms and a common per-band 10% flux scale shared by HCN and HCO+). Derived from q20b_eta.py (same 1.7" data, same nuisances g,h).
# usage: q22_joint.py <model> <region N|S> <data a|ab> <grid suffix ''|_tbg5> [seed]
#  model: M1  core only (Gaussian-core dilution at every aperture: the single-phase core of Q20 asked to explain all apertures)
#         M1E core + extended zone with the SAME (T,n,N_H,N_O) and its own filling (filling differs, physics identical)
#         M2T core + envelope, shared Tkin ;  M2 core + envelope, independent Tkin ;  M2H = M2 with Tc in [10,25] K and Te in [50,150] K (the hypothesis of the item)
#  data : a = r8a (1.7" ladders + Tint, 10 data) ; ab = a + 10" ladders (6 data) = 16
#  env: Q22_LR=1 -> envelope filling at 10" is we*10^lr, lr in [-1,0] free ; Q22_SG (10" band scale, 0.10) ; Q22_APS (aperture error scale) ; Q20_TSC (core size scale) ; Q12_NP particles
# Q12b (revised after review; superset of q12_twocomp.py, env-var variants Q12_*): true profile chi2 + pulls, prior/error/geometry sensitivity, u & size nuisances.
# Q12: two-component (compact + extended) non-LTE model for the N and S nuclei vs the one-zone model, fitted with a tempered SMC sampler (gives chi2_min, posteriors, log-evidence).
# usage: q12_twocomp.py <model> <region N|S> <data set> <grid suffix ''|_tbg5> [seed]
#  model: '1z' one-zone (T,n,NH,NO,fc) ; '2cT' two comp, SHARED Tkin ; '2cI' two comp, independent Tkin ; '2cTX' shared T and shared X=N(HCN)/N(HCO+)
#  data : 'r8'  = 8 ratio data of q4c_joint 'all' mode (coarse-beam deblended ladders d0-d3; 1.7" R32, R10, 3-2/1-0 HCN, HCO+ with flux-scale nuisance g)
#         'r8a' = r8 + compact intrinsic Tb (Tint, HCN & HCO+ 3-2; out/q4b/dilution.csv)                      (10 data)
#         'full'= r8a + extended-component 3-2 Tb at 5.7x3.5" (= Tb57 - compact prediction; 40% error)         (12 data)
# Model (per molecule m, line J; all RADEX radiation temperatures T_R from out/q4c/grid4_*.npz, dv=150 km/s, N/dv scaling NOT applied):
#   compact (size s_m from dilution table, intrinsic filling fc<=1):   Tb_c(m,J,beam b) = fc * TR_c(m,J) * D_m(b),  D(5.7x3.5")=Tb57_pred/Tint, D(1.7")=s^2/(s^2+1.7^2)
#   extended (filling we<=1 in the coarse beam, same at 1.7" up to nuisance u): Tb_e(m,J) = we * TR_e(m,J)
#   ladder ratios are formed from coarse-beam totals (c+e); 1.7" ratios from Tb_c(1.7")+u*Tb_e.  Nuisances (9x9 grid, marginalised): g (3-2 vs 1-0 scale, lognormal 0.10), h (compact 1.7" weight, lognormal 0.30)
#   priors: flat in log T (10-150 K), log n (3-8), log N_m (13-18.5), N_m<=n+12.5; n_e<n_c; log fc, log we flat in [-3,0]; u fixed 1.
import sys,os,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.ndimage import map_coordinates
model,reg,dset,suf=sys.argv[1:5]; dset={'a':'r8a','ab':'r8a+10'}.get(dset,dset); suf='' if suf=='g273' else suf; seed=int(sys.argv[5]) if len(sys.argv)>5 else 1
ENV=lambda k,dflt: float(os.environ.get(k,dflt))
NCMAX=ENV('Q12_NCMAX',8); FLO=ENV('Q12_FLO',-3); FLOOR=ENV('Q12_FLOOR',0.14); HSIG=ENV('Q12_HSIG',0.30); USIG=ENV('Q12_U',0); SKSIG=ENV('Q12_SK',0)
EBOUND=ENV('Q12_EMAXN',99); TAG=os.environ.get('Q12_TAG','base'); SIM=os.environ.get('Q12_SIM','')   # SIM='<truthfile.npz>:<seed>' -> simulated ratio data
LQ=ENV('Q20_LQ',0); TMIN=ENV('Q20_TMIN',10); TSC=ENV('Q20_TSC',1.0); OUTD='/workspace/out/q22/runs/'; LRF=ENV('Q22_LR',0); SG10=ENV('Q22_SG',0.10); APS=ENV('Q22_APS',1.0)
O='/workspace/out/q4c/'; rng=np.random.default_rng(seed)
G={m:np.load(O+f'grid4_{m}{suf}.npz') for m in ('hcn','hco+')}
TK,LN,LNC=G['hcn']['TK'],G['hcn']['LN'],G['hcn']['LNC']
def lg(a): a=a.astype(float).copy(); bad=~(a>0)|(a>1e5); a[bad]=np.nan; return np.log10(a)
LT={m:np.stack([lg(G[m]['TR'][...,j]) for j in range(4)]) for m in G}
TAU={m:G[m]['tau'][...,2].astype(float) for m in G}
LTK=np.log10(TK)
def idx(lT,ln,lN):
    return np.array([np.interp(lT,LTK,np.arange(len(TK))),np.interp(ln,LN,np.arange(len(LN))),np.interp(lN,LNC,np.arange(len(LNC)))])
def TR(m,lT,ln,lN):   # (4,P) radiation temperatures, NaN where invalid
    c=idx(lT,ln,lN); return np.stack([10**map_coordinates(LT[m][j],c,order=1,mode='nearest') for j in range(4)])
def TAU32(m,lT,ln,lN): return map_coordinates(TAU[m],idx(lT,ln,lN),order=1,mode='nearest')
# ---------------- data
I10={'N':{'hcn':[12.38,4.81,1.62],'hco+':[20.27,10.04,4.51]},'S':{'hcn':[11.66,4.26,1.35],'hco+':[20.89,9.95,4.28]}}   # team-b p0028 10" K km/s (2-1,3-2,4-3)
D10=np.log(np.array(I10[reg]['hcn']+I10[reg]['hco+']))
_ap=np.tile(np.hypot(0.03,np.array([0.10,0.03,0.05])*APS),2); C10=np.diag(_ap**2)
for _b in range(3):
    _a=np.zeros(6); _a[_b]=1; _a[3+_b]=1; C10+=SG10**2*np.outer(_a,_a)
CI10=np.linalg.inv(C10)
L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); DIL=pd.read_csv('/workspace/out/q4b/dilution.csv'); LIT=pd.read_csv(O+'lit_vs_ours.csv')
def lad(c,sp,q): r=L[(L.comp==c)&(L.species==sp)&(L.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
def lit(q): return LIT[(LIT.region==reg)&(LIT.quantity==q)].iloc[0].value
d=[]
for sp in ('HCN','HCO+'):
    for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'):
        v,s=lad(reg,sp,q); d.append((v,np.hypot(FLOOR*v,s)))
r32=0.5*(lit('H R32 (K km/s)')+lit('ours R32 (K km/s, 1.7")')); d.append((r32,0.07*r32))
r10=lit('H R10 (K km/s)'); d.append((r10,0.05*r10))
for sp in ('HCN','HCO+'):
    a=lit(f'H {sp} 3-2/1-0 (K km/s)'); b=lit(f'ours3-2/H1-0 {sp} (K km/s)'); m=0.5*(a+b); d.append((m,np.hypot(0.5*abs(a-b),0.05*m)))
d=np.array(d)   # HCN32/21, HCN43/32, HCO32/21, HCO43/32, R32, R10, HCN32/10, HCO32/10
dl=DIL[DIL.region==reg].set_index('line')
Tint={'hcn':dl.loc['HCN 3-2','Tint_K'],'hco+':dl.loc['HCO+ 3-2','Tint_K']}
eTint={m:0.5*(dl.loc[k,'Tint_84']-dl.loc[k,'Tint_16']) for m,k in (('hcn','HCN 3-2'),('hco+','HCO+ 3-2'))}
S={'hcn':dl.loc['HCN 3-2','s_fwhm_arcsec'],'hco+':dl.loc['HCO+ 3-2','s_fwhm_arcsec']}
D57={m:dl.loc[k,'Tb57_pred_compact']/dl.loc[k,'Tint_K'] for m,k in (('hcn','HCN 3-2'),('hco+','HCO+ 3-2'))}
D17={m:S[m]**2/(S[m]**2+1.7**2) for m in S}
EXC={m:dl.loc[k,'Tb57_excess_over_compact'] for m,k in (('hcn','HCN 3-2'),('hco+','HCO+ 3-2'))}
if TSC!=1.0:   # Q20: alternative core size: Tint -> TSC*Tint, size from dil_new=dil_old/TSC (beam 0.24"), D57/D17 rescaled consistently
    for m in S:
        so=S[m]; do=so**2/(so**2+0.24**2); dn=min(do/TSC,0.99); sn=0.24*np.sqrt(dn/(1-dn))
        D57[m]*=(sn**2/(sn**2+4.47**2))/(so**2/(so**2+4.47**2)); S[m]=sn; Tint[m]*=TSC; eTint[m]*=TSC
SG,SH=0.10,HSIG; GG=np.exp(np.linspace(-3*SG,3*SG,9)); HH=np.exp(np.linspace(-3*SH,3*SH,9))
BF=2.5   # beam geometric-mean (arcsec) for size rescaling of D57
# ---------------- parameter layout
two=model!='M1'; shT=model=='M2T'; shX=False; same=model=='M1E'; HYP=model=='M2H'; WITH10=dset.endswith('+10')
names=['lTc','nc','NHc','NOc']
if two and not same: names+=(['ne','NHe','NOe'] if shT else ['lTe','ne','NHe','NOe'])
if two: names+=['lwe']
if two and LRF: names+=['lr']
if WITH10: names+=['lw']
names+=['lfc']
if LQ>0: names+=['lq']
if two and USIG>0: names+=['lu']
if SKSIG>0: names+=['lk2','lk4']
lo={'lT':1.,'n':3.,'N':13.,'lfc':FLO,'lwe':FLO}; hi={'lT':np.log10(150),'n':NCMAX,'N':18.5,'lfc':0.,'lwe':0.}; lo['lT']=np.log10(TMIN)
def bnd(nm):
    if nm=='lq': return -LQ,LQ
    if nm=='lr': return -1.,0.
    if nm=='lw': return np.log(100.),np.log(300.)
    if nm=='lTc' and HYP: return np.log10(10),np.log10(25)
    if nm=='lTe' and HYP: return np.log10(50),np.log10(150)
    if nm=='lu': return -3*USIG,3*USIG
    if nm in('lk2','lk4'): return -3*SKSIG,3*SKSIG
    k='lT' if nm.startswith('lT') else ('N' if nm.startswith('N') else ('lfc' if nm=='lfc' else ('lwe' if nm=='lwe' else 'n')))
    return lo[k],hi[k]
LO=np.array([bnd(n)[0] for n in names]); HI=np.array([bnd(n)[1] for n in names]); ND=len(names)
def unpack(x):
    p={n:x[:,i] for i,n in enumerate(names)}
    if two:
        if shT or same: p['lTe']=p['lTc']
        if same: p['ne'],p['NHe'],p['NOe']=p['nc'],p['NHc'],p['NOc']
    return p
def model_obs(p):
    """returns dict of predictions and validity mask"""
    out={}; ok=np.ones(len(p['lTc']),bool)
    fc=10**p['lfc']
    tr={}
    for m,NN in (('hcn','NH'),('hco+','NO')):
        tr[(m,'c')]=TR(m,p['lTc'],p['nc'],p[NN+'c'])
        if two: tr[(m,'e')]=TR(m,p['lTe'],p['ne'],p[NN+'e'])
    we=10**p['lwe'] if two else 0.
    Tb57={};Tb17={};Tb10={}
    P=len(p['lTc']); kk=np.zeros((4,P)); kk[1]=p.get('lk2',0*p['lTc']); kk[3]=p.get('lk4',0*p['lTc']); kk[0]=kk[1]
    u=np.exp(p['lu']) if 'lu' in p else 1.
    for m in ('hcn','hco+'):
        c=tr[(m,'c')]; e=tr[(m,'e')] if two else 0.*c
        sj=S[m]*np.exp(kk); f0=S[m]**2/(S[m]**2+BF**2*0+4.47**2); fj=sj**2/(sj**2+4.47**2)
        d57=D57[m]*fj/f0; d17=sj**2/(sj**2+1.7**2)
        Tb57[m]=fc*d57*c+we*e        # (4,P) coarse beam totals
        Tb17[m]=(fc*d17*c,u*we*e)       # compact, extended at 1.7"
        d10=sj**2/(sj**2+10.**2); we10=we*(10**p['lr'] if 'lr' in p else 1.)
        Tb10[m]=(fc*d10*c, we10*e)
    return tr,Tb57,Tb17,fc,we,Tb10
NGAUSS=(1 if (two and USIG>0) else 0)+(2 if SKSIG>0 else 0)
LOGZ_CORR=NGAUSS*np.log(6/np.sqrt(2*np.pi))   # truncated +-3 sigma uniform+Gaussian-penalty -> proper Gaussian prior
def mk_grid(n):
    g=np.exp(np.linspace(-3*SG,3*SG,n)); h=np.exp(np.linspace(-3*SH,3*SH,n)); G_,H_=np.meshgrid(g,h); return G_.ravel(),H_.ravel()
GRID9=mk_grid(9)
def pred(p,gv,hv):
    """8 predicted ratios (8,K,P) for particles p and nuisance values gv,hv (K,)"""
    tr,Tb57,Tb17,fc,we,Tb10=model_obs(p)
    H,O_=Tb57['hcn'],Tb57['hco+']
    pr=[H[2]/H[1],H[3]/H[2],O_[2]/O_[1],O_[3]/O_[2]]
    K,P=len(gv),H.shape[1]; g=gv[:,None]; h=hv[:,None]
    qv=10**p['lq'][None,:] if 'lq' in p else 1.
    T={m:{J:h*Tb17[m][0][J][None,:]+Tb17[m][1][J][None,:] for J in (0,2)} for m in Tb17}
    for J in (0,2): T['hco+'][J]=T['hco+'][J]*qv
    q=[np.broadcast_to(a[None,:],(K,P)) for a in pr]+[T['hcn'][2]/T['hco+'][2],T['hcn'][0]/T['hco+'][0],g*T['hcn'][2]/T['hcn'][0],g*T['hco+'][2]/T['hco+'][0]]
    return np.array(q),tr,fc,we,Tb10
def chi2(x,mode='marg',grid=None):
    p=unpack(x); P=x.shape[0]; bad=np.zeros(P,bool)
    for nn in ('NHc','NOc')+(('NHe','NOe') if two else ()):
        nk='nc' if nn.endswith('c') else 'ne'; bad|=p[nn]>p[nk]+12.5
    if two and not same: bad|=p['ne']>=p['nc']
    if two and not same: bad|=p['ne']>EBOUND
    gv,hv=grid if grid is not None else GRID9
    q,tr,fc,we,Tb10=pred(p,gv,hv)
    c2=(((q-d[:,0][:,None,None])/d[:,1][:,None,None])**2).sum(0)+((np.log(gv)/SG)**2+(np.log(hv)/SH)**2)[:,None]
    if mode=='marg':
        mx=(-0.5*c2).max(0); c2=-2*(mx+np.log(np.exp(-0.5*c2-mx).sum(0)/len(gv)))
    else: c2=c2.min(0)
    extra=0.
    if two and USIG>0: extra=extra+(p['lu']/USIG)**2
    if SKSIG>0: extra=extra+(p['lk2']/SKSIG)**2+(p['lk4']/SKSIG)**2
    c2=c2+extra
    if WITH10:
        lm=np.log(np.stack([Tb10['hcn'][0][j]+Tb10['hcn'][1][j] for j in (1,2,3)]+[Tb10['hco+'][0][j]+Tb10['hco+'][1][j] for j in (1,2,3)]))+np.log(1.0645)+p['lw'][None,:]
        r10=D10[:,None]-lm; bad10=~np.isfinite(r10).all(0)
        c2=c2+np.where(bad10,1e9,np.einsum('ip,ij,jp->p',np.nan_to_num(r10),CI10,np.nan_to_num(r10)))
    if dset.startswith('r8a') or dset=='full':
        for m in ('hcn','hco+'): c2=c2+((fc*(10**p['lq'] if (m=='hco+' and 'lq' in p) else 1.)*tr[(m,'c')][2]-Tint[m])/np.hypot(eTint[m],0.10*Tint[m]))**2
    if dset=='full' and two:
        for m in ('hcn','hco+'):
            e=we*tr[(m,'e')][2]; c2=c2+((e-EXC[m])/(0.4*EXC[m]))**2
    return np.where(bad|~np.isfinite(c2),1e9,c2)

def smc(NP=6000,nmcmc=12):
    X=LO+(HI-LO)*rng.random((NP,ND)); c=chi2(X)
    # reject invalid priors (truncation counts as prior mass): keep them with c=1e9
    beta=0.; logZ=0.; hist=[]
    while beta<1:
        lo_,hi_=0.,1-beta
        def ess(db): w=np.exp(-0.5*db*(c-c.min())); return w.sum()**2/(w**2).sum()
        if ess(1-beta)>NP/2: db=1-beta
        else:
            a,b=0.,1-beta
            for _ in range(50):
                mid=0.5*(a+b); 
                if ess(mid)>NP/2: a=mid
                else: b=mid
            db=max(a,1e-12)
        w=np.exp(-0.5*db*(c-c.min())); logZ+=np.log(w.mean())-0.5*db*c.min(); w/=w.sum()
        beta+=db; ii=rng.choice(NP,NP,p=w); X=X[ii]; c=c[ii]
        cov=np.cov(X.T)+1e-9*np.eye(ND); Lc=np.linalg.cholesky(cov); sc=2.38/np.sqrt(ND)
        acc=0
        for it in range(nmcmc):
            Y=X+sc*(rng.standard_normal((NP,ND))@Lc.T)
            inb=np.all((Y>=LO)&(Y<=HI),1); Yc=np.where(inb[:,None],Y,X); cy=chi2(Yc); cy=np.where(inb,cy,1e9)
            a=np.log(rng.random(NP))<-0.5*beta*(cy-c); a&=inb; X[a]=Y[a]; c[a]=cy[a]; acc+=a.mean()
            if it==3 and acc/4<0.15: sc*=0.5
        hist.append((beta,acc/nmcmc))
        if len(hist)>300: break
    return X,c,logZ,hist

if SIM:
    tf,ss=SIM.split(':'); qt=np.load(tf)['qbest']; d=d.copy(); d[:,0]=qt+np.random.default_rng(int(ss)).standard_normal(8)*d[:,1]
NPART=int(os.environ.get('Q12_NP',6000)); X,c,logZ,hist=smc(NPART); logZ+=LOGZ_CORR
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2d
G13=mk_grid(13); G61=mk_grid(61)
order=np.argsort(c)[:20]; xbest=X[order[0]]; cm=1e9
for i in order[:8]:
    r=minimize(lambda z:float(chi2(np.clip(z,LO,HI)[None,:],'prof',G13)[0]),X[i],method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-4,maxiter=1500))
    if r.fun<cm: cm=r.fun; xbest=np.clip(r.x,LO,HI)
cprof=float(chi2(xbest[None,:],'prof',G61)[0])
# pulls at the profile optimum
pb=unpack(xbest[None,:]); 
G61g,G61h=G61; q,_,_,_,_=pred(pb,G61g,G61h)
c2all=(((q-d[:,0][:,None,None])/d[:,1][:,None,None])**2).sum(0)[:,0]+(np.log(G61g)/SG)**2+(np.log(G61h)/SH)**2
k=int(np.argmin(c2all)); gb,hb=G61g[k],G61h[k]; qb=q[:,k,0]; pulls=(qb-d[:,0])/d[:,1]
ndat={'r8':8,'r8a':10,'r8a+10':16,'full':12 if two else 10}[dset]
NEFF=ND-(1 if dset=='r8' else 0)-(1 if (LQ>0 and dset=='r8') else 0)-(NGAUSS)          # r8: fc (and we) only enter via ratios -> one fewer constrained dof; Gaussian-prior nuisances not counted
dof=ndat-NEFF
res=dict(model=model,region=reg,data=dset,grid=suf or '_tbg2.73',seed=seed,tag=TAG,ndata=ndat,npar=ND,neff=NEFF,dof=dof,chi2_prof=cprof,chi2_ratios=float((pulls**2).sum()),
    p_prof=float(chi2d.sf(cprof,dof)) if dof>0 else np.nan,logZ=logZ,g_best=gb,h_best=hb,AIC=cprof+2*NEFF,BIC=cprof+NEFF*np.log(ndat))
for i,pl in enumerate(pulls): res[f'pull{i}']=float(pl)
if WITH10:
    _p=pb; _tr,_,_,_,_,_T10=model_obs(_p)
    _lm=np.log(np.stack([_T10['hcn'][0][j]+_T10['hcn'][1][j] for j in (1,2,3)]+[_T10['hco+'][0][j]+_T10['hco+'][1][j] for j in (1,2,3)]))[:,0]+np.log(1.0645)+_p['lw'][0]
    _r=D10-_lm
    for i in range(6): res[f'pull10_{i}']=float(_r[i]/np.sqrt(C10[i,i]))
    res['chi2_10']=float(_r@CI10@_r)
for n in pb: res['best_'+n]=float(pb[n][0])
p=unpack(X); tr,Tb57,Tb17,fc,we,Tb10=model_obs(p)
def q_(a,nm):
    for pp,qq in ((16,.16),(50,.5),(84,.84)): res[f'{nm}_{pp}']=float(np.nanpercentile(a,pp))
for n in ('lTc','nc','NHc','NOc'): q_(p[n],n)
res['Tc_16'],res['Tc_50'],res['Tc_84']=[10**res[f'lTc_{k}'] for k in (16,50,84)]
q_(p['NHc']-p['NOc'],'logXc'); q_(p['lfc'],'lfc')
if LQ>0: q_(p['lq'],'lq')
q_(p['lfc']+np.log10(1.),'lfc2')
q_(np.log10(fc*np.minimum(1,10**(p['lq'] if LQ>0 else 0))),'lfc_hco')
q_(np.log10(fc*tr[('hcn','c')][2]),'lTbc_hcn'); q_(np.log10(fc*tr[('hco+','c')][2]),'lTbc_hco')
if two:
    q_(10**p['lTe'],'Te'); q_(p['ne'],'ne'); q_(p['NHe'],'NHe'); q_(p['lwe'],'lwe'); res['P(Tc>Te)']=float((p['lTc']>p['lTe']).mean())
q_(TAU32('hcn',p['lTc'],p['nc'],p['NHc']),'tau32H_c')
res['P(Tc<=40)']=float((10**p['lTc']<=40).mean()); res['P(Tc>=80)']=float((10**p['lTc']>=80).mean())
os.makedirs(OUTD,exist_ok=True)
tag=f"{TAG}_{model}_{reg}_{dset}_{suf or 'tbg2.73'}_s{seed}"+('_sim'+SIM.split(':')[1]+'_'+os.path.basename(SIM.split(':')[0])[:-4] if SIM else '')
pd.DataFrame([res]).to_csv(OUTD+tag+'.csv',index=False)
np.savez(OUTD+tag+'.npz',X=X,names=names,qbest=qb,pulls=pulls,dobs=d)
print(tag,'chi2prof %.2f dof %d p %.3f logZ %.2f Tc %.0f [%.0f,%.0f] g %.2f h %.2f'%(cprof,dof,res['p_prof'],logZ,res['Tc_50'],res['Tc_16'],res['Tc_84'],gb,hb),np.round(pulls,2))
