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
model,reg,dset,suf=sys.argv[1:5]; suf='' if suf=='g273' else suf; seed=int(sys.argv[5]) if len(sys.argv)>5 else 1
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
L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); DIL=pd.read_csv('/workspace/out/q4b/dilution.csv'); LIT=pd.read_csv(O+'lit_vs_ours.csv')
def lad(c,sp,q): r=L[(L.comp==c)&(L.species==sp)&(L.quantity==q)].iloc[0]; return r.value,(r.spread_hi-r.spread_lo)/2
def lit(q): return LIT[(LIT.region==reg)&(LIT.quantity==q)].iloc[0].value
d=[]
for sp in ('HCN','HCO+'):
    for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'):
        v,s=lad(reg,sp,q); d.append((v,np.hypot(0.14*v,s)))
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
GG=np.exp(np.linspace(-0.3,0.3,9)); HH=np.exp(np.linspace(-0.9,0.9,9)); SG,SH=0.10,0.30
# ---------------- parameter layout
two=model!='1z'; shT=model in ('2cT','2cTX'); shX=model=='2cTX'
names=['lTc','nc','NHc','NOc']
if two: names+=(['ne','NHe','NOe'] if shT else ['lTe','ne','NHe','NOe'])+['lwe']
names+=['lfc']
if shX: names.remove('NOe')
lo={'lT':1.,'n':3.,'N':13.,'lfc':-3.,'lwe':-3.}; hi={'lT':np.log10(150),'n':8.,'N':18.5,'lfc':0.,'lwe':0.}; lo['lT']=np.log10(10.)
def bnd(nm):
    k='lT' if nm.startswith('lT') else ('N' if nm.startswith('N') else ('lfc' if nm=='lfc' else ('lwe' if nm=='lwe' else 'n')))
    return lo[k],hi[k]
LO=np.array([bnd(n)[0] for n in names]); HI=np.array([bnd(n)[1] for n in names]); ND=len(names)
def unpack(x):
    p={n:x[:,i] for i,n in enumerate(names)}
    if two:
        if shT: p['lTe']=p['lTc']
        if shX: p['NOe']=p['NHe']-(p['NHc']-p['NOc'])
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
    Tb57={};Tb17={}
    for m in ('hcn','hco+'):
        c=tr[(m,'c')]; e=tr[(m,'e')] if two else 0.
        Tb57[m]=fc*D57[m]*c+we*e        # (4,P) coarse beam totals
        Tb17[m]=(fc*D17[m]*c,we*e)       # compact, extended at 1.7"
    return tr,Tb57,Tb17,fc,we
def chi2(x):
    p=unpack(x); P=x.shape[0]; bad=np.zeros(P,bool)
    for nn in ('NHc','NOc')+(('NHe','NOe') if two else ()): 
        nk='nc' if nn.endswith('c') else 'ne'; bad|=p[nn]>p[nk]+12.5
    if two: bad|=p['ne']>=p['nc']
    tr,Tb57,Tb17,fc,we=model_obs(p)
    H,O_=Tb57['hcn'],Tb57['hco+']
    pr=[H[2]/H[1],H[3]/H[2],O_[2]/O_[1],O_[3]/O_[2]]
    c17,e17={m:Tb17[m][0] for m in Tb17},{m:Tb17[m][1] for m in Tb17}
    best=None; lg_=[]
    for g in GG:
        for h in HH:
            T17={m:h*c17[m]+e17[m] for m in c17}
            q=pr+[T17['hcn'][2]/T17['hco+'][2],T17['hcn'][0]/T17['hco+'][0],g*T17['hcn'][2]/T17['hcn'][0],g*T17['hco+'][2]/T17['hco+'][0]]
            c2=sum(((q[i]-d[i,0])/d[i,1])**2 for i in range(8))+(np.log(g)/SG)**2+(np.log(h)/SH)**2
            lg_.append(-0.5*c2)
    lg_=np.array(lg_); mx=lg_.max(0); c2=-2*(mx+np.log(np.exp(lg_-mx).sum(0)/len(lg_)))
    if dset in ('r8a','full'):
        for m in ('hcn','hco+'): c2=c2+((fc*tr[(m,'c')][2]-Tint[m])/np.hypot(eTint[m],0.10*Tint[m]))**2
    if dset=='full' and two:
        for m in ('hcn','hco+'):
            e=we*tr[(m,'e')][2]; c2=c2+((e-EXC[m])/(0.4*EXC[m]))**2
    c2=np.where(bad|~np.isfinite(c2),1e9,c2)
    # tau (3-2) diagnostics are computed on demand
    return c2
# ---------------- tempered SMC
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
NPART=int(os.environ.get('Q12_NP',6000)); X,c,logZ,hist=smc(NPART)
# multiple restarts for the minimum: local optimisation from best particles
from scipy.optimize import minimize
order=np.argsort(c)[:20]; cm=c.min(); xbest=X[order[0]]
for i in order[:8]:
    r=minimize(lambda z:float(chi2(np.clip(z,LO,HI)[None,:])[0]),X[i],method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-4,maxiter=1500))
    if r.fun<cm: cm=r.fun; xbest=np.clip(r.x,LO,HI)
ndat={'r8':8,'r8a':10,'full':12 if two else 10}[dset]
if dset=='full' and not two: ndat=10
p=unpack(X); tr,Tb57,Tb17,fc,we=model_obs(p)
res=dict(model=model,region=reg,data=dset,grid=suf or '_tbg2.73',seed=seed,ndata=ndat,npar=ND,chi2_min=cm,logZ=logZ,BIC=cm+ND*np.log(ndat),nmcmc_levels=len(hist))
pb=unpack(xbest[None,:]); 
for n in pb: res['best_'+n]=float(pb[n][0])
res['dof']=ndat-ND   # nominal (nuisances g,h not counted)
def q(a,nm):
    for pp,qq in ((16,.16),(50,.5),(84,.84)): res[f'{nm}_{pp}']=float(np.nanpercentile(a,pp))
for n in ('lTc','nc','NHc','NOc'): q(p[n],n)
res['Tc_16'],res['Tc_50'],res['Tc_84']=[10**res[f'lTc_{k}'] for k in (16,50,84)]
for n in ('NHc','NOc'): pass
q(p['NHc']-p['NOc'],'logXc')
q(p['lfc'],'lfc')
if two:
    q(10**p['lTe'],'Te'); q(p['ne'],'ne'); q(p['NHe'],'NHe'); q(p['NOe'],'NOe'); q(p['lwe'],'lwe'); q(p['NHe']-p['NOe'],'logXe')
    res['P(Tc>Te)']=float((p['lTc']>p['lTe']).mean())
q(TAU32('hcn',p['lTc'],p['nc'],p['NHc']),'tau32H_c'); q(TAU32('hco+',p['lTc'],p['nc'],p['NOc']),'tau32O_c')
res['P(Tc<=40)']=float((10**p['lTc']<=40).mean()); res['P(Tc>=80)']=float((10**p['lTc']>=80).mean())
# compact fraction of coarse-beam 3-2 flux
fcH=(fc*D57['hcn']*tr[('hcn','c')][2])/Tb57['hcn'][2]; q(fcH,'compact_frac32_HCN')
os.makedirs('/workspace/out/q12/runs',exist_ok=True)
tag=f"{model}_{reg}_{dset}_{suf or 'tbg2.73'}_s{seed}"
pd.DataFrame([res]).to_csv(f'/workspace/out/q12/runs/{tag}.csv',index=False)
np.savez(f'/workspace/out/q12/runs/{tag}.npz',X=X,names=names)
print(tag,'chi2min %.2f ndata %d npar %d logZ %.2f  Tc %.0f [%.0f,%.0f] best:'%(cm,ndat,ND,logZ,res['Tc_50'],res['Tc_16'],res['Tc_84']),{k:round(v,2) for k,v in res.items() if k.startswith('best_')})
