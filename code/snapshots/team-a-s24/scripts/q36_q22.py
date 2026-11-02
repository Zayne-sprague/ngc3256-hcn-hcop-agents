# Q22 round 2: JOINT N+S fit with a SHARED extended zone (reviewer fix 1) and cross-nucleus leakage into the 10" apertures.
# usage: q22b_shared.py <env> <tbg ''|_tbg5> <seed>      env: sh (one shared envelope T,n,N,filling) | shf (shared T,n,N, separate fillings) | sep (independent envelopes) |
#        same (envelope = own-core physics with own filling, M1E-like) | shH (sh with Tc in [10,25], Te in [50,150] K) | none (cores only)
# env vars: Q22_CROSS (weight of the other nucleus' core in a 10" aperture, default 0.5: it sits at the aperture edge, nuclei ~5" apart), Q22_SG (10" band-scale sigma, default 0.10),
#           Q22_DROP (comma list of 10" datum indices 0-5 = HCN21,32,43,HCO21,32,43 to leave out), Q22_TMODE ('log' default | 'lin' = prior flat in T, applied by importance reweighting in post),
#           Q22_TAG, Q22_NP, Q22_LW (free 10" envelope-filling tweak not used), Q22_TMIN (Tk floor), Q12_EMAXN
# Q36 variant: arg2 'q36<scen>' -> core grids from out/q36 (freq-dependent background, per nucleus), envelope grid = RADEX Tbg=Q36_ENVSUF ('' =2.73 K).
# Data/likelihood per nucleus identical to q22_joint.py (1.7" ladders, 8 ratios + Tint, 10" ladders with common per-band scale and free line width).
import sys,os,warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.ndimage import map_coordinates
env,suf0,seed=sys.argv[1],sys.argv[2],int(sys.argv[3]); SCEN=suf0[3:] if suf0.startswith('q36') else None; suf=os.environ.get('Q36_ENVSUF','') if SCEN else ('' if suf0 in('','g273') else suf0); rng=np.random.default_rng(seed)
E=lambda k,d: float(os.environ.get(k,d))
CROSS=E('Q22_CROSS',0.5); SG10=E('Q22_SG',0.10); TAG=os.environ.get('Q22_TAG','b'); TMIN=E('Q20_TMIN',10); EMAXN=E('Q12_EMAXN',99)
DROP=[int(i) for i in os.environ.get('Q22_DROP','').split(',') if i!='']; KEEP=[i for i in range(6) if i not in DROP]
FLOOR,HSIG,SG=0.14,0.30,0.10; GG=np.exp(np.linspace(-3*SG,3*SG,9)); HH=np.exp(np.linspace(-3*HSIG,3*HSIG,9))
G_,H_=np.meshgrid(GG,HH); GV,HV=G_.ravel(),H_.ravel()
O='/workspace/out/q4c/'
G={m:np.load(O+f'grid4_{m}{suf}.npz') for m in ('hcn','hco+')}
TK,LN,LNC=G['hcn']['TK'],G['hcn']['LN'],G['hcn']['LNC']; LTK=np.log10(TK)
def lg(a): a=a.astype(float).copy(); a[~(a>0)|(a>1e5)]=np.nan; return np.log10(a)
LTS={'env':{m:np.stack([lg(G[m]['TR'][...,j]) for j in range(4)]) for m in G}}
ENVS=os.environ.get('Q36_ENVSCEN')   # envelope grid from out/q36 (freq-dependent background), else RADEX Tbg=Q36_ENVSUF
LTE=LTS['env'] if not ENVS else {m:np.stack([lg(np.load(f'/workspace/out/q36/grid_{m}_E_{ENVS}.npz')['TR'][...,j]) for j in range(4)]) for m in ('hcn','hco+')}
for r_ in 'NS':
    LTS['core_'+r_]=LTE if SCEN is None else {m:np.stack([lg(np.load(f'/workspace/out/q36/grid_{m}_{r_}_{SCEN}.npz')['TR'][...,j]) for j in range(4)]) for m in ('hcn','hco+')}
def idx(lT,ln,lN): return np.array([np.interp(lT,LTK,np.arange(len(TK))),np.interp(ln,LN,np.arange(len(LN))),np.interp(lN,LNC,np.arange(len(LNC)))])
LTS['env']=LTE
def TR(m,lT,ln,lN,key='env'): c=idx(lT,ln,lN); return np.stack([10**map_coordinates(LTS[key][m][j],c,order=1,mode='nearest') for j in range(4)])
I10={'N':{'hcn':[12.38,4.81,1.62],'hco+':[20.27,10.04,4.51]},'S':{'hcn':[11.66,4.26,1.35],'hco+':[20.89,9.95,4.28]}}
L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv'); DIL=pd.read_csv('/workspace/out/q4b/dilution.csv'); LIT=pd.read_csv(O+'lit_vs_ours.csv')
class R: pass
def setup(reg):
    r=R(); r.reg=reg
    r.D10=np.log(np.array(I10[reg]['hcn']+I10[reg]['hco+']))[KEEP]
    C=np.diag(np.tile(np.hypot(0.03,np.array([0.10,0.03,0.05])),2)**2)
    for b in range(3): a=np.zeros(6); a[b]=1; a[3+b]=1; C+=SG10**2*np.outer(a,a)
    r.C10=C[np.ix_(KEEP,KEEP)]; r.CI10=np.linalg.inv(r.C10)
    def lad(c,sp,q): x=L[(L.comp==c)&(L.species==sp)&(L.quantity==q)].iloc[0]; return x.value,(x.spread_hi-x.spread_lo)/2
    lit=lambda q: LIT[(LIT.region==reg)&(LIT.quantity==q)].iloc[0].value
    d=[]
    for sp in ('HCN','HCO+'):
        for q in ('3-2 [Xb4f]/2-1 (Tb)','4-3/3-2 [Xb4f 3-2] (Tb)'): v,s=lad(reg,sp,q); d.append((v,np.hypot(FLOOR*v,s)))
    r32=0.5*(lit('H R32 (K km/s)')+lit('ours R32 (K km/s, 1.7")')); d.append((r32,0.07*r32)); r10=lit('H R10 (K km/s)'); d.append((r10,0.05*r10))
    for sp in ('HCN','HCO+'):
        a=lit(f'H {sp} 3-2/1-0 (K km/s)'); b=lit(f'ours3-2/H1-0 {sp} (K km/s)'); m=0.5*(a+b); d.append((m,np.hypot(0.5*abs(a-b),0.05*m)))
    r.d=np.array(d); dl=DIL[DIL.region==reg].set_index('line'); K={'hcn':'HCN 3-2','hco+':'HCO+ 3-2'}
    r.Tint={m:dl.loc[k,'Tint_K'] for m,k in K.items()}; r.eT={m:0.5*(dl.loc[k,'Tint_84']-dl.loc[k,'Tint_16']) for m,k in K.items()}
    r.S={m:dl.loc[k,'s_fwhm_arcsec'] for m,k in K.items()}; r.D57={m:dl.loc[k,'Tb57_pred_compact']/dl.loc[k,'Tint_K'] for m,k in K.items()}
    return r
RG={'N':setup('N'),'S':setup('S')}; REGS=('N','S')
# ---- parameter layout
names=[]
for r in REGS: names+=[f'lTc_{r}',f'nc_{r}',f'NHc_{r}',f'NOc_{r}',f'lfc_{r}']+([f'lw_{r}'])
if env in('sh','shf','shH'): names+=['lTe','ne','NHe','NOe']+(['lwe'] if env!='shf' else ['lwe_N','lwe_S'])
elif env=='sep':
    for r in REGS: names+=[f'lTe_{r}',f'ne_{r}',f'NHe_{r}',f'NOe_{r}',f'lwe_{r}']
elif env=='same': names+=['lwe_N','lwe_S']
elif env!='none': raise SystemExit('env?')
HYP=env=='shH'
def bnd(n):
    b=n.split('_')[0]
    if b=='lw': return np.log(100.),np.log(300.)
    if b=='lTc': return (np.log10(10),np.log10(25)) if HYP else (np.log10(TMIN),np.log10(150))
    if b=='lTe': return (np.log10(50),np.log10(150)) if HYP else (np.log10(TMIN),np.log10(150))
    if b in('lfc','lwe'): return -3.,0.
    if b in('nc','ne'): return 3.,8.
    return 13.,18.5
LO=np.array([bnd(n)[0] for n in names]); HI=np.array([bnd(n)[1] for n in names]); ND=len(names); ix={n:i for i,n in enumerate(names)}
def unpack(x): return {n:x[:,i] for n,i in ix.items()}
def zones(p,r):
    """core and envelope physics (lT,n,NH,NO,we) for nucleus r"""
    c=dict(lT=p['lTc_'+r],n=p['nc_'+r],NH=p['NHc_'+r],NO=p['NOc_'+r])
    if env=='none': e=None
    elif env in('sh','shH'): e=dict(lT=p['lTe'],n=p['ne'],NH=p['NHe'],NO=p['NOe'],we=10**p['lwe'])
    elif env=='shf': e=dict(lT=p['lTe'],n=p['ne'],NH=p['NHe'],NO=p['NOe'],we=10**p['lwe_'+r])
    elif env=='sep': e=dict(lT=p['lTe_'+r],n=p['ne_'+r],NH=p['NHe_'+r],NO=p['NOe_'+r],we=10**p['lwe_'+r])
    else: e=dict(lT=c['lT'],n=c['n'],NH=c['NH'],NO=c['NO'],we=10**p['lwe_'+r])
    return c,e
def chi2_reg(p,r,rg,mode,full=False):
    o='S' if r=='N' else 'N'; ro=RG[o]
    c,e=zones(p,r); co,_=zones(p,o); fc=10**p['lfc_'+r]; fco=10**p['lfc_'+o]; P=len(fc)
    tr={}; tro={}; te={}
    for m,NN in (('hcn','NH'),('hco+','NO')):
        tr[m]=TR(m,c['lT'],c['n'],c[NN],'core_'+r); tro[m]=TR(m,co['lT'],co['n'],co[NN],'core_'+o)
        te[m]=TR(m,e['lT'],e['n'],e[NN]) if e else 0*tr[m]
    we=e['we'] if e else 0.
    Tb57={m:fc*(rg.D57[m]*1.)*tr[m]+we*te[m] for m in tr}
    d17={m:rg.S[m]**2/(rg.S[m]**2+1.7**2) for m in tr}
    T17={m:(fc*d17[m]*tr[m][[0,2]]) for m in tr}; E17={m:we*te[m][[0,2]] for m in tr}
    H,Oo=Tb57['hcn'],Tb57['hco+']
    pr=[H[2]/H[1],H[3]/H[2],Oo[2]/Oo[1],Oo[3]/Oo[2]]
    g=GV[:,None]; h=HV[:,None]
    T={m:h[None]*T17[m][:,None,:]+E17[m][:,None,:] for m in tr}   # (2,K,P): [1-0, 3-2]; h scales compact flux (as q22_joint)
    q=np.array([np.broadcast_to(a[None,:],(len(GV),P)) for a in pr]+[T['hcn'][1]/T['hco+'][1],T['hcn'][0]/T['hco+'][0],g*T['hcn'][1]/T['hcn'][0],g*T['hco+'][1]/T['hco+'][0]])
    c2=(((q-rg.d[:,0][:,None,None])/rg.d[:,1][:,None,None])**2).sum(0)+((np.log(GV)/SG)**2+(np.log(HV)/HSIG)**2)[:,None]
    if mode=='marg': mx=(-0.5*c2).max(0); c2=-2*(mx+np.log(np.exp(-0.5*c2-mx).sum(0)/len(GV)))
    else: c2=c2.min(0)
    for m in tr: c2=c2+((fc*tr[m][2]-rg.Tint[m])/np.hypot(rg.eT[m],0.10*rg.Tint[m]))**2
    # 10" aperture: own core + CROSS * other core (dilution at 10" with the other nucleus' own size) + envelope (filling we) ; line width lw free
    lm=[]
    for m in ('hcn','hco+'):
        d10=rg.S[m]**2/(rg.S[m]**2+100.); d10o=ro.S[m]**2/(ro.S[m]**2+100.)
        T10=fc*d10*tr[m]+CROSS*fco*d10o*tro[m]+we*te[m]
        lm+= [np.log(T10[j]) for j in (1,2,3)]
    lm=np.stack(lm)[KEEP]+np.log(1.0645)+p['lw_'+r][None,:]
    res=rg.D10[:,None]-lm; bad=~np.isfinite(res).all(0); res=np.nan_to_num(res)
    c10=np.einsum('ip,ij,jp->p',res,rg.CI10,res); c2=c2+np.where(bad,1e9,c10)
    return (c2,c10) if full else c2
def chi2(x,mode='marg'):
    p=unpack(x); bad=np.zeros(len(x),bool)
    for r in REGS:
        for nn,nk in (('NHc_','nc_'),('NOc_','nc_')): bad|=p[nn+r]>p[nk+r]+12.5
    if env in('sh','shf','shH'):
        bad|=(p['NHe']>p['ne']+12.5)|(p['NOe']>p['ne']+12.5)|(p['ne']>=p['nc_N'])|(p['ne']>=p['nc_S'])|(p['ne']>EMAXN)
    elif env=='sep':
        for r in REGS: bad|=(p['NHe_'+r]>p['ne_'+r]+12.5)|(p['NOe_'+r]>p['ne_'+r]+12.5)|(p['ne_'+r]>=p['nc_'+r])|(p['ne_'+r]>EMAXN)
    c2=sum(chi2_reg(p,r,RG[r],mode) for r in REGS)
    return np.where(bad|~np.isfinite(c2),1e9,c2)
def smc(NP,nmcmc=12):
    X=LO+(HI-LO)*rng.random((NP,ND)); c=chi2(X); beta=0.; logZ=0.; n=0
    while beta<1:
        def ess(db): w=np.exp(-0.5*db*(c-c.min())); return w.sum()**2/(w**2).sum()
        if ess(1-beta)>NP/2: db=1-beta
        else:
            a,b=0.,1-beta
            for _ in range(50):
                mid=0.5*(a+b); a,b=(mid,b) if ess(mid)>NP/2 else (a,mid)
            db=max(a,1e-12)
        w=np.exp(-0.5*db*(c-c.min())); logZ+=np.log(w.mean())-0.5*db*c.min(); w/=w.sum(); beta+=db
        ii=rng.choice(NP,NP,p=w); X=X[ii]; c=c[ii]
        Lc=np.linalg.cholesky(np.cov(X.T)+1e-9*np.eye(ND)); sc=2.38/np.sqrt(ND); acc=0
        for it in range(nmcmc):
            Y=X+sc*(rng.standard_normal((NP,ND))@Lc.T); inb=np.all((Y>=LO)&(Y<=HI),1)
            cy=np.where(inb,chi2(np.where(inb[:,None],Y,X)),1e9); a=(np.log(rng.random(NP))<-0.5*beta*(cy-c))&inb
            X[a]=Y[a]; c[a]=cy[a]; acc+=a.mean()
            if it==3 and acc/4<0.15: sc*=0.5
        n+=1
        if n>400: break
    return X,c,logZ
NP=int(E('Q22_NP',8000)); X,c,logZ=smc(NP)
o=np.argsort(c)[:300]; cp=chi2(X[o],'prof'); xb=X[o[np.argmin(cp)]][None,:]; cprof=float(cp.min())
pb=unpack(xb); out=dict(env=env,grid=(suf0 if SCEN else (suf or '_tbg2.73')),seed=seed,tag=TAG,cross=CROSS,sg10=SG10,drop=''.join(map(str,DROP)),tmin=TMIN,npar=ND,ndat=32-2*len(DROP),logZ=logZ,chi2_prof=cprof)
out['BIC']=cprof+ND*np.log(out['ndat'])
for r in REGS:
    c2,c10=chi2_reg(pb,r,RG[r],'prof',True); out['chi2_'+r]=float(c2[0]); out['chi2_10_'+r]=float(c10[0])
p=unpack(X)
for n_ in names:
    for k,pc in ((16,16),(50,50),(84,84)): out[f'{n_}_{k}']=float(np.percentile(p[n_],pc))
if 'lTe' in p: out['P(Te>=50)']=float((10**p['lTe']>=50).mean()); out['P(Te>Tc_N&Tc_S)']=float(((p['lTe']>p['lTc_N'])&(p['lTe']>p['lTc_S'])).mean())
tag=f"{TAG}_{env}_{suf0 if SCEN else (suf or 'tbg2.73')}_s{seed}"; os.makedirs('/workspace/out/q36/runs',exist_ok=True)
pd.DataFrame([out]).to_csv(f'/workspace/out/q36/runs/{tag}.csv',index=False); np.savez(f'/workspace/out/q36/runs/{tag}.npz',X=X,names=names,xb=xb)
print(tag,'chi2prof %.1f (N %.1f S %.1f; 10arcsec N %.1f S %.1f) k %d logZ %.2f BIC %.1f'%(cprof,out['chi2_N'],out['chi2_S'],out['chi2_10_N'],out['chi2_10_S'],ND,logZ,out['BIC']),
  ('Te %.0f P(Te>=50) %.2f'%(10**out['lTe_50'],out['P(Te>=50)']) if 'lTe_50' in out else ''))
