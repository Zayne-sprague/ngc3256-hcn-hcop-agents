# Q22 post-processing of the SMC posteriors (out/q22/runs/*.npz): P(T), physical quantities, abundance/path-length cut, H13CN(4-3) limit, core/envelope masses.
# Assumptions (labelled): X(HCN)=N(HCN)/N(H2) same in core and envelope; mu=2.8 m_H per H2; D=40 Mpc (194 pc/"); core path length L_c=N_H2/n <= 66 pc (item) ; envelope L_e<=500 pc;
# H13CN: HCN/H13CN(4-3) flux ratio of each zone (1-e^-t)/(1-e^-t/R), same excitation, R=12C/13C; limits S>4, N>8.3 (Q12), applied to the summed 1.7" 4-3 flux.
import numpy as np, pandas as pd, glob, sys
from scipy.ndimage import map_coordinates
O='/workspace/out/q4c/'
G={m:np.load(O+f'grid4_{m}.npz') for m in ('hcn','hco+')}
TK,LN,LNC=G['hcn']['TK'],G['hcn']['LN'],G['hcn']['LNC']; LTK=np.log10(TK)
def interp(a,lT,ln,lN):
    c=np.array([np.interp(lT,LTK,np.arange(len(TK))),np.interp(ln,LN,np.arange(len(LN))),np.interp(lN,LNC,np.arange(len(LNC)))])
    return map_coordinates(a.astype(float),c,order=1,mode='nearest')
mu=2.8*1.6605e-24; pc=3.086e18; Msun=1.989e33; PC=194.
SIZE={'S':0.455,'N':0.507}   # HCN 3-2 pedestal-inclusive core FWHM (out/q4b/dilution.csv), used by the fits at TSC=1
rows=[]
for tag in ('base','tmin30','tmin40','tmin20'):
  for model in ('M1E','M2T','M2','M2H'):
    for reg in ('S','N'):
        fs=sorted(glob.glob(f'/workspace/out/q22/runs/{tag}_{model}_{reg}_r8a+10_tbg2.73_s*.npz'))
        if not fs: continue
        Z=[np.load(f,allow_pickle=True) for f in fs]; names=list(Z[0]['names']); X=np.vstack([z['X'] for z in Z])
        p={n:X[:,i] for i,n in enumerate(names)}
        if model=='M2T' : p['lTe']=p['lTc']
        if model=='M1E': p['lTe']=p['lTc']; p['ne']=p['nc']; p['NHe']=p['NHc']; p['NOe']=p['NOc']
        fc=10**p['lfc']; we=10**p['lwe']
        s=SIZE[reg]
        tauH43c=interp(G['hcn']['tau'][...,3],p['lTc'],p['nc'],p['NHc']); tauH43e=interp(G['hcn']['tau'][...,3],p['lTe'],p['ne'],p['NHe'])
        tauH32c=interp(G['hcn']['tau'][...,2],p['lTc'],p['nc'],p['NHc']); tauH32e=interp(G['hcn']['tau'][...,2],p['lTe'],p['ne'],p['NHe'])
        TRc=interp(G['hcn']['TR'][...,3],p['lTc'],p['nc'],p['NHc']); TRe=interp(G['hcn']['TR'][...,3],p['lTe'],p['ne'],p['NHe'])
        d17=s**2/(s**2+1.7**2)
        lim=8.3 if reg=='N' else 4.0
        def h13ok(R):
            a=fc*d17*TRc; b=we*TRe
            rc=(1-np.exp(-tauH43c/R))/(1-np.exp(-tauH43c)); re_=(1-np.exp(-tauH43e/R))/(1-np.exp(-tauH43e))
            return (a+b)/(a*rc+b*re_)>lim
        q=lambda a,w=None:np.nanpercentile(a,[16,50,84]) if w is None else wq(a,w)
        def wq(a,w):
            i=np.argsort(a); c=np.cumsum(w[i])/w.sum(); return np.interp([.16,.5,.84],c,a[i])
        base=dict(tag=tag,model=model,region=reg,N=len(X),Tc=q(10**p['lTc']),Te=q(10**p['lTe']),P_Te_ge50=float((10**p['lTe']>=50).mean()),P_Tc_le25=float((10**p['lTc']<=25).mean()),P_Tc_le30=float((10**p['lTc']<=30).mean()),
                  nc=q(p['nc']),ne=q(p['ne']),logNHc=q(p['NHc']),logNHe=q(p['NHe']),eta_core=q(fc),we=q(we),tau32c=q(tauH32c),tau32e=q(tauH32e),tau43c=q(tauH43c),
                  P_H13=[float(h13ok(R).mean()) for R in (30,50,100)])
        rows.append(base); print({k:(np.round(v,2) if isinstance(v,np.ndarray) else v) for k,v in base.items()})
        for Xab in (1e-9,3e-9,1e-8,3e-8):
            NH2c=10**p['NHc']/Xab; Lc=NH2c/10**p['nc']/pc; NH2e=10**p['NHe']/Xab; Le=NH2e/10**p['ne']/pc
            ok=(Lc<=66.)&(Le<=500.) if model!='M1E' else (Lc<=66.)
            f=ok.mean()
            Mc=fc*NH2c*mu*1.133*(s*PC*pc)**2/Msun; Me=we*NH2e*mu*np.pi*(5*PC*pc)**2/Msun
            if model=='M1E': pass
            row=dict(tag=tag,model=model,region=reg,X=Xab,f_pass_path=float(f),dlnZ=float(np.log(max(f,1e-6))),Lc=q(Lc),Le=q(Le))
            if f>0.005:
                w=ok.astype(float)
                row.update(Tc=wq(10**p['lTc'],w),Te=wq(10**p['lTe'],w),Mc=wq(np.log10(Mc),w),Me=wq(np.log10(Me),w),fm_core=wq(Mc/(Mc+Me),w),P_Te_ge50=float(((10**p['lTe']>=50)*w).sum()/w.sum()),P_H13_R50=float((h13ok(50)*w).sum()/w.sum()))
            rows.append(row); print('   X=%g'%Xab,{k:(np.round(v,2) if isinstance(v,np.ndarray) else v) for k,v in row.items() if k not in('tag','model','region')})
pd.DataFrame(rows).to_csv('/workspace/out/q22/q22_post.csv',index=False)
