# Q8e: pooled, uv-corrected N/S contrast with the full budget; physical mapping (T_N/T_S, compact-component limits); T43/T21 N/S (no 3-2 array); equal-alpha continuum estimate.
# inputs: out/q8c/fac_samples.pkl (q8d_predict.py), out/q8b/mc_NS_samples.csv (q8b_diagnostics.py), out/q8/*.csv, RADEX grids out/q4c/grid4_*.npz
import numpy as np, pandas as pd, pickle
rng=np.random.default_rng(11); O='/workspace/out/q8c/'
FAC,CMP=pickle.load(open(O+'fac_samples.pkl','rb')); MC=pd.read_csv('/workspace/out/q8b/mc_NS_samples.csv')
PHYS=['taper_phys','taper_raw_phys']; FREE=['taper','taper_raw','sharp']; NR=40000
def pool(l,fam,q=None):
    out=[];cols=[]
    for V in fam:
        for qq in ((1.0,0.7) if q is None else (q,)):
            f,w,rc,reN,reS=FAC[(V,l,qq)]; 
            if w.sum()<30: continue
            k=rng.choice(len(f),NR//(len(fam)*(2 if q is None else 1)),p=w/w.sum()); out.append(f[k]); cols.append((reN[k],reS[k],rc[k]))
    return np.concatenate(out),[np.concatenate([c[i] for c in cols]) for i in range(3)]
rows=[]
def summ(name,x,ref=1.0):
    p=np.percentile(x,[2.5,16,50,84,97.5]); return dict(estimate=name,p2_5=p[0],p16=p[1],p50=p[2],p84=p[3],p97_5=p[4],P_gt1=float(np.mean(x>ref)))
POOLED={}
for l,col in (('HCN','NS_HCN'),('HCO+','NS_HCOp')):
    obs=MC[col].values
    for name,fam in (('physical Bmin>=15 m (alpha-bounded)',PHYS),('unbounded L (data-preferred, Xb4f ~ lossless)',FREE),('pooled 50/50',None)):
        if fam is None:
            fp,cp=pool(l,PHYS); ff,cf=pool(l,FREE); n=min(len(fp),len(ff)); fac=np.concatenate([fp[:n],ff[:n]]); cc=[np.concatenate([a[:n],b[:n]]) for a,b in zip(cp,cf)]
        else: fac,cc=pool(l,fam)
        o=rng.choice(obs,len(fac)); corr=o*fac; POOLED[(l,name)]=(corr,cc)
        r=summ(f'{l} N/S T43/T32 corrected: {name}',corr); r['uncorrected_MC_median']=np.median(o); r['factor_median']=np.median(fac); rows.append(r)
    r=summ(f'{l} N/S uncorrected (MC over beam, position, astrometry, noise)',obs); rows.append(r)
S=pd.DataFrame(rows); S.to_csv(O+'NS_corrected_pooled.csv',index=False); pd.set_option('display.width',250); print(S.round(3).to_string())
# continuum-calibrated estimate under equal spectral index: line N/S contrast / continuum Xb47/Xb4f N/S double ratio (1.198), with delta-alpha budget
cont=1.198; rows=[]
for l,col in (('HCN','NS_HCN'),('HCO+','NS_HCOp')):
    o=rng.choice(MC[col].values,NR)
    for da in (0.0,0.3,0.5):
        dal=rng.normal(0,da,NR) if da>0 else 0  # alpha_N-alpha_S
        corr=o/(cont/np.exp(dal*np.log(346.5/257.7)))
        rows.append(dict(line=l,sigma_dalpha=da,**{k:v for k,v in summ('x',corr).items() if k!='estimate'}))
E2=pd.DataFrame(rows); E2.to_csv(O+'NS_continuum_equal_alpha.csv',index=False); print(E2.round(3).to_string())
# T43/T21 (no 3-2 array involved): N/S = (T32/T21)(T43/T32) at the Xb55 beam
L=pd.read_csv('/workspace/out/q8b/ladder_Xb55beam.csv'); rows=[]
for sp in ('HCN','HCO+'):
    g=L[(L.line==sp)&(L.arr32=='3-2 Xb4f')]
    N=g[g.region=='N'].iloc[0]; S_=g[g.region=='S'].iloc[0]; B=g[g.region=='Bridge'].iloc[0]
    # T43/T21 = T43/T21 directly
    r=lambda x: (x.T43/x.T21, x.T43/x.T21*np.hypot(0.003/x.T43*0+0.0,0.0))
    ns=(N.T43/N.T21)/(S_.T43/S_.T21); nb=(N.T43/N.T21)/(B.T43/B.T21)
    rows.append(dict(line=sp,T43_T21_N=N.T43/N.T21,T43_T21_S=S_.T43/S_.T21,T43_T21_Bridge=B.T43/B.T21,NS=ns,N_over_Bridge=nb))
Q=pd.DataFrame(rows); Q.to_csv(O+'T43_T21_NS.csv',index=False); print(Q.round(3).to_string())
# physical mapping: RADEX one-zone grids (T,n) for the corrected extended-component ratios
res=[]
for m,mm in (('hcn','HCN'),('hco+','HCO+')):
    g=np.load(f'/workspace/out/q4c/grid4_{m}.npz'); TR=g['TR']; TK=g['TK']; LN=g['LN']; j=np.argmin(abs(g['LNC']-15.5)); R=TR[:,:,j,3]/TR[:,:,j,2]
    for name in ('physical Bmin>=15 m (alpha-bounded)','unbounded L (data-preferred, Xb4f ~ lossless)'):
        corr,cc=POOLED[(mm,name)]; reN,reS,rc=cc
        for ln in (5.0,5.5,6.0):
            i=np.argmin(abs(LN-ln)); col=R[:,i]
            def T(r): 
                ok=(r>=col.min())&(r<=col.max()); t=np.full(len(r),np.nan); t[ok]=10**np.interp(r[ok],col,np.log10(TK)); return t
            tN,tS=T(reN),T(reS); ok=np.isfinite(tN)&np.isfinite(tS)
            if ok.sum()<100: res.append(dict(mol=mm,family=name,logn=ln,note='ratio outside RADEX range',frac_ok=ok.mean())); continue
            x=tN[ok]/tS[ok]; res.append(dict(mol=mm,family=name,logn=ln,frac_ok=ok.mean(),rN_ext_p50=np.median(reN),rS_ext_p50=np.median(reS),TN_p50=np.median(tN[ok]),TS_p50=np.median(tS[ok]),TNoverTS_p16=np.percentile(x,16),TNoverTS_p50=np.median(x),TNoverTS_p84=np.percentile(x,84),P_TN_gt_TS=float(np.mean(x>1))))
P=pd.DataFrame(res); P.to_csv(O+'physical_mapping_extended.csv',index=False); print(P.round(2).to_string())
# compact-component-only limits when the extended gas has a common T43/T32 at N and S
rows=[]
for (V,l,q,re0),(rcN,rcS,aN,aS) in CMP.items():
    rows.append(dict(variant=V,line=l,q_compact=q,r_ext_common=re0,n_draw=len(rcN),rcN_p16=np.percentile(rcN,16),rcN_p50=np.median(rcN),rcN_p84=np.percentile(rcN,84),rcS_p16=np.percentile(rcS,16),rcS_p50=np.median(rcS),rcS_p84=np.percentile(rcS,84),rcN_over_rcS_p50=np.median(rcN/rcS),compact_share_43_N=aN,compact_share_43_S=aS,frac_rc_le_1=float(np.mean((rcN<=1.05)&(rcS<=1.05)))))
C=pd.DataFrame(rows); C.to_csv(O+'compact_component_limits.csv',index=False); print(C.round(2).to_string())
