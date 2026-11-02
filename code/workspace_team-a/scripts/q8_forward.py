# Q8-4: forward models. (A) size-only: can a compact+extended mixture with IDENTICAL intrinsic excitation (r_c, r_e) at N and S reproduce the N/S T43/T32 contrast?
#   Sizes/Tb of the compact and extended Gaussians from the 3-2 growth-curve posteriors (same importance sampler as scripts/q12c_geometry.py; re-run here to keep samples).
#   Ratio at beam B = sum_i a_i r_i (a_i = 3-2 flux fraction of component i in B; 4-3 size = q*3-2 size, q in {1,0.7}, for the 4-3 beam fractions).
#   r_c in {0.6,0.8,1.0}: 1.0 = RJ-thick thermal limit (T>>hnu/k); r_e from the Bridge/outer ratio 0.36 (HCN) / 0.48 (HCO+) at Xb4f (Bridge, Xb4f 3-2) -- and scanned.
# (B) Tkin vs density: invert the measured Tb ratios with the RADEX grid (out/q4c/grid4_*.npz, logN=15.5, 15.0, 16.0): dT at equal n, dlogn at equal T.
import numpy as np, pandas as pd
T=pd.read_csv('/workspace/out/q4/peak_tb.csv'); O17=pd.read_csv('/workspace/out/q4c/ours_1p7arcsec.csv'); rng=np.random.default_rng(8)
def F(s,b): return s**2/np.sqrt((s**2+b[0]**2)*(s**2+b[1]**2))
BXB4F=(6.2,3.7)  # exact Xb4f beam used in q8_matched (axis lengths; PA ignored: N-S axis along minor axis -> use a mean)
BMEAN=(np.sqrt(6.2*3.7),)*2
rows=[]; post={}
for reg in ('N','S'):
  for line in ('HCN 3-2','HCO+ 3-2'):
    pts=[]
    for o in ('Xb4b','Xb4d','Xb4f'):
        r=T[(T.region==reg)&(T.line==line)&(T.ous==o)].iloc[0]; b=tuple(float(x) for x in r.beam.split('x')); pts.append((o,b,r.Tb_peak_K,np.hypot(r.eTb_K,0.1*r.Tb_peak_K)))
    r=O17[(O17.region==reg)&(O17.line==line)].iloc[0]; pts.append(('s1.7',(1.7,1.7),r.Tpeak_K,np.hypot(r.rms_K,0.1*r.Tpeak_K)))
    N=3000000; s=np.exp(rng.uniform(np.log(.08),np.log(1.5),N)); se=np.exp(rng.uniform(np.log(1.5),np.log(25),N))
    Tc=np.exp(rng.uniform(np.log(1e-2),np.log(1e2),N)); Te=np.exp(rng.uniform(np.log(1e-2),np.log(1e2),N)); f=np.exp(np.log(0.85)+0.2*rng.standard_normal(N))
    c2=0
    for o,b,tb,e in pts:
        m=Tc*F(s,b)+Te*F(se,b)
        if o=='Xb4f': m=m*f
        c2=c2+((m-tb)/e)**2
    w=np.exp(-0.5*(c2-c2.min())); w/=w.sum(); idx=rng.choice(N,4000,p=w)
    post[(reg,line)]=dict(s=s[idx],se=se[idx],Tc=Tc[idx],Te=Te[idx])
def ac(reg,line,q,b=BMEAN):
    p=post[(reg,line)]; c=p['Tc']*F(p['s']*q,b); e=p['Te']*F(p['se']*q,b)   # 3-2 compact fraction at beam b if sizes scaled by q (q=1 for 3-2)
    return c/(c+e)
out=[]
for line,re_ in (('HCN',0.36),('HCO+',0.48)):
  for q in (1.0,0.7):
    n=min(len(post[('N',line+' 3-2')]['s']),len(post[('S',line+' 3-2')]['s']))
    # 3-2 fraction defines the 3-2 flux; 4-3 compact fraction uses sizes scaled by q (more compact 4-3 -> smaller compact share at coarse beam; extended also scaled by q)
    for rc in (0.6,0.8,1.0):
      for rext in (re_,0.25,0.45):
        aN3=ac('N',line+' 3-2',1.0); aS3=ac('S',line+' 3-2',1.0)
        aN4=ac('N',line+' 3-2',q); aS4=ac('S',line+' 3-2',q)
        # ratio of Tb(4-3)/Tb(3-2) in a mixture: numerator = Tb43_c*a4... with intrinsic peak ratio r_c (compact) and r_e (extended):
        # 3-2: Tc*Fc + Te*Fe ; 4-3: r_c*Tc*Fc' + r_e*Te*Fe' where ' are 4-3 filling factors (sizes*q)
        def ratio(reg):
            p=post[(reg,line+' 3-2')]; b=BMEAN
            n3=p['Tc']*F(p['s'],b)+p['Te']*F(p['se'],b); n4=rc*p['Tc']*F(p['s']*q,b)+rext*p['Te']*F(p['se']*q,b); return n4/n3
        rN,rS=ratio('N'),ratio('S'); c=rN/rS
        out.append(dict(line=line,q43=q,r_compact=rc,r_ext=rext,rN_med=np.median(rN),rS_med=np.median(rS),contrast_med=np.median(c),c16=np.percentile(c,16),c84=np.percentile(c,84),
            P_contrast_ge_obs=np.mean(c>=(1.33 if line=='HCN' else 1.13)),aN3_med=np.median(aN3),aS3_med=np.median(aS3)))
A=pd.DataFrame(out); A.to_csv('/workspace/out/q8/size_only_forward.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',200); print(A.round(3).to_string())
# (B) RADEX inversion
res=[]
obs={'HCN':{'Xb4f':(0.425,0.319),'Xb4d':(0.486,0.443),'Xb4b':(0.537,0.531)},'HCO+':{'Xb4f':(0.547,0.483),'Xb4d':(0.627,0.619),'Xb4b':(0.687,0.743)}}
for m,mm in (('hcn','HCN'),('hco+','HCO+')):
    g=np.load(f'/workspace/out/q4c/grid4_{m}.npz'); TR=g['TR']; TK=g['TK']; LN=g['LN']
    for lnc in (15.0,15.5,16.0):
        j=np.argmin(abs(g['LNC']-lnc)); R=TR[:,:,j,3]/TR[:,:,j,2]   # (T, logn)
        for arr,(rN,rS) in obs[mm].items():
            for ln in (5.0,5.5,6.0):
                i=np.argmin(abs(LN-ln)); col=R[:,i]
                Tof=lambda r: np.interp(r,col,np.log10(TK)) if col.min()<=r<=col.max() else np.nan
                res.append(dict(mol=mm,logN=lnc,arr32=arr,logn=ln,T_N=10**Tof(rN),T_S=10**Tof(rS),kind='T at equal n'))
            for Tk in (25,50):
                it=np.argmin(abs(TK-Tk)); row=R[it]
                nof=lambda r: np.interp(r,row,LN) if np.all(np.diff(row)>0) and row.min()<=r<=row.max() else np.nan
                res.append(dict(mol=mm,logN=lnc,arr32=arr,logn=nof(rN),T_N=Tk,T_S=Tk,kind=f'logn_N (S:{nof(rS):.2f}) at equal T',dlogn=nof(rN)-nof(rS)))
B=pd.DataFrame(res); B.to_csv('/workspace/out/q8/tkin_vs_density_inversion.csv',index=False); print(B.round(2).to_string())
# (C) beam filling of the compact (dense) component and compact share of the 3-2 peak Tb, per native beam (posterior median [16,84])
BEAMSN={'Xb4b 0.25"':(0.25,0.23),'Xb4d 0.93x0.74"':(0.93,0.74),'1.7" (smoothed)':(1.7,1.7),'Xb47 4.2x2.7"':(4.2,2.7),'Xb4f 5.7x3.5"':(5.7,3.5),'Xb55 9x6"':(9.0,6.0)}
rows=[]
for (reg,line),p in post.items():
    for bn,b in BEAMSN.items():
        fc=F(p['s'],b); fe=F(p['se'],b); share=p['Tc']*fc/(p['Tc']*fc+p['Te']*fe)
        q=lambda a: np.percentile(a,[16,50,84])
        rows.append(dict(region=reg,line=line,beam=bn,fill_compact_16=q(fc)[0],fill_compact_50=q(fc)[1],fill_compact_84=q(fc)[2],compact_share_16=q(share)[0],compact_share_50=q(share)[1],compact_share_84=q(share)[2]))
    rows.append(dict(region=reg,line=line,beam='size FWHM s (")',fill_compact_16=np.percentile(p['s'],16),fill_compact_50=np.median(p['s']),fill_compact_84=np.percentile(p['s'],84)))
    rows.append(dict(region=reg,line=line,beam='size FWHM se (")',fill_compact_16=np.percentile(p['se'],16),fill_compact_50=np.median(p['se']),fill_compact_84=np.percentile(p['se'],84)))
C=pd.DataFrame(rows); C.to_csv('/workspace/out/q8/filling_by_beam.csv',index=False); print(C.round(3).to_string())
