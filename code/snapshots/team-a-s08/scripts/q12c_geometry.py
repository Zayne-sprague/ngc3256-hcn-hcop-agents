# Q12c-1: matched-beam geometry from the data. Growth curve of the peak Tb of HCN/HCO+ 3-2 at N and S across four beams
# (Xb4b 0.25", Xb4d 0.93x0.74", Xb4d smoothed to 1.7" round [out/q4c/ours_1p7arcsec.csv], Xb4f 5.73x3.54"; out/q4/peak_tb.csv).
# Model: compact Gaussian (FWHM s, intrinsic Tc) + extended Gaussian (FWHM se, Te); peak Tb in an elliptical beam = T*s^2/sqrt((s^2+bmaj^2)(s^2+bmin^2)).
# Errors: channel rms + 10% abs. flux scale per OUS (independent; 1.7" point shares the Xb4d scale -> 10% added in quadrature only, correlation ignored as small vs the others);
# flux-recovery nuisance f on the Xb4f point (lognormal, median 0.85, sigma 0.2; team-b p0003/p0008 ~0.7 centre, assumption).
# Prior: log-uniform s in [0.08,1.5]", se in [1.5,25]", Tc, Te in [1e-2,1e2] K; importance sampling 4e6 draws. Derived: u = Tb_e(1.7")/Tb_e(4.47" round) (the model's extended 1.7"/coarse filling ratio),
# ext fraction at 1.7", size s (for per-line compact size prior), lognormal summary of u. Output out/q12c/geometry.csv
import numpy as np, pandas as pd
T=pd.read_csv('/workspace/out/q4/peak_tb.csv'); O17=pd.read_csv('/workspace/out/q4c/ours_1p7arcsec.csv')
rng=np.random.default_rng(3)
def F(s,b): return s**2/np.sqrt((s**2+b[0]**2)*(s**2+b[1]**2))
rows=[]
for reg in ('N','S'):
  for line in ('HCN 3-2','HCO+ 3-2'):
    pts=[]
    for o in ('Xb4b','Xb4d','Xb4f'):
        r=T[(T.region==reg)&(T.line==line)&(T.ous==o)].iloc[0]; b=tuple(float(x) for x in r.beam.split('x')); pts.append((o,b,r.Tb_peak_K,np.hypot(r.eTb_K,0.1*r.Tb_peak_K)))
    r=O17[(O17.region==reg)&(O17.line==line)].iloc[0]; pts.append(('s1.7',(1.7,1.7),r.Tpeak_K,np.hypot(r.rms_K,0.1*r.Tpeak_K)))
    N=4000000; ls=rng.uniform(np.log(.08),np.log(1.5),N); lse=rng.uniform(np.log(1.5),np.log(25),N)
    lTc=rng.uniform(np.log(1e-2),np.log(1e2),N); lTe=rng.uniform(np.log(1e-2),np.log(1e2),N); lf=np.log(0.85)+0.2*rng.standard_normal(N)
    s,se,Tc,Te,f=np.exp(ls),np.exp(lse),np.exp(lTc),np.exp(lTe),np.exp(lf)
    c2=0
    for o,b,tb,e in pts:
        m=Tc*F(s,b)+Te*F(se,b); 
        if o=='Xb4f': m=m*f
        c2=c2+((m-tb)/e)**2
    w=np.exp(-0.5*(c2-c2.min())); w/=w.sum(); ess=1/(w**2).sum()
    u=(F(se,(1.7,1.7))/F(se,(4.47,4.47)))
    Tb17c=Tc*F(s,(1.7,1.7)); Tb17e=Te*F(se,(1.7,1.7)); fe17=Tb17e/(Tb17c+Tb17e)
    def q(a,p=(16,50,84)):
        i=np.argsort(a); cw=np.cumsum(w[i]); return [a[i][np.searchsorted(cw,x/100)] for x in p]
    lu=np.log(u); mu=np.sum(w*lu); sd=np.sqrt(np.sum(w*(lu-mu)**2))
    rows.append(dict(region=reg,line=line,chi2min=c2.min(),ess=ess,s_16=q(s)[0],s_50=q(s)[1],s_84=q(s)[2],se_16=q(se)[0],se_50=q(se)[1],se_84=q(se)[2],
        u_16=q(u)[0],u_50=q(u)[1],u_84=q(u)[2],lnu_mean=mu,lnu_sd=sd,fext17_16=q(fe17)[0],fext17_50=q(fe17)[1],fext17_84=q(fe17)[2],
        Tc_50=q(Tc)[1],Te_50=q(Te)[1],f_50=q(f)[1],dlns_sd=np.sqrt(np.sum(w*(ls-np.sum(w*ls))**2))))
R=pd.DataFrame(rows); pd.set_option('display.width',250,'display.max_columns',40); R.round(3).to_csv('/workspace/out/q12c/geometry.csv',index=False); print(R.round(3).T.to_string())
