# Q9b fix 3: 4-3/3-2 across the line profile at the common 6.5" beam (Xb47 4-3, Xb4f 3-2; merged N+S), un-masked 40 km/s bins and windows,
# plus a forward model: N/S weights per channel from the Xb4d 3-2 cubes (convolved to 6.5" about the pointing centre, N = y>0, S = y<0) applied to the
# observed Xb4f 3-2 spectrum, I43(v)=I32obs(v)*[R_N f_N(v)+R_S (1-f_N(v))] with constant R_N,R_S (Tb ratios) -> does a constant-per-nucleus ratio reproduce the red-wing drop?
# python3 -I q9b_ladder.py
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
from matched import smoothed_spectrum
from q9b_lib import cube,RA0 as _r,DE0 as _d
import pandas as pd
from scipy.optimize import least_squares
TB=(6.5,6.5,0.); ra,de=RA0,DE0
LINES={'HCN':{'32':('Xb4f','spw18',265.886434),'43':('Xb47','spw16',354.505477)},'HCO+':{'32':('Xb4f','spw16',267.557626),'43':('Xb47','spw18',356.734223)}}
spec={}
for mol,d in LINES.items():
    for J,(ous,spw,f0) in d.items():
        p=f(ous,spw,kind='cube.I.pbcor'); sm,nat,nb=smoothed_spectrum(p,TB,ra,de); hd=fits.open(p,memmap=True)[0].header
        v=(1-freqaxis(hd)/1e9/f0)*C-VSYS; m=(np.abs(v)>600)&(np.abs(v)<2000); pf=np.polyfit(v[m],sm[m],1); y=sm-np.polyval(pf,v)
        rms=1.4826*np.median(np.abs(y[m]-np.median(y[m]))); spec[(mol,J)]=(v,y,rms,f0)
# Xb4d N/S-weighted 3-2 spectra convolved to 6.5" about the pointing centre
def nsweights(mol):
    c=cube('Xb4d',mol); r2=c.X**2+c.Y**2; G=np.exp(-4*np.log(2)*r2/6.5**2)
    fl=c.cube/(OM(c.beam)/c.pix**2)           # Jy per pixel per channel
    PN=(fl*(G*(c.Y>0))[None]).sum((1,2)); PS=(fl*(G*(c.Y<=0))[None]).sum((1,2)); return c.v,PN,PS
NB=40.; edges=np.arange(-450,451,NB)
def binsum(v,y,lo,hi):
    dv=abs(np.median(np.diff(v))); w=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,None)/dv; return (y*w).sum()*dv
out=[]; fits_=[]
for mol in LINES:
    v43,y43,r43,f43=spec[(mol,'43')]; v32,y32,r32,f32=spec[(mol,'32')]; fac=(f32/f43)**2
    d43=abs(np.median(np.diff(v43))); d32=abs(np.median(np.diff(v32)))
    vd,PN,PS=nsweights(mol)
    cen=(edges[:-1]+edges[1:])/2
    A=np.array([binsum(v43,y43,a,b) for a,b in zip(edges[:-1],edges[1:])]); B=np.array([binsum(v32,y32,a,b) for a,b in zip(edges[:-1],edges[1:])])
    eA=np.full(len(cen),r43*d43*np.sqrt(NB/d43)); eB=np.full(len(cen),r32*d32*np.sqrt(NB/d32))
    pn=np.array([binsum(vd,PN,a,b) for a,b in zip(edges[:-1],edges[1:])]); ps=np.array([binsum(vd,PS,a,b) for a,b in zip(edges[:-1],edges[1:])])
    fN=pn/(pn+ps)
    for ci,a,ea,b,eb,fn,p_,s_ in zip(cen,A,eA,B,eB,fN,pn,ps):
        R=a/b*fac if b!=0 else np.nan; out.append(dict(mol=mol,v=ci,I43=a,e43=ea,I32=b,e32=eb,R_Tb=R,eR=abs(R)*np.hypot(ea/a,eb/b) if a*b!=0 else np.nan,fN_Xb4d=fn,PN=p_,PS=s_,snr43=a/ea,snr32=b/eb))
    g=pd.DataFrame([o for o in out if o['mol']==mol]); ok=(np.abs(g.v)<=150)&(g.snr43>3)&(g.snr32>3)
    # fit: I43*fac' = I32obs*(R_N fN+R_S(1-fN)) in Tb, fit over bins with I32 S/N>3 (|v|<=250), no mask on 4-3 (linear, errors propagated)
    sel=(g.snr32>3)&(np.abs(g.v)<=250)
    X=np.vstack([g.I32[sel]*g.fN_Xb4d[sel],g.I32[sel]*(1-g.fN_Xb4d[sel])]).T/fac; yv=g.I43[sel].values; w=1/g.e43[sel].values
    # include the 3-2 error propagated (small): use e43 only plus 3% of model
    mod_err=np.sqrt(g.e43[sel].values**2+(0.0*yv)**2)
    P,res,_,_=np.linalg.lstsq(X/mod_err[:,None],yv/mod_err,rcond=None); cov=np.linalg.inv((X/mod_err[:,None]).T@(X/mod_err[:,None]))
    chi=(((yv-X@P)/mod_err)**2).sum(); dof=sel.sum()-2
    c1=np.linalg.lstsq((X.sum(1))[:,None]/mod_err[:,None],yv/mod_err,rcond=None)[0][0]; chi1=(((yv-X.sum(1)*c1)/mod_err)**2).sum()
    print(f'\n{mol}: 4-3/3-2 Tb per 40 km/s bin at 6.5" (no S/N mask; fN = N fraction of the Xb4d 3-2 flux in the 6.5" beam)')
    print(g[['v','I43','e43','I32','e32','R_Tb','eR','fN_Xb4d','snr43','snr32']].round(3).to_string(index=False))
    print(f'  forward model (bins with 3-2 S/N>3, |v|<=250): R_N={P[0]:.3f}+-{np.sqrt(cov[0,0]):.3f}  R_S={P[1]:.3f}+-{np.sqrt(cov[1,1]):.3f}  chi2/dof={chi:.1f}/{dof}; single constant R: {c1:.3f} chi2/dof={chi1:.1f}/{dof+1}')
    gg=g[sel].copy(); gg['model']=(X@P)*fac/gg.I32; gg['model_single']=X.sum(1)*c1*fac/gg.I32
    print(gg[['v','R_Tb','eR','model','model_single','fN_Xb4d']].round(3).to_string(index=False))
    fits_.append(dict(mol=mol,R_N=P[0],eR_N=np.sqrt(cov[0,0]),R_S=P[1],eR_S=np.sqrt(cov[1,1]),chi2=chi,dof=dof,R_single=c1,chi2_single=chi1))
    # three-component variant: + extended remainder E(v)=I32obs-(PN+PS) (Xb4d-missing flux; ratio R_E free), bins with 3-2 S/N>10
    s3=(g.snr32>10)&(np.abs(g.v)<=250); E=(g.I32-g.PN-g.PS)[s3].values
    X3=np.vstack([g.PN[s3],g.PS[s3],E]).T/fac; y3=g.I43[s3].values; e3=g.e43[s3].values
    P3=np.linalg.lstsq(X3/e3[:,None],y3/e3,rcond=None)[0]; cv3=np.linalg.inv((X3/e3[:,None]).T@(X3/e3[:,None])); ch3=(((y3-X3@P3)/e3)**2).sum()
    print(f'  3-component (N,S from Xb4d + extended remainder; bins 3-2 S/N>10): R_N={P3[0]:.3f}+-{np.sqrt(cv3[0,0]):.3f} R_S={P3[1]:.3f}+-{np.sqrt(cv3[1,1]):.3f} R_E={P3[2]:.3f}+-{np.sqrt(cv3[2,2]):.3f} chi2/dof={ch3:.1f}/{s3.sum()-3}')
    fits_.append(dict(mol=mol,model='3comp',R_N=P3[0],eR_N=np.sqrt(cv3[0,0]),R_S=P3[1],eR_S=np.sqrt(cv3[1,1]),R_E=P3[2],eR_E=np.sqrt(cv3[2,2]),chi2=ch3,dof=s3.sum()-3))
    gg3=g[s3].copy(); gg3['model3']=(X3@P3)*fac/gg3.I32; print(gg3[['v','R_Tb','eR','model3']].round(3).to_string(index=False))
    # window-integrated ladder ratios (unmasked)
    for lo,hi in ((-250,-100),(-100,0),(0,100),(-100,100),(100,250),(-250,250)):
        a=binsum(v43,y43,lo,hi); b=binsum(v32,y32,lo,hi); ea=r43*d43*np.sqrt((hi-lo)/d43); eb=r32*d32*np.sqrt((hi-lo)/d32); R=a/b*fac
        pnw=binsum(vd,PN,lo,hi); psw=binsum(vd,PS,lo,hi)
        print(f'  window {lo}..{hi}: I43={a:.3f}+-{ea:.3f} I32={b:.3f}+-{eb:.3f} -> R_Tb={R:.3f}+-{abs(R)*np.hypot(ea/a,eb/b):.3f}; N fraction of 3-2 (Xb4d) {pnw/(pnw+psw):.2f}; model (R_N,R_S fit) {(P[0]*pnw+P[1]*psw)/(pnw+psw):.3f}')
pd.DataFrame(out).to_csv('/workspace/out/q9b/ladder_bins.csv',index=False); pd.DataFrame(fits_).to_csv('/workspace/out/q9b/ladder_fit.csv',index=False)
# recovery check: Xb4d-based 3-2 spectrum (convolved to 6.5") vs Xb4f observed, |v|<100
for mol in LINES:
    v32,y32,r32,f32=spec[(mol,'32')]; vd,PN,PS=nsweights(mol)
    a=binsum(vd,PN+PS,-100,100); b=binsum(v32,y32,-100,100); print(f'{mol} 3-2 |v|<100: Xb4d->6.5" {a:.3f} vs Xb4f {b:.3f} Jy km/s (Xb4d/Xb4f={a/b:.2f})')
