# Q9b: velocity-resolved 4-3/3-2 ratio per molecule at a common 6.5" circular beam (Xb47 vs Xb4f), spectrum at the N-S midpoint (merged nuclei).
# T_b ratio = (S43/S32)*(nu32/nu43)^2 ; 40-km/s bins; noise from line-free channels (|v|>600).  python3 -I q9_ladder.py
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
from matched import smoothed_spectrum
import pandas as pd
TB=(6.5,6.5,0.)
ra=RA0; de=DE0-0.0/3600
LINES={'HCN':{'32':('Xb4f','spw18',265.886434),'43':('Xb47','spw16',354.505477)},
       'HCO+':{'32':('Xb4f','spw16',267.557626),'43':('Xb47','spw18',356.734223)}}
NB=9; out=[]
spec={}
for mol,d in LINES.items():
    for J,(ous,spw,f0) in d.items():
        p=f(ous,spw,kind='cube.I.pbcor'); sm,nat,nb=smoothed_spectrum(p,TB,ra,de); hd=fits.open(p,memmap=True)[0].header
        v=(1-freqaxis(hd)/1e9/f0)*C-VSYS
        # baseline: linear fit to |v|>600 within +-2000
        m=(np.abs(v)>600)&(np.abs(v)<2000); pf=np.polyfit(v[m],sm[m],1); y=sm-np.polyval(pf,v)
        rms=1.4826*np.median(np.abs(y[m]-np.median(y[m])))
        spec[(mol,J)]=(v,y,rms,f0); print(mol,J,ous,spw,'chan rms mJy/beam %.2f'%(rms*1e3),'native beam',nb)
def binned(v,y,rms):
    dv=abs(np.median(np.diff(v))); edges=np.arange(-450,451,40)
    c=[];S=[];E=[]
    for a,b in zip(edges[:-1],edges[1:]):
        k=(v>=a)&(v<b); c.append((a+b)/2); S.append(y[k].sum()*dv/ (b-a)*(b-a)); E.append(rms*dv*np.sqrt(k.sum()))
    return np.array(c),np.array(S),np.array(E)
for mol in LINES:
    c,A,eA=binned(*spec[(mol,'43')][:3]); _,B,eB=binned(*spec[(mol,'32')][:3])
    fac=(spec[(mol,'32')][3]/spec[(mol,'43')][3])**2
    for ci,a,ea,b,eb in zip(c,A,eA,B,eB):
        ok=(a>3*ea)&(b>3*eb); R=a/b*fac
        out.append(dict(mol=mol,v=ci,I43=a,e43=ea,I32=b,e32=eb,R_Tb=R if ok else np.nan,eR=abs(R)*np.hypot(ea/a,eb/b) if ok else np.nan))
    g=pd.DataFrame([o for o in out if o['mol']==mol]); gg=g.dropna()
    w=1/gg.eR**2; Rm=(gg.R_Tb*w).sum()/w.sum(); chi=((gg.R_Tb-Rm)**2*w).sum()
    print('\n',mol,'4-3/3-2 (Tb) per 40 km/s bin at 6.5" beam'); print(g.round(3).to_string(index=False))
    print(f'  mean {Rm:.3f}  chi2/dof={chi:.1f}/{len(gg)-1}')
    # integrated +-250
    for lo,hi in ((-250,250),(-250,-100),(-100,100),(100,250)):
        k=lambda v:(v>=lo)&(v<hi)
        ia=spec[(mol,'43')]; ib=spec[(mol,'32')]
        a=ia[1][k(ia[0])].sum()*abs(np.median(np.diff(ia[0]))); b=ib[1][k(ib[0])].sum()*abs(np.median(np.diff(ib[0])))
        ea=ia[2]*abs(np.median(np.diff(ia[0])))*np.sqrt(k(ia[0]).sum()); eb=ib[2]*abs(np.median(np.diff(ib[0])))*np.sqrt(k(ib[0]).sum())
        R=a/b*fac; print(f'  window {lo}..{hi}: I43={a:.3f}+-{ea:.3f} I32={b:.3f}+-{eb:.3f} Jy/beam km/s  R_Tb={R:.3f}+-{abs(R)*np.hypot(ea/a,eb/b):.3f}')
pd.DataFrame(out).to_csv('/workspace/out/q9/q9_ladder_vel.csv',index=False)
