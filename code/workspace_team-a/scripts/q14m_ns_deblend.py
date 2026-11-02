# Q14m (fix 2): S-N contrast of CH3OH 254/HCO+ 3-2 with the template-deblended (CH3OH+SiS+30SiO, HCO+ line-shape) whole-line amplitudes; null = same fit at 80 random off-source position pairs (9-16" ring; template = HCO+ spectrum of N / S),
# plus the paired shifted-velocity null is not possible for a template fit (lines everywhere), so the ring null (baseline/sidelobe/primary-beam structure + thermal noise) is the empirical error; error floor = fit error x1.3.
# Output out/q14e/ns_deblend.csv. python3 -I q14m_ns_deblend.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
from scipy.interpolate import interp1d
INF=1.3
def template(v,y,v0,fw):
    side=(np.abs(v-v0)>2.0*fw)&(np.abs(v-v0)<900); p=np.polyfit(v[side]-v0,y[side],1); t=y-np.polyval(p,v-v0); t=np.where(np.abs(v-v0)<2.2*fw,t,0.); dv=abs(np.median(np.diff(v))); return interp1d(v-v0,t/(t.sum()*dv),bounds_error=False,fill_value=0.)
COMPS=[254.015377,254.103211,254.216656]
def amp(fr,y,v0,T,vlo=-380,vhi=800):
    v1=vel(fr,COMPS[0]); sel=(v1>v0+vlo)&(v1<v0+vhi); f=fr[sel]; yy=y[sel]
    A=np.vstack([np.ones(sel.sum()),(v1[sel]-v0)/100.]+[T(vel(f,r)-v0) for r in COMPS]).T; use=np.ones(len(yy),bool)
    for _ in range(3):
        p=np.linalg.lstsq(A[use],yy[use],rcond=None)[0]; r=yy-A@p; sg=robust_sigma(r[use]); use=np.abs(r)<4*sg
    cov=np.linalg.pinv(A[use].T@A[use])*sg**2*INF**2; return p[2],np.sqrt(cov[2,2])
SPW={'Xb4d':('spw25','spw29'),'Xb4f':('spw16','spw20')}
BEAMS=[('Xb4d','1.7',(1.7,1.7,0.)),('Xb4d','3.0',(3.0,3.0,0.)),('Xb4d','4.45x2.95',(4.45,2.95,87.)),('Xb4d','6.3x3.8',(6.3,3.8,82.)),('Xb4f','6.3x3.8',(6.3,3.8,82.))]
rng=np.random.default_rng(21); out=[]
def P(r,t): return (CEN[0]+r*np.sin(t)/3600/np.cos(np.deg2rad(CEN[1])),CEN[1]+r*np.cos(t)/3600)
for ous,bn,tb in BEAMS:
    fh,ah,cwh,*_=cut(ous,SPW[ous][0],tb,half=20.); fc,ac,cwc,*_=cut(ous,SPW[ous][1],tb,half=20.); vh=vel(fh,267.557626)
    T={};H={};v0={};
    for reg in ('N','S'):
        s=sample(ah,cwh,REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']); p=hco_fit(vh,s); v0[reg]=p[1]; T[reg]=template(vh,s,p[1],p[2]); H[reg]=p[0]
    xN,eN=amp(fc,sample(ac,cwc,REG.loc['N','RA_deg'],REG.loc['N','Dec_deg']),v0['N'],T['N']); xS,eS=amp(fc,sample(ac,cwc,REG.loc['S','RA_deg'],REG.loc['S','Dec_deg']),v0['S'],T['S'])
    dR=xS/H['S']-xN/H['N']; nulN=[];nulS=[]
    for k in range(80):
        for reg,lst in (('N',nulN),('S',nulS)):
            r,t=rng.uniform(9,16),rng.uniform(0,2*np.pi); lst.append(amp(fc,sample(ac,cwc,*P(r,t)),v0[reg],T[reg])[0])
    nulN=np.array(nulN);nulS=np.array(nulS); nd=nulS/H['S']-nulN[::-1]/H['N']
    sN,sS=max(np.std(nulN),eN),max(np.std(nulS),eS); edR=np.hypot(sS/H['S'],sN/H['N'])
    out.append(dict(ous=ous,beam=bn,A_N=xN,eA_N=eN,A_S=xS,eA_S=eS,H_N=H['N'],H_S=H['S'],R_N=xN/H['N'],R_S=xS/H['S'],dR=dR,null_rms_N=np.std(nulN),null_rms_S=np.std(nulS),null_mean_N=nulN.mean(),null_mean_S=nulS.mean(),
                    dR_err=edR,sig_dR=dR/edR,sig_dR_ringpairs=dR/np.std(nd),S_det_sig=xS/sS,N_det_sig=xN/sN,S_over_N=xS/H['S']/(xN/H['N']) if xN>0 else np.nan))
    print(ous,bn,flush=True)
O=pd.DataFrame(out); O.to_csv('/workspace/out/q14e/ns_deblend.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_columns',30); print(O.round(4).to_string())
