# Q14h (fix 2,3): template deblend of the two line clusters, all at matched beams per column.
#  (A) Xb4d spw29 / Xb4f spw20: CH3OH 2(0)-1(-1)E 254.015377 + SiS 14-13 254.103211 (+104 km/s) + 30SiO 6-5 254.216656 (+237 km/s; NB the reviewer quoted 254.76: wrong, CDMS 254.2167, Eu 42.7 K, logA -3.07).
#  (B) Xb47 spw22: CS 7-6 342.883 + H2CS 10(0,10)-9(0,9) 342.946424 (+56 km/s) + 29SiO 8-7 342.980843 (+85 km/s).
# Line profile = the HCO+ 3-2 (A) / HCO+ 4-3 (B) spectrum of the same pixel and beam (data-driven template, baseline removed), alternatives: HCN spectrum, Gaussian(v0,FWHM of HCO+).
# Linear least squares: baseline (linear) + sum_k A_k T(v_k - shift); errors from the residual rms (robust, channel-correlation inflation 1.3) + template spread. Output out/q14e/deblend_*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
from scipy.interpolate import interp1d
INF=1.3
def template(v,y,v0,fw):
    side=(np.abs(v-v0)>2.0*fw)&(np.abs(v-v0)<900); p=np.polyfit(v[side]-v0,y[side],1) if side.sum()>8 else [0,0]
    t=y-np.polyval(p,v-v0); t=np.where(np.abs(v-v0)<2.2*fw,t,0.); dv=abs(np.median(np.diff(v))); t=t/(t.sum()*dv)
    return interp1d(v-v0,t,bounds_error=False,fill_value=0.)
def gtemp(v0,fw):
    s=fw/2.3548; return lambda x: np.exp(-0.5*(x/s)**2)/(s*2.5066)
def fitcluster(fr,y,v0,T,comps,vlo,vhi):
    """comps: list of (name,rest GHz). window: v relative to comps[0] in [vlo,vhi] around v0. returns amplitudes, errors, cov, rms, residual"""
    v1=vel(fr,comps[0][1]); sel=(v1>v0+vlo)&(v1<v0+vhi); f=fr[sel]; yy=y[sel]; dv=abs(np.median(np.diff(v1)))
    cols=[np.ones(sel.sum()),(v1[sel]-v0)/100.]+[T(vel(f,r)-v0) for _,r in comps]; A=np.vstack(cols).T
    use=np.ones(len(yy),bool)
    for _ in range(3):
        p=np.linalg.lstsq(A[use],yy[use],rcond=None)[0]; r=yy-A@p; sg=robust_sigma(r[use]); use=np.abs(r)<4*sg
    cov=np.linalg.pinv(A[use].T@A[use])*sg**2*INF**2
    return p[2:],np.sqrt(np.diag(cov)[2:]),cov[2:,2:],sg,len(yy)
CL={'A':dict(comps=[('CH3OH',254.015377),('SiS',254.103211),('30SiO',254.216656)],ref='HCO+ 3-2',alt='HCN 3-2',f0={'HCO+ 3-2':267.557626,'HCN 3-2':265.886434},
   spw={'Xb4d':{'cl':'spw29','HCO+ 3-2':'spw25','HCN 3-2':'spw27'},'Xb4f':{'cl':'spw20','HCO+ 3-2':'spw16','HCN 3-2':'spw18'}},vlo=-380,vhi=+800),
    'B':dict(comps=[('CS',342.883000),('H2CS',342.946424),('29SiO',342.980843)],ref='HCO+ 4-3',alt='HCN 4-3',f0={'HCO+ 4-3':356.734223,'HCN 4-3':354.505478},
   spw={'Xb47':{'cl':'spw22','HCO+ 4-3':'spw18','HCN 4-3':'spw16'}},vlo=-380,vhi=+800)}
JOBS=[('A','Xb4d','1.7',(1.7,1.7,0.)),('A','Xb4d','4.45x2.95',(4.45,2.95,87.)),('A','Xb4d','6.3x3.8',(6.3,3.8,82.)),('A','Xb4f','6.3x3.8',(6.3,3.8,82.)),
      ('B','Xb47','4.45x2.95',(4.45,2.95,87.)),('B','Xb47','6.3x3.8',(6.3,3.8,82.))]
RG=['N','S','Bridge','Wclump','Soutflow']; rows=[]
for cl,ous,bn,tb in JOBS:
    P=CL[cl]; D={}
    for key in ('cl',P['ref'],P['alt']):
        f0=P['f0'].get(key,None); fr,a,cw,pix,beam,nb=cut(ous,P['spw'][ous][key],tb,half=20.); D[key]=(fr,a,cw)
    for reg in RG:
        ra,de=REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']
        fr_r,a_r,cw_r=D[P['ref']]; v=vel(fr_r,P['f0'][P['ref']]); s_r=sample(a_r,cw_r,ra,de); p=hco_fit(v,s_r); v0,fw=p[1],p[2]
        fr_a,a_a,cw_a=D[P['alt']]; va=vel(fr_a,P['f0'][P['alt']]); s_a=sample(a_a,cw_a,ra,de)
        frc,ac,cwc=D['cl']; yc=sample(ac,cwc,ra,de)
        res={}; comps_used={}
        for tn,T in (('HCO+',template(v,s_r,v0,fw)),('alt',template(va,s_a,v0,fw)),('gauss',gtemp(v0,fw))):
            for fixed in ('all','no_blend','noH2CS','noSiO'):
                if cl=='A' and fixed in('noH2CS','noSiO'): continue
                comps={'all':P['comps'],'no_blend':P['comps'][:1],'noH2CS':[P['comps'][0],P['comps'][2]] if cl=='B' else None,'noSiO':P['comps'][:2]}[fixed]
                try: A,eA,cov,sg,n=fitcluster(frc,yc,v0,T,comps,P['vlo'],P['vhi'])
                except Exception as e: continue
                res[(tn,fixed)]=(A,eA,cov,sg,n); comps_used[(tn,fixed)]=comps
        for (tn,fixed),(A,eA,cov,sg,n) in res.items():
            for k,(nm,_) in enumerate(comps_used[(tn,fixed)]):
                rows.append(dict(cluster=cl,ous=ous,beam=bn,region=reg,template=tn,model=fixed,comp=nm,I=A[k],eI=eA[k],snr=A[k]/eA[k],v0=v0,fwhm=fw,rms_K=sg,nchan=n,
                                 corr_with_next=(cov[k,k+1]/np.sqrt(cov[k,k]*cov[k+1,k+1]) if k+1<len(A) else np.nan)))
    print(cl,ous,bn,flush=True)
R=pd.DataFrame(rows); R.to_csv('/workspace/out/q14e/deblend_all.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',400)
print(R[(R.template=='HCO+')&(R.model=='all')][['cluster','ous','beam','region','comp','I','eI','snr','corr_with_next']].round(3).to_string())
