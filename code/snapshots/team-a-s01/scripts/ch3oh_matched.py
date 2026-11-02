# Matched-beam (9.36x6.28" PA -88, the Band-5 beam) CH3OH comparison at the Band-5 peak pixel (phase centre): Band 6 Xb4f and Xb4d cubes convolved to that beam.
# LTE ladder fit (A,Trot) over Q-branch groups J1-3, J4, J5, J6, J7, J8 (Band 5) + 2(0)-1(-1)E, 5(-2,4)-4(-1,3)E, 6(-1,5)-5(-2,4)E (Band 6). Output out/ch3oh_matched.csv/.png
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from matched import smoothed_spectrum
from common import f as fpath, freqaxis
from ch3oh import c,toK,C
from astropy.io import fits
RA,DE=156.9634583333,-43.90461111111; TB=(9.357,6.284,-88.45); V0,FW=2777.,159.
def ltew(f0,T): 
    l=c[(c.f0-f0).abs()<2e-4].iloc[0]; return l.f0*l.sijmu2*np.exp(-l.Eu/T),l.Eu,l.sijmu2
def fit_one(fr,y_K,fr0s,half=450.,T=10.,v0=V0,fw=FW,order=1):
    fc=fr0s[0]*(1-v0/C); sel=np.abs((1-fr/fc)*C)<half+60; f=fr[sel]; y=y_K[sel]
    ws=np.array([ltew(a,T)[0] for a in fr0s]); ws/=ws.sum(); t=np.zeros_like(f)
    for fi,wi in zip(fr0s,ws):
        fl=fi*(1-v0/C); sg=fw/2.3548/C*fl; t+=wi*np.exp(-0.5*((f-fl)/sg)**2)/(fw/2.3548*2.5066)
    x=(f-fc)/fc*C/100.; M=np.vstack([t]+[x**k for k in range(order+1)]).T; use=np.ones(len(f),bool)
    for _ in range(4):
        p=np.linalg.lstsq(M[use],y[use],rcond=None)[0]; r=y-M@p; far=np.abs((1-f/fc)*C-(v0-v0))>fw
        sg=1.4826*np.median(np.abs(r[far&use])); use=np.abs(r)<4*sg
    cv=np.linalg.pinv(M[use].T@M[use])*sg**2
    return p[0],np.sqrt(cv[0,0]),sg
GROUPS=[('Q J=1-3',(165.050175,165.061130,165.099240),'Xb55','spw22'),('Q J=4',(165.190475,),'Xb55','spw22'),('Q J=5',(165.369341,),'Xb55','spw22'),
 ('Q J=6',(165.678649,),'Xb55','spw22'),('Q J=7',(166.169098,),'Xb55','spw20'),('Q J=8',(166.898566,),'Xb55','spw20'),
 ('2(0)-1(-1) E',(254.015377,),'Xb4f','spw20'),('2(0)-1(-1) E',(254.015377,),'Xb4d','spw29'),
 ('5(-2,4)-4(-1,3) E',(266.838148,),'Xb4f','spw18'),('5(-2,4)-4(-1,3) E',(266.838148,),'Xb4d','spw27'),('6(-1,5)-5(-2,4) E',(265.289562,),'Xb4f','spw18')]
if __name__=='__main__':
    rows=[]
    for name,fr0s,ous,spw in GROUPS:
        p=fpath(ous,spw); hd=fits.getheader(p); fr=freqaxis(hd)/1e9
        if ous=='Xb55':
            h=fits.open(p,memmap=True); sp=np.nan_to_num(np.array(h[0].data[0][:,int(round(float(__import__("astropy.wcs",fromlist=["WCS"]).WCS(hd).celestial.world_to_pixel_values(RA,DE)[1]))),int(round(float(__import__("astropy.wcs",fromlist=["WCS"]).WCS(hd).celestial.world_to_pixel_values(RA,DE)[0])))],float))
            bm=(hd['BMAJ']*3600,hd['BMIN']*3600)
        else:
            sp,_,nb=smoothed_spectrum(p,TB,RA,DE); bm=TB[:2]
        yK=toK(sp,fr,bm)
        for T in (8.,):
            A,eA,sg=fit_one(fr,yK,fr0s,T=T)
            w,Eu,smu=ltew(fr0s[0],T)
            rows.append(dict(line=name,ous=ous,spw=spw,nu_rest=fr0s[0],Eu=Eu,Smu2=smu,I_Kkms=A,eI_Kkms=eA,snr=A/eA,sigma_ch_K=sg,nblend=len(fr0s),beam='9.36x6.28 matched' if ous!='Xb55' else 'native 9.36x6.28'))
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/ch3oh_matched.csv',index=False)
    pd.set_option('display.width',200); print(df.round(3).to_string())
