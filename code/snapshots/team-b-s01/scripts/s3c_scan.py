# Template (boxcar) line scan in velocity space on all extracted spectra; identification with Splatalogue list at VSYS.
import sys,glob; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
LL=pd.read_csv('/workspace/out/linelist.csv')
TIER1={'CO','13CO','C18O','HCN','H13CN','HC15N','HCO+','H13CO+','HC18O+','HNC','HN13C','CS','13CS','C34S','SO','34SO','SiO','29SiO','30SiO','CN','13CN','CH3OH','HC3N','CCH','H2CO','H2CS','c-HCCCH','CH3CCH','N2H+','H2S','HOC+','NS','H-RRL','OCS','HNCO','SO2','SiS','HCS+','NO','HCO','CH3CN','13CCH','C13CH','H213CO'}
def vib(n): return (' v' in n) and ('v=0' not in n and 'v = 0' not in n)
LL['tier']=np.where(LL.sp.isin(TIER1)&~LL.name.map(vib),1,np.where(LL.sp.isin(TIER1),2,3))
LL['proxy']=LL.sijmu2.fillna(10)*np.exp(-LL.upper_state_energy_K/40.)   # crude LTE-ish intensity proxy (Tex=40 K), only for ranking
def robust_sigma(x):
    x=x[np.isfinite(x)]
    for _ in range(6):
        med=np.median(x); s=1.4826*np.median(np.abs(x-med)); x=x[np.abs(x-med)<3.5*s]
    return s
def scan(nu,y,f0_guess,wv,thr):
    """boxcar of full width wv km/s; returns candidates (peak freq, snr, integrated, err)"""
    v=C*(1-nu/f0_guess); dv=abs(np.median(np.diff(v))); w=max(1,int(round(wv/dv)))
    y=np.nan_to_num(y); y=y-np.median(y)
    k=np.ones(w)/w; s=np.convolve(y,k,'same'); sig=robust_sigma(s)
    snr=s/sig; idx=np.where(snr>thr)[0]; idx=idx[(idx>w//2)&(idx<len(y)-w//2-1)]
    if idx.size==0: return []
    groups=np.split(idx,np.where(np.diff(idx)>w)[0]+1); res=[]
    for g in groups:
        j=g[np.argmax(snr[g])]; res.append((nu[j],snr[j],s[j],sig,w,dv))
    return res
def identify(nu_obs,tol=60.):
    fr=rest_of(nu_obs)/1e9
    dvs=(fr-LL.freq_GHz)/fr*C
    c=LL[np.abs(dvs)<tol].copy(); c['dv']=dvs[np.abs(dvs)<tol]
    c=c.sort_values(['tier','proxy'],ascending=[True,False])
    return fr,c
if __name__=='__main__':
    thr=float(sys.argv[1]) if len(sys.argv)>1 else 5.0
    rows=[]
    for fn in sorted(glob.glob('/workspace/out/spectra/*.npz')):
        z=np.load(fn); nm=os.path.basename(fn)[:-4]
        if 'cont' in nm or 'repBW' in nm or 'manual.image.' in nm and nm.endswith('pbcor') is False: pass
        if 'cont' in nm: continue
        nu=z['nu']; f0=float(np.mean(nu))*1.0093
        for key in ['N_pix','S_pix','peak_pix','tot_r8']:
            for wv in (300.,100.):
                for nu0,snr,amp,sig,w,dv in scan(nu,z[key],f0,wv,thr):
                    fr,c=identify(nu0)
                    ids=[f"{r['name']} {str(r.resolved_QNs).strip()[:22]} {r.freq_GHz:.3f}({r.dv:+.0f})" for _,r in c.head(3).iterrows()]
                    rows.append(dict(file=nm,spec=key,box_kms=wv,nu_obs=nu0/1e9,nu_rest=fr,snr=snr,mean_mJy=amp*1e3,sig_mJy=sig*1e3,wint_Jykms=amp*(w*dv),nid=len(c),best=' | '.join(ids)))
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/s3c_scan.csv',index=False); print(len(df))
