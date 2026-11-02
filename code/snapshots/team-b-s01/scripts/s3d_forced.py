# Forced photometry on every catalogue line (tier 1-2) in the extracted spectra: W over +-150 km/s of VSYS, empirical noise from line-free windows
import sys,glob; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
from s3c_scan import LL
HALF=150.   # km/s half window (HCN/HCO+ FWZI ~ 400 km/s, FWHM~150-300) - narrower than full profile -> conservative for weak lines
def wsum(nu,y,f0,half=HALF):
    v=C*(1-nu/f0); dv=abs(np.median(np.diff(v)))
    m=abs(v)<=half; n=m.sum()
    if n<0.8*(2*half/dv) : return None   # window not fully covered
    # local baseline: |v| in [half+50, half+700] with iterative clipping
    b=(abs(v)>half+50)&(abs(v)<half+800)
    if b.sum()<20: return None
    yb=y[b]; 
    for _ in range(4):
        med=np.median(yb); s=1.4826*np.median(np.abs(yb-med)); yb=yb[abs(yb-med)<3*s]
    base=np.median(yb); rms=np.std(yb)
    W=np.nansum(y[m]-base)*dv
    # empirical sigma of an n-channel boxcar sum: scan boxes within baseline-region (clipped)
    yy=y-base; ok=(abs(v)>half+50)&(abs(v)<half+800)&(abs(yy)<4*rms)
    sums=[np.sum(yy[i:i+n]) for i in range(0,len(yy)-n) if ok[i:i+n].all()]
    sigW=(np.std(sums) if len(sums)>5 else rms*np.sqrt(n))*dv
    return W,sigW,rms,n
rows=[]
cat=LL[(LL.tier<=2)&(LL.upper_state_energy_K<600)].copy()
for fn in sorted(glob.glob('/workspace/out/spectra/*.npz')):
    nm=os.path.basename(fn)[:-4]
    if 'cont' in nm or 'repBW' in nm or 'tapered' in nm or 'manual' in nm: continue
    z=np.load(fn); nu=z['nu']; fr_lo,fr_hi=rest_of(nu.min())/1e9,rest_of(nu.max())/1e9
    c=cat[(cat.freq_GHz>fr_lo)&(cat.freq_GHz<fr_hi)]
    # merge hyperfine/duplicate entries within 1 MHz
    c=c.sort_values('proxy',ascending=False).drop_duplicates(subset=['sp',c.freq_GHz.round(3).name]) if False else c
    for key in ['N_ap','S_ap','tot_r8']:
        for _,r in c.iterrows():
            f0=r.freq_GHz*1e9*(1-VSYS/C)   # observed freq at systemic
            o=wsum(nu,z[key],f0)
            if o is None: continue
            W,sW,rms,n=o
            if W/sW>=3.0: rows.append(dict(file=nm,spec=key,line=r['name'],qn=str(r.resolved_QNs).strip()[:25],rest_GHz=r.freq_GHz,Eu=r.upper_state_energy_K,tier=r.tier,proxy=r.proxy,W_Jykms=W,sigW=sW,snr=W/sW,rms_mJy=rms*1e3))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s3d_forced.csv',index=False); print(len(df))
