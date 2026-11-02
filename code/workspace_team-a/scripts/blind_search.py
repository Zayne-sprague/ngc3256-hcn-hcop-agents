# Blind line search: boxcar matched filter (FWHM-ish 150 km/s) on N/S spectra; robust noise (MAD of filtered spectrum, iterative clip).
import sys,warnings; warnings.filterwarnings('ignore')
import numpy as np, glob, os, json
from scipy.ndimage import uniform_filter1d
C=299792.458
def search(fr,s,wkms=150.,thr=4.0):
    dv=abs(np.median(np.diff(fr)))/np.mean(fr)*C
    nb=max(1,int(round(wkms/dv)))
    ok=np.isfinite(s); s=np.where(ok,s,0)
    sm=uniform_filter1d(s,nb,mode='nearest')   # mean over nb channels
    # noise of smoothed: robust, iterative sigma clip
    m=np.ones_like(sm,bool)
    for _ in range(8):
        sg=1.4826*np.median(np.abs(sm[m]-np.median(sm[m]))); m=np.abs(sm-np.median(sm[m]))<2.5*sg
    sg=np.std(sm[m]); base=np.median(sm[m])
    snr=(sm-base)/sg
    pk=[i for i in range(1,len(snr)-1) if snr[i]>thr and snr[i]>=snr[i-1] and snr[i]>=snr[i+1]]
    # merge peaks within 1.5 boxes
    out=[]
    for i in pk:
        if out and abs(i-out[-1])<1.5*nb:
            if snr[i]>snr[out[-1]]: out[-1]=i
        else: out.append(i)
    return [(fr[i]/1e9,float(sm[i]*1e3),float(snr[i])) for i in out],sg*1e3,dv,nb
if __name__=='__main__':
    res={}
    for k in sorted(glob.glob('/workspace/out/spec_*.npz')):
        d=np.load(k); key=os.path.basename(k)[5:-4]
        for n in [x for x in d.files if x!='freq']:
            pk,sg,dv,nb=search(d['freq'],d[n])
            res[f'{key}|{n}']=dict(sigma_smoothed_mJy=sg,dv=dv,nbox=nb,peaks=pk)
            print(key,n,f'sig(box{nb})={sg:.2f}mJy',' '.join(f'{f:.3f}GHz({s:.1f}mJy,{q:.0f}s)' for f,s,q in pk))
    json.dump(res,open('/workspace/out/blind_peaks.json','w'),indent=1)
