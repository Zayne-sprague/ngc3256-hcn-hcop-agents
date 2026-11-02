# Blind matched-filter line search on every extracted spectrum; match to Splatalogue list at the systemic velocity
import sys,glob; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
Z=float(sys.argv[1]) if len(sys.argv)>1 else ZSYS
LL=pd.read_csv('/workspace/out/linelist_splatalogue.csv')
LL['name']=LL.name.str.replace(r'<[^>]+>','',regex=True).str.strip()
LL['f']=LL.orderedfreq/1e3
def robust_sigma(x):
    m=np.isfinite(x); x=x[m]
    for _ in range(5):
        med=np.median(x); s=1.4826*np.median(np.abs(x-med)); x=x[np.abs(x-med)<3*s]
    return s
def search(nu,y,widths=(3,6,12,24,48,96),thr=5.0):
    y=np.nan_to_num(y); cands=[]
    base=y-np.median(y)
    for w in widths:
        if w>=len(y)//3: continue
        k=np.ones(w)/w; s=np.convolve(base,k,'same')
        # remove local baseline? use global robust noise of smoothed series (iteratively clipped)
        sig=robust_sigma(s[::max(1,w//2)]) if w>1 else robust_sigma(s)
        sig=robust_sigma(s)
        snr=s/sig
        idx=np.where(snr>thr)[0]
        if idx.size==0: continue
        groups=np.split(idx,np.where(np.diff(idx)>w)[0]+1)
        for g in groups:
            j=g[np.argmax(snr[g])]; cands.append((nu[j],w,snr[j],s[j],sig))
    # merge: keep best snr within 1.5*widths
    cands.sort(key=lambda c:-c[2]); keep=[]
    for c in cands:
        if all(abs(c[0]-k[0])>1.5*max(c[1],k[1])*abs(nu[1]-nu[0]) for k in keep): keep.append(c)
    return sorted(keep)
def match(nu_obs,z=Z,tol_kms=150.,maxn=4):
    fr=nu_obs*(1+z)/1e9
    d=np.abs(LL.f-fr)/fr*C
    c=LL[d<tol_kms].copy(); c['dv']=(fr-c.f)/fr*C
    c=c.sort_values('sijmu2',ascending=False).head(maxn)
    return fr,[f"{r['name']} {r.resolved_QNs.strip()} {r.f:.3f} (dv={r.dv:+.0f}, Eu={r.upper_state_energy_K:.0f}K, Sμ²={r.sijmu2:.1f})" for _,r in c.iterrows()]
if __name__=='__main__':
    rows=[]
    for fn in sorted(glob.glob('/workspace/out/spectra/*.npz')):
        z=np.load(fn); nm=os.path.basename(fn)[:-4]
        if 'cont' in nm: continue
        for key in ['N_pix','S_pix','peak_pix','tot_r8']:
            y=z[key]; 
            for nu0,w,snr,amp,sig in search(z['nu'],y):
                fr,m=match(nu0)
                rows.append(dict(file=nm,spec=key,nu_obs=nu0/1e9,nu_rest=fr,width_ch=w,snr=snr,amp_mJy=amp*1e3,sig_mJy=sig*1e3,matches=' | '.join(m[:3])))
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/s3_candidates.csv',index=False)
    pd.set_option('display.width',250,'display.max_colwidth',120,'display.max_rows',500)
    print(df[df.spec.isin(['N_pix','S_pix'])].round(2).to_string())
