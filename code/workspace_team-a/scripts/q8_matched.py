# Q8-1: matched-beam T43/T32 at N, S, Bridge, no deblend: moment-0 (|v-VSYS|<350 km/s) maps smoothed EXACTLY (beam-difference kernels) to common beams,
# sampled at the N/S/Bridge positions (q1_deblend POS) and at the local peak within 0.8". 3-2 from three independent arrays (Xb4b, Xb4d, Xb4f) -> also a recovery test.
# Tb ratio = Jy ratio * (nu32/nu43)^2 (same beam). Errors: smoothed-null-window rms at the position + 10% flux scale per OUS (quadrature).
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15)}
BEAMS={'B43(4.4x3.0,87)':(4.4,3.0,87.),'Xb4f(6.2x3.7,82)':(6.2,3.7,82.),'circ6.5':(6.5,6.5,0.),'Xb55(9.5x6.4,-88)':(9.5,6.4,-88.)}
FR={'HCN':(265.886434,354.505478),'HCO+':(267.557626,356.734223)}
def smooth_to(d,tb):
    nb=d['beam']; k=kernel(d['hd'],tb,nb)
    if k is None: return None
    return fftconvolve(np.nan_to_num(d['m0']),k,mode='same')*OM(tb)/OM(nb),k
def noise(d,tb,sm):
    # robust sigma of smoothed null maps within 12" of centre (per position-independent)
    out=[]
    for _,n in d['nulls']:
        nb=d['beam']; k=kernel(d['hd'],tb,nb); s=fftconvolve(np.nan_to_num(n),k,mode='same')*OM(tb)/OM(nb)
        r=np.hypot(d['dx'],d['dy']); v=s[(r<12)&(s!=0)]; out.append(1.4826*np.median(np.abs(v-np.median(v))))
    return np.percentile(out,25)
rows=[]
for sp in ('HCN','HCO+'):
    keys={'3-2 Xb4b':f'{sp} 3-2@Xb4b','3-2 Xb4d':f'{sp} 3-2@Xb4d','3-2 Xb4f':f'{sp} 3-2@Xb4f','4-3 Xb47':f'{sp} 4-3','2-1 Xb55':f'{sp} 2-1'}
    D={k:load(v,crop=14 if 'Xb4b' in v else None) for k,v in keys.items()}
    for bn,tb in BEAMS.items():
        for k,d in D.items():
            r=smooth_to(d,tb)
            if r is None: continue
            img,_=r; sg=noise(d,tb,img); nu=FR[sp][1] if '4-3' in k else (FR[sp][0] if '3-2' in k else 177.26 if sp=='HCN' else 178.375)
            fac=1.222e6/(nu**2*tb[0]*tb[1])  # K per Jy/beam, elliptical beam
            for reg,(x,y) in POS.items():
                v=sample(img,d,x,y)
                # local peak within 0.8"
                m=np.hypot(d['dx']-x,d['dy']-y)<0.8; pk=np.nanmax(np.where(m,img,np.nan))
                rows.append(dict(line=sp,transition=k,beam=bn,region=reg,m0_pos=v,m0_peak=pk,sigma=sg,Tb_pos_Kkms=v*fac,Tb_peak_Kkms=pk*fac,eTb=sg*fac,fscale_pct=10))
R=pd.DataFrame(rows); R.to_csv('/workspace/out/q8/matched_m0.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500); print(R.round(3).to_string())
