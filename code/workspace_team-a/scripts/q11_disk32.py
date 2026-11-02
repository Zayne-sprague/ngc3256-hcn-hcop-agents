# Q11: like-for-like disk annuli HCN/HCO+ at 3-2 (Xb4f, Xb4d) vs 2-1 (Xb55), same beam, annuli about N-S midpoint, zero-level variants, null-window noise.
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import json
TB={'10':(10.,10.,0.),'Xb55':(9.36,6.28,-88.)}
KEYS={'2-1@Xb55':('HCN 2-1','HCO+ 2-1'),'3-2@Xb4f':('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'),'3-2@Xb4d':('HCN 3-2@Xb4d','HCO+ 3-2@Xb4d')}
ANN=[(0,4),(4,8),(8,12)]
out={}
def smooth_tb(d,tb,arr=None):
    k=kernel(d['hd'],tb,d['beam'])
    if k is None: return None
    a=np.nan_to_num(d['m0'] if arr is None else arr); return fftconvolve(a,k,mode='same')*OM(tb)/OM(d['beam'])
def xy_mid(d): return -0.075,0.0     # midpoint of N (0,+2.55) and S (-0.15,-2.55)
for tbn,tb in TB.items():
  for tag,(kh,ko) in KEYS.items():
    D=[load(k,nullspec=True) for k in (kh,ko)]
    S=[smooth_tb(d,tb) for d in D]
    if S[0] is None or S[1] is None: print(tbn,tag,'native beam larger'); continue
    res={}
    for d,s,k in zip(D,S,(kh,ko)):
        r=np.hypot(d['dx']-xy_mid(d)[0],d['dy']-xy_mid(d)[1]); d['r']=r
        d['s']=s; d['ns']=[smooth_tb(d,tb,n) for _,n in d['nulls']]
    rmax=D[0]['r'].max(); 
    # zero-level options: none; far-field annulus (as available); median of null-window maps in annuli (should be ~0 baseline)
    far=[(30,45)] if rmax>40 else [(20,32)]
    for fl in [None]+far:
      for a,b in ANN:
        row=[]
        for d in D:
            r=d['r']; m=(r>=a)&(r<b)
            z=np.nanmean(d['s'][(r>=fl[0])&(r<fl[1])]) if fl else 0.
            omp=OM(tb)/d['pix']**2
            f=(d['s'][m]-z).sum()/omp
            # noise: scatter of smoothed null maps (annulus-summed), empirical
            nv=[(n[m]).sum()/omp for n in d['ns']]
            row.append((f,np.std(nv) if len(nv)>2 else np.nan,m.sum()*d['pix']**2/OM(tb)))
        (fH,eH,nb),(fO,eO,_)=row
        R=fH/fO; eR=abs(R)*np.hypot(eH/fH if fH else 0,eO/fO)
        print(f'{tbn:5s} {tag:9s} zero={fl} {a:2d}-{b:2d}": HCN {fH:7.2f}+-{eH:.2f} HCO+ {fO:7.2f}+-{eO:.2f} Jykm/s  R={R:.3f}+-{eR:.3f} (nbeams {nb:.1f}, rmax {rmax:.0f}")')
        out[f'{tbn}|{tag}|{fl}|{a}-{b}']=dict(HCN=fH,eHCN=eH,HCO=fO,eHCO=eO,R=R,eR=eR)
json.dump(out,open('/workspace/out/q11/q11_disk32.json','w'),indent=1,default=float)
