# Q23: measured low-edge profile (stack of N,S,Bridge spectrum/continuum, 6-ch bins) for the three ramp spws used by the ladder lines, vs shapes.
# Tests whether the ramp has died away before the line window (red wing starts at v=+350): window start ch = chl - 350/dv.
import sys, json, numpy as np; sys.path.insert(0,'/workspace/scripts')
from common import vrad, VSYS
S131=json.load(open('/workspace/out/q22/s131.json'))
F0={'HCN':{3:265.88618e9,4:354.50547e9},'HCO+':{2:178.37506e9}}
for ous,spw,mol,J in (('Xb55',18,'HCO+',2),('Xb4f',18,'HCN',3),('Xb47',16,'HCN',4)):
    z=np.load(f'/workspace/out/regspec/{ous}.spw{spw}.cube.pbcor.npz'); v=vrad(z['nu'],F0[mol][J])-VSYS; dv=abs(np.median(np.diff(v))); chl=int(np.argmin(abs(v))); w0=chl-int(350/dv)
    R=['N','S','Bridge']; cont={r:S131[f'{ous}.spw{spw}.{r}']['cont'] for r in R}
    rms=np.mean([np.std(np.nan_to_num(z[r+'_pix'])[v<-450]*1e3) for r in R]); c=np.mean([cont[r] for r in R])
    bw=6; nb=w0//bw+2
    st=[];er=[]
    for k in range(nb):
        sl=slice(k*bw,(k+1)*bw)
        a=np.sum([np.nan_to_num(z[r+'_pix'])[sl].mean()*1e3*cont[r] for r in R])/np.sum([cont[r]**2 for r in R])  # least-squares amp per cont
        e=rms*1.3/np.sqrt(bw)/np.sqrt(np.sum([cont[r]**2 for r in R])); st.append(a); er.append(e)
    print(f'\n{ous} spw{spw} ({mol}{J}-{J-1}) line ch {chl}, +350 km/s window starts ch {w0}; ramp amplitude (x cont) per {bw}-ch bin from low edge:')
    print('  bins ch  :',' '.join(f'{k*bw:5d}' for k in range(nb)))
    print('  f +-     :',' '.join(f'{a:+5.2f}' for a in st),'  err ~',round(np.mean(er),3))
    lo=np.average(st[:3],weights=1/np.square(er[:3])); mid=np.average(st[3:6],weights=1/np.square(er[3:6])); 
    wn=[k for k in range(nb) if w0-12<=k*bw<w0+bw]; hi=np.average([st[k] for k in wn],weights=[1/er[k]**2 for k in wn])
    print(f'  mean first 18 ch {lo:+.3f}; ch 18-36 {mid:+.3f}; last 12 ch before/at window start (ch {w0-12}-{w0+bw}) {hi:+.3f} +- {np.mean([er[k] for k in wn])/np.sqrt(len(wn)):.3f}')
