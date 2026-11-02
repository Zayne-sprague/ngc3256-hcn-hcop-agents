# Q5a: cache native-resolution (0.25") Xb4b HCN/HCO+ 3-2 cubes, +-6" about phase centre, per-pixel linear baseline from the Q3 line-free channels.
# python3 -I scripts/q5a_cache.py -> out/q5/cube_native_Xb4b.npz (also Xb4d native for the astrometry check, +-9")
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
for ous,crop in (('Xb4b',6.),('Xb4d',9.)):
    d={}
    for sp in ('HCN','HCO+'):
        c=vcube(sp,ous,None,crop); r=pixbase(c,FREE[sp]); v=c['v']
        fm=np.zeros(len(v),bool)
        for a,b in FREE[sp]: fm|=(v>=a)&(v<=b)
        rms=1.4826*np.median(np.abs(r[fm]-np.median(r[fm],0)),0)   # per-pixel per-channel rms (pbcor units)
        k='H' if sp=='HCN' else 'O'
        d[k]=r.astype(np.float32); d['v'+k]=v; d['rms'+k]=rms.astype(np.float32); d['beam'+k]=np.array(c['beam'])
        d['dx']=c['dx']; d['dy']=c['dy']; d['pix']=c['pix']; d['dv'+k]=c['dv']
        print(ous,sp,r.shape,c['pix'],c['beam'],c['dv'],'median rms mJy',np.median(rms)*1e3)
    np.savez(f'/workspace/out/q5/cube_native_{ous}.npz',**d)
