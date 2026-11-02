# reviewer check: is ridge wing flux (and lobe) beyond what core brightness produces? wing vs core over the whole region, red and blue windows
import sys; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']
    c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
    def Wm(lo,hi):
        k=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
    h=(np.hypot(*np.mgrid[-8:9,-8:9])*ps<=1.0).astype(float)/R['ba']
    A={k:fftconvolve(Wm(*w),h,mode='same') for k,w in dict(core=(-150,150),red=(180,320),blue=(-320,-180),red2=(-320+0,-180)).items()}
    # control: wing at |v| 400-540 (beyond ridge), measures residual leakage scale
    A['redfar']=fftconvolve(Wm(400,540),h,mode='same'); A['bluefar']=fftconvolve(Wm(-540,-400),h,mode='same')
    fm=field_pixels(R,rmin=3.5,rmax=22.0)
    print(mol,'field sigma red/blue',np.std(A['red'][fm]),np.std(A['blue'][fm]),'far',np.std(A['redfar'][fm]),np.std(A['bluefar'][fm]))
    core=A['core']; 
    for lo,hi in ((0.5,1),(1,2),(2,4),(4,8),(8,100)):
        m=(core>lo)&(core<=hi)&(np.hypot(*np.mgrid[:n,:n].__sub__(n//2) if False else np.mgrid[:n,:n])>=0)
        # restrict to central 8" box
        yy,xx=np.mgrid[:n,:n]; cy,cx=R['pxy'](*S0)[::-1] if False else (None,None)
        print(f' core bin {lo}-{hi}: npix {m.sum()}; red mean {A["red"][m].mean():.3f} sd {A["red"][m].std():.3f}; blue mean {A["blue"][m].mean():.3f} sd {A["blue"][m].std():.3f}; redfar {A["redfar"][m].mean():.3f}; bluefar {A["bluefar"][m].mean():.3f}')
    for nm,pos in (('ridge',offset_pos(S0,.15,3.3)),('lobe',offset_pos(S0,.15,-1.85)),('N',N0),('S',S0)):
        m=apmask(n,ps,R['pxy'](*pos),1.0)
        print(' ',nm,'core(mean ap)',A['core'][m].mean() if False else float(Wm(-150,150)[m].sum()/R['ba']),'red',float(Wm(180,320)[m].sum()/R['ba']),'blue',float(Wm(-320,-180)[m].sum()/R['ba']),'redfar',float(Wm(400,540)[m].sum()/R['ba']),'bluefar',float(Wm(-540,-400)[m].sum()/R['ba']))
