# Q19 r2: direct bound on the 3-2 recovery differential: Xb4d J3 (0.94x0.74") convolved to the Xb4f beam vs Xb4f J3, at N, S, Bridge, HCN and HCO+.
# R3(x)=S(4f native)/S(4d smoothed); N/S ratio of R3 -> part of D from the 3-2 side.  usage: python3 -I s102_q19_xb4d_vs_xb4f_32.py -> out/q19/s102.txt
import sys; sys.path.insert(0,'/workspace/scripts')
from q16b_lib import *
import s30_matchedJ_ratios as S
from scipy.signal import fftconvolve
log=[]
def P(s=''): print(s,flush=True); log.append(s)
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def smooth(m,ps,b0,tgt):
    C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2); C2=U@np.diag(np.clip(w_,1e-6,None))@U.T
    n=int(8*np.sqrt(w_.max())/ps)|1; ky,kx=np.mgrid[:n,:n]-n//2; q=np.stack([(-kx*ps).ravel(),(ky*ps).ravel()])
    ker=np.exp(-0.5*np.einsum('ij,ij->j',q,np.linalg.inv(C2)@q)).reshape(n,n); ker/=ker.sum()
    return fftconvolve(np.nan_to_num(m),ker,mode='same')*(tgt[0]*tgt[1])/(b0[0]*b0[1])
for mol in ('HCN','HCO+'):
    b=M[f'4d3|{mol}|beam']; ps=float(M[f'4d3|{mol}|ps']); P(f'# {mol}: Xb4d beam array {b} ps {ps}')
    m4d=M[f'4d3|{mol}|m0']; b0=tuple(float(x) for x in b)
    sm=smooth(m4d,ps,b0,B0)
    # sample both at positions (use q16b_lib.sample on a temporary entry)
    out={}
    for p in ('N','S','Bridge'):
        def smp(arr,tag):
            e,n_=grid(tag,mol)
            from scipy.interpolate import RegularGridInterpolator
            ra=M[f'{tag}|{mol}|ra']; de=M[f'{tag}|{mol}|dec']
            ee=(ra-RA0)*CD*3600; nn=(de-DE0)*3600
            d2=(ee-POS[p][0])**2+(nn-POS[p][1])**2
            # bilinear via map_coordinates on pixel grid: nearest-pixel mean in 0.5" radius
            msk=d2<0.5**2; return float(arr[msk].mean())
        s4d=smp(sm,'4d3'); s4f=smp(M[f'4f3|{mol}|m0'],'4f3'); out[p]=s4f/s4d
        P(f'  {p:7s} S(4f native)={s4f:.3f}  S(4d->4f beam)={s4d:.3f} Jy/b km/s  R3=4f/4d {s4f/s4d:.3f}')
    P(f'  {mol} R3 N/S {out["N"]/out["S"]:.3f}  N/Bridge {out["N"]/out["Bridge"]:.3f}  S/Bridge {out["S"]/out["Bridge"]:.3f}  (R3<1: Xb4f loses vs Xb4d; N/S>1: S loses more 3-2 in Xb4f)')
open('/workspace/out/q19/s102.txt','w').write('\n'.join(log))
