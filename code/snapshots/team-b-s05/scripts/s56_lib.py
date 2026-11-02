# Q12 shared helpers: convolution of native moment-0 maps (s56a) to circular beams and resampling onto the Xb4f 0.7" grid (+-19" crop) for Fourier-domain work.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from s30_matchedJ_ratios import transfer
from scipy.ndimage import map_coordinates
CEN=[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]
SPW={'HCO+':{'Xb4f':'spw16','Xb4d':'spw25','Xb4b':'spw25'},'HCN':{'Xb4f':'spw18','Xb4d':'spw27','Xb4b':'spw27'}}
LAMBDA={'Xb4f':1.1296e-3}   # m, 265.4 GHz
LAM=2.99792458e8/265.7e9    # observing wavelength of HCO+/HCN 3-2 [m]
ARC=206265.0
def q2L(q): return q*ARC*LAM            # spatial frequency [cycles/arcsec] -> baseline [m]
def L2q(L): return L/(ARC*LAM)
def conv(m,ps,b0,tgt):
    """m: (n,n) or (k,n,n) native-beam Jy/beam maps -> target beam, flux rescaled by beam-area ratio (as s26/s30)"""
    m=np.asarray(m,np.float32); m=m[None] if m.ndim==2 else m
    n=m.shape[-1]; P=2*n; T=transfer((P,P),ps,b0,tgt)[:,:P//2+1]; fac=tgt[0]*tgt[1]/(b0[0]*b0[1]); out=[]
    for a2 in m:
        a=np.zeros((P,P),np.float32); a[:n,:n]=a2; out.append(np.fft.irfft2(np.fft.rfft2(a)*T,s=(P,P))[:n,:n]*fac)
    return np.array(out)
class Arr:
    def __init__(s,ous,line):
        d=np.load(f'/workspace/out/q12/maps_{ous}_{line}.npz'); s.m0=d['m0']; s.noise=d['noise']; s.ps=float(d['ps']); s.k=int(d['k']); s.x0=int(d['x0']); s.y0=int(d['y0']); s.b0=tuple(d['beam'])
        hd=fits.getheader(files(f'*{ous}.ngc_3256_sci.{SPW[line][ous]}.cube.I.pbcor.fits')[0]); s.wc=WCS(hd).celestial; s.ous=ous
    def sky2idx(s,ra,dec):
        px,py=s.wc.world_to_pixel_values(ra,dec); return (px-s.x0-(s.k-1)/2)/s.k,(py-s.y0-(s.k-1)/2)/s.k
    def at_beam(s,tgt):
        M=conv(np.concatenate([s.m0[None],s.noise]),s.ps,s.b0,tgt); return M[0],M[1:]
_ref={}
def common_grid(line):
    """sky coordinates of the Xb4f crop pixels (the common 0.7" grid)"""
    if line in _ref: return _ref[line]
    f=Arr('Xb4f',line); n=f.m0.shape[0]; Y,X=np.mgrid[:n,:n]
    ra,dec=f.wc.pixel_to_world_values(f.x0+X*f.k,f.y0+Y*f.k); _ref[line]=(ra,dec,f.ps); return _ref[line]
def on_grid(a,M,line):
    ra,dec,_=common_grid(line); x,y=a.sky2idx(ra,dec); return map_coordinates(M,[y,x],order=1,mode='nearest')
def maps_on_grid(ous,line,tgt,A=None):
    a=A or Arr(ous,line); M,N=a.at_beam(tgt); return on_grid(a,M,line),np.array([on_grid(a,q,line) for q in N]),a
