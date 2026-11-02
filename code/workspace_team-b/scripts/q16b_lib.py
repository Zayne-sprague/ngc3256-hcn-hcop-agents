# shared helpers for Q16 round 2 (template deblend, injection-recovery). Maps come from out/q16b/maps.npz (s76).
import sys,numpy as np; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
from regions import REGIONS
from scipy.ndimage import map_coordinates
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import nnls
M=np.load('/workspace/out/q16b/maps.npz')
RA0,DE0=156.96343,-43.90461; CD=np.cos(np.radians(DE0))
NU={'HCN':{2:177.26111,3:265.88618,4:354.50547},'HCO+':{2:178.37506,3:267.55763,4:356.73424}}
B0=(5.72,3.54,82.0)
def off(ra,dec): return (ra-RA0)*CD*3600,(dec-DE0)*3600     # (east, north) arcsec
POS={k:off(REGIONS[k][0],REGIONS[k][1]) for k in REGIONS}
def grid(tag,mol):
    ra,dec=M[f'{tag}|{mol}|ra'],M[f'{tag}|{mol}|dec']; e,n=off(ra,dec); return e,n
def sample(tag,mol,e_t,n_t,field='m0'):
    """bilinear sample of a map on its own regular grid at offsets (e_t,n_t) in arcsec"""
    m=M[f'{tag}|{mol}|{field}']; e,n=grid(tag,mol)
    # pixel coordinates via least-squares affine fit of (e,n) vs (x,y)
    Y,X=np.mgrid[:m.shape[0],:m.shape[1]]; A=np.c_[X.ravel(),Y.ravel(),np.ones(X.size)]
    ce=np.linalg.lstsq(A,e.ravel(),rcond=None)[0]; cn=np.linalg.lstsq(A,n.ravel(),rcond=None)[0]
    Mx=np.array([[ce[0],ce[1]],[cn[0],cn[1]]]); b=np.array([ce[2],cn[2]])
    xy=np.linalg.solve(Mx,np.vstack([np.ravel(e_t)-b[0],np.ravel(n_t)-b[1]]))
    return map_coordinates(np.nan_to_num(m),[xy[1],xy[0]],order=1,mode='nearest').reshape(np.shape(e_t))
def tb(flux,mol,J): return flux/NU[mol][J]**2     # RJ brightness T (K km/s up to a constant, same beam): ratios only
def ratio43(a4,a3,mol): return (a4/a3)*(NU[mol][3]/NU[mol][4])**2
def conv_map(img,ps,b0,tgt):
    return S.conv(img[None].astype(np.float32),ps,b0,tgt)[0]
