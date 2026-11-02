# Q28 round 2 library: joint N+S two-phase model (compact C per nucleus + ONE shared extended E), fine grids, beam overlap, two-sided T_b.
# Conventions: data vector (14) per nucleus, all in "W*T_R" units (K km/s): 10" [HCN 2-1,3-2,4-3, HCO+ 2-1,3-2,4-3], 2.2" [HCN 1-0, HCN 3-2, HCO+ 1-0, HCO+ 3-2],
#  1.29" [HCN 3-2, HCO+ 3-2] (T_b x159.7), 0.24" [HCN 3-2, HCO+ 3-2] peak T_b x159.7 (two-sided data now, not floors).
# Model C (compact): Gaussian FWHM theta, brightness filling inside = fC (clumping), coupling fC*theta^2/(theta^2+b^2); at the OTHER nucleus (d=5.2") the 10" beam sees it with
#   extra factor exp(-4ln2 d^2/(theta^2+b^2)). Model E: smooth, uniform over both nuclei (shared), filling fE. Shared line width dv=150 km/s.
import sys,json,numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L
K=1.00927; W=159.7; SEP=5.2; PCAS=205.   # pc per arcsec (49 pc / 0.24")
BEAMS=np.array([10.]*6+[2.2]*4+[1.29]*2+[0.24]*2); SEL=[1,2,3,5,6,7, 0,2,4,6, 2,6, 2,6]
d=np.load('/workspace/out/q3d/grid_lvg500.npz'); LTK=np.log10(d['TK']); LN=d['LN']; LNC=d['LNC']
def mk(a): return RGI((LTK,LN,LNC),a,bounds_error=False,fill_value=None)
_tr={m:np.log(np.clip(np.nan_to_num(d[m][...,0],nan=0.),1e-6,1e3)) for m in ('hcn','hcop')}
_ta={m:np.log(np.clip(np.nan_to_num(d[m][...,2],nan=0.),1e-8,1e4)) for m in ('hcn','hcop')}
ITR={m:[mk(_tr[m][...,j]) for j in range(4)] for m in _tr}; ITA={m:[mk(_ta[m][...,j]) for j in range(4)] for m in _ta}
def models(lt,ln,lN,A):
    """grid of (logTk, logn, logN(HCO+), logA) -> params [M,4], T_R [M,8] (HCN 1-4, HCO+ 1-4), tau [M,8]"""
    g=np.array(np.meshgrid(lt,ln,lN,A,indexing='ij')).reshape(4,-1).T
    ok=(g[:,2]+g[:,3]<=LNC[-1]+1e-6)&(g[:,2]>=LNC[0])
    g=g[ok]; ph=np.c_[g[:,0],g[:,1],g[:,2]+g[:,3]]; po=g[:,:3]
    tr=np.stack([np.exp(ITR['hcn'][j](ph)) for j in range(4)]+[np.exp(ITR['hcop'][j](po)) for j in range(4)],1)
    ta=np.stack([np.exp(ITA['hcn'][j](ph)) for j in range(4)]+[np.exp(ITA['hcop'][j](po)) for j in range(4)],1)
    return g,np.minimum(np.nan_to_num(tr),1e3),np.nan_to_num(ta,posinf=1e6)
r183=json.load(open('/workspace/out/q29/s183.json')); r192=json.load(open('/workspace/out/q34/s192.json'))['res']; s160=json.load(open('/workspace/out/q6r2/s160.json'))
def data(pos,sig=(0.10,0.15,0.15,0.15)):
    I10,e10=L.DATA[pos]; a=r183['2015.1.00993 X961|2.2|'+pos]; c=r192['2p2|comb|'+pos]
    I22=np.array([a['HCN']['h150_ap'][0],c['HCN32|ap|h150'][0],a['HCOp']['h150_ap'][0],c['HCOp32|ap|h150'][0]])*K
    I13=np.array([s160['HCN3-2'][pos]['250pc(1.29")']['Tb'],s160['HCO+3-2'][pos]['250pc(1.29")']['Tb']])*W
    I02=np.array([s160['HCN3-2'][pos]['native 0.24"']['Tb'],s160['HCO+3-2'][pos]['native 0.24"']['Tb']])*W
    y=np.concatenate([I10,I22,I13,I02])
    sg=np.concatenate([y[:6]*np.sqrt(np.asarray(e10)**2+sig[0]**2),I22*sig[1],I13*sig[2],I02*sig[3]])
    return y,sg
cp=lambda th,b:th**2/(th**2+b**2)
THS=np.array([0.15,0.2,0.3,0.5,0.8,1.2]); FC=np.array([0.1,0.3,1.0])
def cscale(th,fc,off=0.):
    """[nth,nfc,14] compact coupling; off>0 -> at the other nucleus (only 10" matters)"""
    s=np.empty((len(th),len(fc),14))
    for i,t in enumerate(th):
        c=cp(t,BEAMS)*np.where(BEAMS==10.,np.exp(-4*np.log(2)*off**2/(t**2+BEAMS**2)),1.) if off>0 else cp(t,BEAMS)
        if off>0: c=np.where(BEAMS==10.,c,0.)
        s[i]=fc[:,None]*c[None,:]
    return s
