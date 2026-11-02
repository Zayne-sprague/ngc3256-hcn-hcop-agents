# Q12 step b: convolve the native moment-0 maps (s56a) to common circular beams (and to the Xb4f header beam), read the flux density at the six positions
# (N,S,Bridge,E1,E2,W + Sblue) and curves of growth of circular aperture sums around the N-S midpoint, for Xb4b/Xb4d/Xb4f, HCO+ and HCN 3-2.
# Conventions: Jy/beam km/s at the target beam (rescaled by beam-area ratio, as s26/s30); aperture sums are divided by the target beam area (Jy km/s);
# errors = rms over 16 sign-flipped line-free pseudo-moment realisations (same convolution), per position/aperture.
# usage: python3 -I s56b_measure.py -> out/q12/measure.csv (long format), out/q12/growth.csv
import sys,csv; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from s30_matchedJ_ratios import transfer
from scipy.ndimage import map_coordinates
CEN=[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]
SPW={'HCO+':{'Xb4f':'spw16','Xb4d':'spw25','Xb4b':'spw25'},'HCN':{'Xb4f':'spw18','Xb4d':'spw27','Xb4b':'spw27'}}
POS=['N','S','Bridge','E1','E2','W','Sblue']
BEAMS=[('1.05',(1.05,1.05,0)),('2.0',(2.0,2.0,0)),('3.0',(3.0,3.0,0)),('4.5',(4.5,4.5,0)),('6.5',(6.5,6.5,0)),('10',(10.,10.,0)),('xb4f',(5.73,3.54,82.0))]
RADII=[1,2,3,4,6,8,10,12,15]
def conv(m,ps,b0,tgt):
    n=m.shape[-1]; P=2*n; T=transfer((P,P),ps,b0,tgt)[:,:P//2+1]; fac=tgt[0]*tgt[1]/(b0[0]*b0[1]); out=[]
    for a2 in np.atleast_3d(m.reshape((-1,n,n)) if m.ndim==3 else m[None]):
        a=np.zeros((P,P),np.float32); a[:n,:n]=a2; out.append(np.fft.irfft2(np.fft.rfft2(a)*T,s=(P,P))[:n,:n]*fac)
    return np.array(out)
rows=[];grow=[]
for line in ('HCO+','HCN'):
    f0={'HCO+':267.55763e9,'HCN':265.88618e9}[line]
    for ous in ('Xb4b','Xb4d','Xb4f'):
        d=np.load(f'/workspace/out/q12/maps_{ous}_{line}.npz'); m0=d['m0']; nz=d['noise']; ps=float(d['ps']); k=int(d['k']); x0=int(d['x0']); y0=int(d['y0']); b0=tuple(d['beam'])
        hd=fits.getheader(files(f'*{ous}.ngc_3256_sci.{SPW[line][ous]}.cube.I.pbcor.fits')[0]); wc=WCS(hd).celestial
        def sky2idx(ra,dec):
            px,py=wc.world_to_pixel_values(ra,dec); return (px-x0-(k-1)/2)/k,(py-y0-(k-1)/2)/k
        n=m0.shape[0]; Y,X=np.mgrid[:n,:n]; cx,cy=sky2idx(*CEN); rr=np.hypot(X-cx,Y-cy)*ps
        for bn,tgt in BEAMS:
            if bn in ('1.05','2.0','3.0','4.5') and ous=='Xb4f': continue
            if bn!='xb4f' and (b0[0]>tgt[0] or b0[1]>tgt[0]): continue
            M=conv(np.concatenate([m0[None],nz]),ps,b0,tgt); mm,nn=M[0],M[1:]
            area=1.1331*tgt[0]*tgt[1]/ps**2
            for p in POS:
                x,y=sky2idx(REGIONS[p][0],REGIONS[p][1])
                val=float(map_coordinates(mm,[[y],[x]],order=1)[0]); er=float(np.std([map_coordinates(q,[[y],[x]],order=1)[0] for q in nn]))
                rows.append(dict(line=line,ous=ous,beam=bn,pos=p,val=val,err=er))
            for R in RADII:
                mk=rr<=R; s=float(mm[mk].sum()/area); e=float(np.std([q[mk].sum()/area for q in nn]))
                grow.append(dict(line=line,ous=ous,beam=bn,R=R,val=s,err=e))
            # annulus 8-20 (flux in the outskirts, per beam area)
        print(line,ous,'done',flush=True)
for fn,r in (('measure',rows),('growth',grow)):
    with open(f'/workspace/out/q12/{fn}.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(r[0])); w.writeheader(); w.writerows(r)
