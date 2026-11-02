# Q5c: are the 0.25" sub-structures real? HCN-HCO+ cross-correlation of high-pass maps and maximum coincidence vs shifted nulls; census of matched clumps.
# python3 -I scripts/q5c_clumps.py -> out/q5/q5c.log
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
from scipy import ndimage as ndi
d=load(); dx=d['dx']; dy=d['dy']; pix=float(d['pix']); log=open(O+'q5c.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
mH,_=win_map(d,'H',*VW); mO,_=win_map(d,'O',*VW)
REG={'N':(0.,2.6,2.2),'S':(-0.15,-2.55,2.0),'Bridge':(0.3,0.15,1.5),'all5':(0,0,5.5)}
def hp(m,s_pix=21):   # remove >1" scale (sigma 1"=21 px)
    return m-ndi.gaussian_filter(m,s_pix)
hH,hO=hp(mH),hp(mO)
# pre-smooth by 0.1" (2 px) sigma (a bit below the beam) so noise speckle is suppressed before correlating
sH,sO=ndi.gaussian_filter(hH,2),ndi.gaussian_filter(hO,2)
def corr(a,b,mask):
    a=a[mask]; b=b[mask]; a=a-a.mean(); b=b-b.mean(); return (a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum())
rng=np.random.default_rng(5)
for nm,(x,y,r) in REG.items():
    msk=np.hypot(dx-x,dy-y)<r
    c0=corr(sH,sO,msk)
    # null: shift HCO+ by random offsets 0.8..2" (keep mask fixed), circular shift edge effects negligible with np.roll on 255px
    nl=[]
    for _ in range(300):
        ang=rng.uniform(0,2*np.pi); rr=rng.uniform(0.8,2.5)/pix; sh=(int(rr*np.sin(ang)),int(rr*np.cos(ang)))
        nl.append(corr(sH,np.roll(sO,sh,axis=(0,1)),msk))
    nl=np.array(nl); pr(f'{nm}: corr(HCN,HCO+ high-pass) lag0 = {c0:.2f}; shifted-null mean {nl.mean():.2f} sd {nl.std():.2f} -> {(c0-nl.mean())/nl.std():.1f} sigma')
# matched maxima census (maxima of 0.1"-smoothed m0 above thr, spacing >= 0.3")
def maxima(m,thr,r=5.5,sep=6):
    s=ndi.gaussian_filter(m,2.0); mx=(s==ndi.maximum_filter(s,size=2*sep+1))&(np.hypot(dx,dy)<r)&(s>thr)
    iy,ix=np.nonzero(mx); return np.array([(dx[a,b],dy[a,b],s[a,b]) for a,b in zip(iy,ix)])
def matched(A,B,tol=0.2):
    return [(a,b) for a in A for b in B if np.hypot(a[0]-b[0],a[1]-b[1])<tol]
pr('\n# matched local maxima (within 0.2") of 0.1"-smoothed m0 (+-350 km/s), thresholds in Jy/b km/s; null = HCO+ map shifted by 1-2.5" in random directions')
for thr in (0.3,0.4,0.5):
    A=maxima(mH,thr); B=maxima(mO,thr); nobs=len(matched(A,B))
    nn=[]
    for _ in range(200):
        ang=rng.uniform(0,2*np.pi); rr=rng.uniform(1.0,2.5)/pix; sh=(int(rr*np.sin(ang)),int(rr*np.cos(ang)))
        Bn=maxima(np.roll(mO,sh,axis=(0,1)),thr); nn.append(len(matched(A,Bn)))
    pr(f'thr {thr}: HCN {len(A)} maxima, HCO+ {len(B)}, matched {nobs}; chance-matches {np.mean(nn):.1f}+-{np.std(nn):.1f} => excess {nobs-np.mean(nn):.1f}')
# list of matched clumps in N, S, Bridge at thr 0.4
A=maxima(mH,0.4); B=maxima(mO,0.4); M=matched(A,B)
pr('\n# matched clumps thr 0.4 (HCN x,y,m0 | HCO+ x,y,m0)')
for a,b in sorted(M,key=lambda t:-t[0][2]):
    reg=[n for n,(x,y,r) in REG.items() if n!='all5' and np.hypot(a[0]-x,a[1]-y)<r]
    pr(f'({a[0]:+.2f},{a[1]:+.2f}) {a[2]:.2f} | ({b[0]:+.2f},{b[1]:+.2f}) {b[2]:.2f}  {reg}')
# --- census on the high-pass maps (clump contrast over the >1" background); spacing >= 0.3", match tol 0.2"
pr('\n# matched maxima of high-pass (>1" removed) 0.1"-smoothed maps; HCO+ shifted null')
for thr in (0.10,0.15,0.20):
    A=maxima(hH,thr); B=maxima(hO,thr); M=matched(A,B)
    nn=[]
    for _ in range(200):
        ang=rng.uniform(0,2*np.pi); rr=rng.uniform(1.0,2.5)/pix; sh=(int(rr*np.sin(ang)),int(rr*np.cos(ang)))
        nn.append(len(matched(A,maxima(np.roll(hO,sh,axis=(0,1)),thr))))
    byreg={n:sum(1 for a,b in M if np.hypot(a[0]-x,a[1]-y)<r) for n,(x,y,r) in REG.items() if n!='all5'}
    pr(f'hp thr {thr}: HCN {len(A)}, HCO+ {len(B)}, matched {len(M)} (N {byreg["N"]}, S {byreg["S"]}, Bridge {byreg["Bridge"]}); chance {np.mean(nn):.1f}+-{np.std(nn):.1f}')
    if thr==0.15:
        for a,b in sorted(M,key=lambda t:-t[0][2]):
            pr(f'   ({a[0]:+.2f},{a[1]:+.2f}) hp {a[2]:.2f} | ({b[0]:+.2f},{b[1]:+.2f}) {b[2]:.2f}')
