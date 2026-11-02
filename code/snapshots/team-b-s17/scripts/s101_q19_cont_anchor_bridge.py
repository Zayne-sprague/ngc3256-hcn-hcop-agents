# Q19 round 2: continuum anchor at the LINE beam (Xb4f native 5.72x3.54 PA82; Xb47 and Xb4f continuum both convolved to it), N, S, Bridge, E1;
# aperture spread -> anchor error; Bridge control; interval-of-scenarios reporting.  usage: python3 -I s101_q19_cont_anchor_bridge.py -> out/q19/s101.txt/json
import sys,json; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
from scipy.signal import fftconvolve
rng=np.random.default_rng(101); OUT='/workspace/out/q19/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
RA0,DE0=156.96343,-43.90461
def cov(a,b,pa):
    th=np.deg2rad(pa); e1=np.array([np.sin(th),np.cos(th)]); e2=np.array([np.cos(th),-np.sin(th)]); return (a/2.3548)**2*np.outer(e1,e1)+(b/2.3548)**2*np.outer(e2,e2)
def convolve_to(d,ps,b0,tgt):
    C2=cov(*tgt)-cov(*b0); w_,U=np.linalg.eigh(C2)
    assert w_.min()>-1e-3*w_.max()
    C2=U@np.diag(np.clip(w_,1e-6,None))@U.T
    n=int(8*np.sqrt(w_.max())/ps)|1
    ky,kx=np.mgrid[:n,:n]-n//2; q=np.stack([(-kx*ps).ravel(),(ky*ps).ravel()]); Ci=np.linalg.inv(C2)
    ker=np.exp(-0.5*np.einsum('ij,ij->j',q,Ci@q)).reshape(n,n); ker/=ker.sum()
    return fftconvolve(np.nan_to_num(d),ker,mode='same')*(tgt[0]*tgt[1])/(b0[0]*b0[1])
TGT=(5.72,3.54,82.0)
CONT={'Xb47':'*Xb47*16_18_20_22.cont.I.pbcor.fits','Xb4f':'*Xb4f*16_18_20_22.cont.I.pbcor.fits'}
PO=['N','S','Bridge','E1']
val={}
for o,pat in CONT.items():
    d,h=load2d(files(pat)[0]); w=WCS(h).celestial; ps=abs(h['CDELT1'])*3600; b0=(h['BMAJ']*3600,h['BMIN']*3600,h['BPA'])
    c=convolve_to(d,ps,b0,TGT) if abs(b0[0]-TGT[0])>0.03*TGT[0] or abs(b0[1]-TGT[1])>0.03*TGT[1] else d.copy()
    P(f'# {o}: native beam {b0[0]:.2f}x{b0[1]:.2f} PA{b0[2]:.0f} -> {TGT}')
    Y,X=np.mgrid[:d.shape[0],:d.shape[1]]
    for p in PO:
        ra,dec=REGIONS[p][0],REGIONS[p][1]; x,y=[float(np.squeeze(t)) for t in w.world_to_pixel_values(ra,dec)]
        r=np.hypot(X-x,Y-y)*ps
        v={'peak':float(np.max(np.where(r<=1.0,c,-1e9))),'pt':float(c[int(round(y)),int(round(x))])}
        for R in (1.5,2.5):   # mean of the convolved map within R (beam-matched smooth aperture)
            v[f'm{R}']=float(np.mean(c[r<=R]))
        val[(o,p)]=v
ap=['peak','pt','m1.5','m2.5']
P('# continuum flux (mJy/beam, Xb4f beam)'); 
for p in PO: P(f'  {p:7s} '+'  '.join(f'{o}:{ap_}={val[(o,p)][ap_]*1e3:.2f}' for o in CONT for ap_ in ap))
Dc={}
def D(a,b,k): return (val[('Xb47',a)][k]/val[('Xb4f',a)][k])/(val[('Xb47',b)][k]/val[('Xb4f',b)][k])
for a,b in (('N','S'),('N','Bridge'),('S','Bridge'),('E1','S')):
    Dc[f'{a}/{b}']={k:D(a,b,k) for k in ap}
    P(f'# continuum (47/4f) ratio-of-ratios {a}/{b}: '+', '.join(f'{k} {v:.3f}' for k,v in Dc[f'{a}/{b}'].items()))
# observed line contrasts (s82, Xb4f beam; stat+centring) HCN:
obs={'N/S':(1.262,0.074),'N/Bridge':(1.25,0.13),'S/Bridge':(0.99,0.11)}
obsO={'N/S':(1.164,0.064),'S/Bridge':(0.98,0.07)}
sig_c_aper=lambda key: max(0.10,np.std(list(Dc[key].values()),ddof=1))
n=200000; res={}
def summ(name,x):
    q=np.percentile(x,[2.5,16,50,84,97.5]); res[name]=dict(p=list(q),P_gt1=float((x>1).mean())); P(f'  {name:58s} {q[2]:.3f} 68% [{q[1]:.3f},{q[3]:.3f}] 95% [{q[0]:.3f},{q[4]:.3f}] P(>1)={(x>1).mean():.2f}')
for key,(m,s) in obs.items():
    L=rng.normal(m,s,n)
    if key=='N/S': L=L*rng.normal(0.98,0.025,n)   # Q23 ramp
    sc=sig_c_aper(key); P(f'\n## HCN {key}: observed {m}+-{s}; continuum anchor mean over apertures {np.mean(list(Dc[key].values())):.3f}, spread sd {np.std(list(Dc[key].values()),ddof=1):.3f}, sig_c used {sc:.3f}')
    P(f'  D needed for unity: {m*(0.98 if key=="N/S" else 1):.2f}')
    summ(f'{key} no recovery',L)
    for ak in ap:
        for da in (0.0,0.3):
            Dv=Dc[key][ak]*np.exp(rng.normal(0,sc,n))/(346.5/257.7)**rng.normal(0,da,n); summ(f'{key} anchor Xb4f beam {ak} dalpha sd {da}',L/Dv)
    Dall=np.exp(rng.normal(np.log(np.mean(list(Dc[key].values()))),sc,n)); summ(f'{key} anchor aperture-mean, sig_c={sc:.2f}',L/Dall)
json.dump(dict(res=res,Dc=Dc),open(OUT+'s101.json','w'),indent=1,default=float); open(OUT+'s101.txt','w').write('\n'.join(log))
# pool sensitivity (HCN N/S): A (Q12 priors) and B (line-beam anchor, aperture mean, sig_c 0.1, dalpha sd 0.3); C (team-a) excluded; weights wA
L=rng.normal(1.262,0.074,n)*rng.normal(0.98,0.025,n); A_=L/rng.uniform(1,0.92/0.89,n)
B_=L/(np.exp(rng.normal(np.log(np.mean(list(Dc['N/S'].values()))),0.10,n))/(346.5/257.7)**rng.normal(0,0.3,n))
P('\n## pool A+B (C excluded) HCN N/S, weight on A:')
for wA in (0.25,0.5,0.75): summ(f'pool A+B wA={wA}',np.concatenate([A_[:int(wA*n)],B_[:int((1-wA)*n)]]))
json.dump(dict(res=res,Dc=Dc),open(OUT+'s101.json','w'),indent=1,default=float); open(OUT+'s101.txt','w').write('\n'.join(log))
