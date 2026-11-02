# Q9b fix 2: beam / spatial-mixing test.
# (a) same aperture, same nuclear centre: Xb4b convolved to 0.30" vs 1.05" vs Xb4d@1.05" -> ratio vs velocity (40-km/s bins, NO S/N mask) and weighted linear slope;
# (b) HCN-vs-HCO+ centroid per 40-km/s bin around N (flux-weighted, r<1.2" at 1.05"; r<0.6" at 0.30") with S/N-based position errors;
# (c) injection: HCN = R(r_N) * HCO+ (Xb4b native cube, ratio declining with radius as measured), identical kinematics, smoothed to 1.05" -> predicted ratio-vs-velocity slope of the beam alone.
# python3 -I q9b_beam.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q9b_lib import *
import pandas as pd
rng=np.random.default_rng(5)
EDG=np.arange(-140,141,40)
def binspec(c,sp,lo,hi): return winflux(c,sp,lo,hi)
def ratio_profile(ous,tb,reg,r1,r2,Hs=None,Os=None,label=''):
    x0,y0=centre(ous,tb,reg) if Hs is None else centre(ous,tb,reg)
    out=[]
    cs={}
    for ln,s in (('HCN',Hs),('HCO+',Os)):
        c=cube(ous,ln); s=smcube(ous,ln,tb) if s is None else s; m=apmask(c,x0,y0,r1,r2); sp=apspec(c,s,m,tb)
        ns,sch=nullstats(c,s,tb,r1,r2,[40],rng,120) if s is not None else (None,None)
        cs[ln]=(c,sp,ns[40],sch)
    rows=[]
    for a,b in zip(EDG[:-1],EDG[1:]):
        F=[winflux(cs[l][0],cs[l][1],a,b) for l in ('HCN','HCO+')]; e=[np.hypot(cs[l][2],cs[l][3]*cs[l][0].dv*np.sqrt(basevar(cs[l][0],a,b))) for l in ('HCN','HCO+')]
        R=F[0]/F[1]; rows.append(dict(v=(a+b)/2,HCN=F[0],HCO=F[1],R=R,eR=abs(R)*np.hypot(e[0]/F[0],e[1]/F[1])))
    g=pd.DataFrame(rows); w=1/g.eR; A=np.vstack([np.ones(len(g)),g.v/100.]).T*w.values[:,None]; yv=g.R.values*w.values
    cov=np.linalg.inv(A.T@A); p=cov@A.T@yv; chi=((yv-A@p)**2).sum(); Rm=np.average(g.R,weights=w**2); chic=(((g.R-Rm)*w)**2).sum()
    return g,dict(label=label,ous=ous,beam=tb[0],reg=reg,r1=r1,r2=r2,Rmean=Rm,chi2_const=chic,dof_const=len(g)-1,slope100=p[1],eslope100=np.sqrt(cov[1,1]),R_at_0=p[0],chi2_lin=chi)
res=[]
for reg in ('N','S'):
    for ous,tb in (('Xb4d',(1.05,)*2+(0.,)),('Xb4b',(1.05,)*2+(0.,)),('Xb4b',(0.30,)*2+(0.,))):
        for r1,r2 in ((0,.5),(0,1.0)) :
            g,d=ratio_profile(ous,tb,reg,r1,r2,label=f'{ous}@{tb[0]}'); res.append(d)
            print(f'\n{reg} {ous}@{tb[0]}" r={r1}-{r2}"'); print(g.round(3).to_string(index=False)); print({k:(round(v,3) if isinstance(v,float) else v) for k,v in d.items()})
pd.DataFrame(res).to_csv('/workspace/out/q9b/beam_slopes.csv',index=False)
# (b) centroids
print('\n== (b) centroids (arcsec offsets of HCN minus HCO+ per bin, N, S) ==')
cen=[]
for ous,tb,rad in (('Xb4d',(1.05,1.05,0.),1.2),('Xb4b',(0.30,0.30,0.),0.6)):
    for reg in ('N','S'):
        x0,y0=centre(ous,tb,reg); P={}
        for ln in ('HCN','HCO+'):
            c=cube(ous,ln); s=smcube(ous,ln,tb); m=np.hypot(c.X-x0,c.Y-y0)<rad
            free=np.where(np.abs(c.v)>550)[0]; 
            def bmap(a,b): k=(c.v>=a)&(c.v<b); return s[k].sum(0)*c.dv
            nb=[]
            for lo in np.arange(-900,-540,40): 
                if lo>=c.v.min()+20: nb.append(bmap(lo,lo+40)[m])
            for lo in np.arange(560,900,40):
                if lo+40<=c.v.max()-20: nb.append(bmap(lo,lo+40)[m])
            sig=1.4826*np.median(np.abs(np.array(nb)-np.median(nb)))     # per-pixel noise of a 40-km/s map (Jy/beam km/s)
            P[ln]=(c,s,m,sig,bmap)
        for a in EDG[:-1]:
            r={}
            for ln,(c,s,m,sig,bmap) in P.items():
                im=bmap(a,a+40); w=np.where(m&(im>3*sig),im,0.); tot=w.sum()
                pk=im[m].max(); snr=pk/sig
                if tot<=0 or snr<8: r[ln]=None; continue
                xc=(w*c.X).sum()/tot; yc=(w*c.Y).sum()/tot; beam=tb[0]; err=0.6*beam/snr   # centroid error ~ FWHM/(1.7 SNR)
                r[ln]=(xc,yc,err,snr)
            if r['HCN'] and r['HCO+']:
                dx=r['HCN'][0]-r['HCO+'][0]; dy=r['HCN'][1]-r['HCO+'][1]; e=np.hypot(r['HCN'][2],r['HCO+'][2])
                cen.append(dict(ous=ous,beam=tb[0],reg=reg,v=a+20,dx=dx,dy=dy,err=e,HCOx=r['HCO+'][0],HCOy=r['HCO+'][1],snrH=r['HCN'][3],snrO=r['HCO+'][3]))
C=pd.DataFrame(cen); C.to_csv('/workspace/out/q9b/centroids.csv',index=False); print(C.round(3).to_string(index=False))
# (c) injection
print('\n== (c) injection: HCN=R(r_N)*HCO+ (Xb4b native), kinematics identical ==')
O=cube('Xb4b','HCO+'); H=cube('Xb4b','HCN'); x0,y0=centre('Xb4b',(0.3,0.3,0.),'N')
rN=np.hypot(O.X-x0,O.Y-y0); nodes=np.array([0.125,0.375,0.75,1.5,3.0]); Rn=np.array([0.98,0.82,0.63,0.48,0.40])
Rmap=np.where(O.Y>0,np.interp(rN,nodes,Rn),0.38)
inj=O.cube*Rmap[None]      # noise-free model: scaled HCO+ cube (includes its noise x R; negligible for bright bins)
tb=(1.05,1.05,0.); sO=smcube('Xb4b','HCO+',tb); sI=np.empty_like(O.cube)
for i in range(len(O.v)): sI[i]=smooth(inj[i],O,tb)
rows=[]
for r1,r2 in ((0,.5),(0,1.0)):
    xx,yy=centre('Xb4b',tb,'N'); m=apmask(O,xx,yy,r1,r2)
    sp1=apspec(O,sI,m,tb); sp2=apspec(O,sO,m,tb)
    for a,b in zip(EDG[:-1],EDG[1:]): rows.append(dict(r2=r2,v=(a+b)/2,Rinj=winflux(O,sp1,a,b)/winflux(O,sp2,a,b)))
G=pd.DataFrame(rows); print(G.round(3).to_string(index=False))
for r2,g in G.groupby('r2'):
    p=np.polyfit(g.v/100.,g.Rinj,1); print(f'r<{r2}": injection slope {p[0]:.3f} per 100 km/s, mean {g.Rinj.mean():.3f}; (flat intrinsic ratio would give 0)')
G.to_csv('/workspace/out/q9b/injection.csv',index=False)
