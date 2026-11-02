# Q13g: global (post-selection) significance of the S-lobe E = HCN - k*HCO+ (Xb4d, 1.05" beam, r=1" aperture, 2375-2575 km/s).
#  (1) blind search: E map (disk-mean r=1") over the whole field (r<=R), rank of the lobe peak among all local maxima; same for red mirror window
#  (2) null from sliding 200 km/s windows in line-free channels (HCN 1380-2100, HCO+ own-frame free ranges): max-over-field of E per window, with k fixed 0.574
#  (3) matched-environment: ring apertures at 4-5" from S (and N) at all azimuths, E and red mirror, plus sliding windows on the same ring
# Outputs: out/q13/q13g.log, q13g_*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve
log=open(O+'q13g.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
rob=lambda a: 1.4826*np.median(np.abs(a-np.median(a)))
K=0.574; RAD=1.0; RMAX=7.5
Mx=wmaps('Xb4d'); dx,dy,dv,pix=Mx['dx'],Mx['dy'],Mx['dv'],Mx['pix']
vh,vo=Mx['v_HCN'],Mx['v_HCO+']; H,Ob=Mx['HCN_spec'],Mx['HCO+_spec']
n=int(np.ceil(RAD/pix)); yy,xx=np.mgrid[-n:n+1,-n:n+1]; ker=(np.hypot(xx,yy)*pix<=RAD).astype(float); ker/=ker.sum()
sm=lambda M: fftconvolve(M,ker,mode='same')
field=np.hypot(dx,dy)<=RMAX
def win(spec,v,a,b): m=(v>=a)&(v<b); return spec[m].sum(0)*dv
def Emap(a,b,k=K): return sm(win(H,vh,a,b)-k*win(Ob,vo,a,b))
def peaks(E,mask,thr=None):
    mx=(E==maximum_filter(E,size=2*n+1))&mask
    return E[mx],dx[mx],dy[mx]
x0,y0=POS['Soutflow']; lob=np.hypot(dx-x0,dy-y0)<=0.3
Eb=Emap(2375,2575); Er=Emap(2975,3175)
E0=Eb[lob].max(); pr(f'lobe E (r=1", k={K}) max within 0.3" of Soutflow: {E0:.4f}; at nominal pos: {Eb.ravel()[np.argmin(np.hypot(dx-x0,dy-y0))]:.4f}')
# (1) blind search
for nm,E in (('E_blue',Eb),('E_red(mirror)',Er)):
    for excl in (0,1.5):
        mask=field.copy()
        if excl:
            for p in ('N','S'): mask&=np.hypot(dx-POS[p][0],dy-POS[p][1])>excl
        pk,px,py=peaks(E,mask); o=np.argsort(-pk)
        pr(f'{nm}: local maxima (field r<={RMAX}", excl {excl}" around N,S): n={len(pk)}; top5 (x,y,E):',[(round(px[i],1),round(py[i],1),round(pk[i],3)) for i in o[:5]])
        if nm=='E_blue': pr(f'   rank of lobe among maxima: {(pk>E0+1e-9).sum()+1} of {len(pk)}')
sE=rob(Eb[field]); sEr=rob(Er[field]); pr(f'robust sigma of E-blue map over field (r<={RMAX}", includes real structure): {sE:.4f}; red mirror {sEr:.4f}; lobe z={E0/sE:.2f}')
from scipy.stats import norm
Nind=np.pi*RMAX**2/(1.1331*1.05**2)
for nm,sg in (('map-MAD E-blue',sE),('random-aperture null (Q13a) 0.0270',0.0270),('map-MAD E-red',sEr)):
    z=E0/sg; p1=norm.sf(z); pg=1-(1-p1)**Nind
    pr(f'  [{nm}] z={z:.2f} p_local={p1:.2e}; N_ind={Nind:.0f}: p_global={pg:.2e} = {norm.isf(pg):.2f} sigma (one-sided)')
# (2) sliding windows null: 200 km/s, start step 25, both HCN and HCO+ free
def free(v,rngs,a,b): return any(a>=r0 and b<=r1 for r0,r1 in rngs)
rows=[]
for a in np.arange(1380,1900,25):
    b=a+200
    if not (free(vh,FREE['HCN'],a,b) and free(vo,FREE['HCO+'],a,b)): continue
    E=Emap(a,b)
    for excl in (0,1.5):
        mask=field.copy()
        if excl:
            for p in ('N','S'): mask&=np.hypot(dx-POS[p][0],dy-POS[p][1])>excl
        pk,_,_=peaks(E,mask)
        rows.append(dict(a=a,excl=excl,maxE=pk.max(),sig=rob(E[mask]),nmax=len(pk),n_ge=int((pk>=E0).sum())))
W=pd.DataFrame(rows); W.round(4).to_csv(O+'q13g_sliding.csv',index=False)
pr('sliding 200-km/s line-free windows (HCN 1380-2100, HCO+ free): max over field of E(r=1")'); pr(W.round(4).to_string(index=False))
pr('rms (robust) of E over field, line-free windows: %.4f ; lobe E / that = %.2f'%(W[W.excl==0].sig.median(),E0/W[W.excl==0].sig.median()))
# independent windows: starts separated by >=200 -> at most 1 in 1380-1680? use non-overlapping subset
ind=W[(W.excl==0)&(W.a.isin([1380,1480,1880]))]
pr('non-overlapping/near-independent windows (starts 1380,1480*,1880): max E',ind.maxE.round(3).tolist(),' (1480 overlaps 1380 by 100)')
# pooled Gaussian-equivalent trials: N_ind positions
Nind=np.pi*RMAX**2/(1.13*1.05**2*np.pi/ (np.pi))  # area/(beam area for FWHM 1.05: 1.1331*theta^2)
Nind=np.pi*RMAX**2/(1.1331*1.05**2); pr(f'N_ind positions (field r<={RMAX}", beam area 1.133*1.05^2) = {Nind:.0f}')
from scipy.stats import norm
sig=W[W.excl==0].sig.median()
z=E0/sig; p1=norm.sf(z); pg=1-(1-p1)**Nind
pr(f'Gaussian: local z={z:.2f}, p_local={p1:.2e}; global (N_ind={Nind:.0f}) p={pg:.2e} -> {norm.isf(pg/2) if pg<1 else 0:.2f} sigma two-sided-equivalent, one-sided {norm.isf(pg):.2f}')
# empirical sig for non-Gaussian tails: use all field E values in sliding windows, pooled z excess
Z=[]
for a in W[W.excl==0].a:
    E=Emap(a,a+200); Z.append(E[field])
Z=np.concatenate(Z); pr(f'pooled field E values in line-free windows: n={Z.size}; robust sigma {rob(Z):.4f}; fraction >= lobe {np.mean(Z>=E0):.2e}; max {Z.max():.3f}; fraction <= -lobe {np.mean(Z<=-E0):.2e}')
# (3) ring null
rows=[]; phis=np.deg2rad(np.arange(0,360,10))
def ring(center,R,E):
    out=[]
    for ph in phis:
        x=center[0]+R*np.sin(ph); y=center[1]+R*np.cos(ph)   # azimuth from north
        if np.hypot(x,y)>RMAX+1: out.append(np.nan); continue
        i=np.argmin(np.hypot(dx-x,dy-y)); out.append(E.ravel()[i])
    return np.array(out)
R_lobe=np.hypot(x0-POS['S'][0],y0-POS['S'][1]); pr(f'lobe distance from S: {R_lobe:.2f}"; azimuth from N {np.rad2deg(np.arctan2(x0-POS["S"][0],y0-POS["S"][1])):.0f} deg')
ringdata={}
for cname in ('S','N'):
    for R in (R_lobe,):
        for nm,(a,b) in (('blue',(2375,2575)),('red',(2975,3175))):
            r=ring(POS[cname],R,Emap(a,b)); ringdata[(cname,nm)]=r
            pr(f'ring around {cname} R={R:.2f}" {nm} window: n={np.isfinite(r).sum()} median {np.nanmedian(r):.3f} robust-sd {rob(r[np.isfinite(r)]):.3f} max {np.nanmax(r):.3f} min {np.nanmin(r):.3f}; #>=lobe {np.sum(r>=E0)}')
# ring on sliding windows: pooled
for cname in ('S','N'):
    pool=np.concatenate([ring(POS[cname],R_lobe,Emap(a,a+200)) for a in W[W.excl==0].a]); pool=pool[np.isfinite(pool)]
    pr(f'ring {cname} R={R_lobe:.2f}" line-free sliding windows pooled: n={pool.size} sd {pool.std():.4f} robust {rob(pool):.4f} max {pool.max():.3f} min {pool.min():.3f}; frac>=lobe {np.mean(pool>=E0):.3e}; lobe z vs ring sd {E0/pool.std():.2f}')
    # same, windows near the lobe velocity from adjacent unused? also wider ring 4-5"
for cname in ('S','N'):
    vals=[]
    for R in (4.0,4.5,5.0):
        for a in W[W.excl==0].a: vals.append(ring(POS[cname],R,Emap(a,a+200)))
    vals=np.concatenate(vals); vals=vals[np.isfinite(vals)]
    pr(f'ring {cname} 4-5" sliding pooled: n={vals.size} sd {vals.std():.4f} max {vals.max():.3f}; frac>=lobe {np.mean(vals>=E0):.2e}')
# peak search on ring-region (S annulus 3-6", all azimuth) for line-free windows: max per window
for a in W[W.excl==0].a:
    E=Emap(a,a+200); r=np.hypot(dx-POS['S'][0],dy-POS['S'][1]); m=(r>=3.5)&(r<=5.5)&field
    pk,_,_=peaks(E,m); rows.append(dict(a=a,maxE_Sannulus=pk.max(),sd=E[m].std()))
R=pd.DataFrame(rows); R.round(4).to_csv(O+'q13g_Sannulus.csv',index=False)
pr('S-annulus (3.5-5.5") max E per sliding window:',R.maxE_Sannulus.round(3).tolist(),' sd(E) in annulus:',R.sd.round(4).mean())
# ring excluding the lobe azimuth (+-45 deg), blue window and red mirror, S ring R=1.93"
phl=np.rad2deg(np.arctan2(x0-POS['S'][0],y0-POS['S'][1]))
for nm,E in (('blue',Eb),('red',Er)):
    r=ringdata[('S',nm)]; az=np.rad2deg(phis); d=np.abs((az-phl+180)%360-180); m=d>45
    pr(f'S ring {nm} excluding lobe +-45deg: n={m.sum()} mean {r[m].mean():.3f} sd {r[m].std():.3f} max {r[m].max():.3f}; lobe z vs this sd: {(E0-r[m].mean())/r[m].std():.1f}')
