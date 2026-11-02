# Q2 step c: wing fluxes, positions (2D Gaussian fits with noise-injection errors), velocity extent, and wing kinematics vs offset along the S outflow axis.
# Conventions: Xb4d 1.05" circular (flux-rescaled convolution), also Xb4b smoothed to 1.05" (recovery) ; pbcor; baseline = per-pixel median v-VSYS in -900..-450;
# windows about VSYS=2780 radio: blue -405..-150, red 150..405 (red candidate 180..320). Fluxes are Jy km/s in r=1" apertures (beam-sum / beam area => Jy km/s of the aperture)
# and 2D-Gaussian integrated fluxes of the moment-0 wing maps (Jy/beam km/s * area/beam area). Noise from random line-free apertures / injected noise cut-outs.
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
from scipy.optimize import curve_fit
OUT='/workspace/out/q2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
S0=REGIONS['S'][:2]; N0=REGIONS['N'][:2]; rng=np.random.default_rng(7)
BLUE=(-405,-150); RED=(150,405); REDC=(180,320)
res={}
def prep(ous,fw,half):
    o={}
    for mol,spw in (('HCN','spw27'),('HCO+','spw25')):
        c,v,ps,wc,x0,y0,b0=load_cube(ous,spw,F0[mol][3],half=half)
        c=conv(c,ps,b0,(fw,fw,0.0)); c=c-np.median(c[(v>-900)&(v<-450)],axis=0)
        o[mol]=dict(c=c,v=v,dv=abs(np.median(np.diff(v))),ps=ps,n=c.shape[1],beam=(fw,fw,0.0),ba=1.1331*fw*fw/ps**2,pxy=pxy_fn(wc,x0,y0))
    return o
def mom(R,a,b): k=(R['v']>=a)&(R['v']<=b); return R['c'][k].sum(0)*R['dv']
def field_mask(R,rmin=3.5,rmax=None):
    # not-bright pixels r>rmin from the N-S midpoint; rmax default 0.45*image width
    m0=mom(R,-350,350); sg=np.median(np.abs(m0-np.median(m0)))*1.4826; n=R['n']; ps=R['ps']
    br=binary_dilation(m0>4*sg,iterations=int(1.5/ps)); Y,X=np.mgrid[:n,:n]; cx,cy=R['pxy'](*CEN)
    rc=np.hypot(X-cx,Y-cy)*ps; return (~br)&(rc>rmin)&(rc<(rmax or n*ps*0.47))
def gauss2d(xy,A,x0,y0,sx,sy,th):
    x,y=xy; c,s=np.cos(th),np.sin(th); u=(x-x0)*c+(y-y0)*s; w=-(x-x0)*s+(y-y0)*c
    return A*np.exp(-0.5*((u/sx)**2+(w/sy)**2))
def fit2d(M,ps,pos_pix,beam,half=2.5,p0sig=1.0):
    """fit elliptical Gaussian (already beam-convolved, free widths) to cut-out around pos_pix (pixels); returns params (arcsec offsets handled by caller)"""
    r=int(half/ps); x0=int(round(pos_pix[0])); y0=int(round(pos_pix[1])); sub=M[y0-r:y0+r+1,x0-r:x0+r+1]
    Y,X=np.mgrid[y0-r:y0+r+1,x0-r:x0+r+1]; pk=np.unravel_index(np.argmax(sub),sub.shape)
    p0=[sub.max(),X[pk],Y[pk],0.5*beam/ps,0.5*beam/ps,0.3]
    lo=[0,x0-r,y0-r,0.35*beam/ps,0.35*beam/ps,-np.pi]; hi=[10,x0+r,y0+r,4*beam/ps,4*beam/ps,np.pi]
    p,_=curve_fit(gauss2d,(X.ravel(),Y.ravel()),sub.ravel(),p0=p0,bounds=(lo,hi),maxfev=20000); return p
def gfit_summary(p,ps,beam,S_pix):
    A,x0,y0,sx,sy,th=p; flux_bm=A*2*np.pi*sx*sy/(1.1331*beam**2/ps**2)   # Jy km/s (A in Jy/b km/s)
    dE=-(x0-S_pix[0])*ps; dN=(y0-S_pix[1])*ps; return dict(A=A,dE=dE,dN=dN,fwhm_a=2.355*max(sx,sy)*ps,fwhm_b=2.355*min(sx,sy)*ps,flux=flux_bm)
def fit_with_errors(R,M,pos_pix,S_pix,nmc=30):
    p=fit2d(M,R['ps'],pos_pix,R['beam'][0]); base=gfit_summary(p,R['ps'],R['beam'][0],S_pix)
    fm=field_mask(R); ys,xs=np.where(fm); r=int(2.5/R['ps']); draws=[]; rr=int(r)
    Y,X=np.mgrid[-r:r+1,-r:r+1]; model=gauss2d((X.ravel()+0.0,Y.ravel()+0.0),p[0],p[1]-int(round(pos_pix[0])),p[2]-int(round(pos_pix[1])),p[3],p[4],p[5]).reshape(X.shape)
    tries=0
    while len(draws)<nmc and tries<2000:
        tries+=1; i=rng.integers(len(ys)); cy,cx=ys[i],xs[i]
        if cy-r<0 or cx-r<0 or cy+r+1>R['n'] or cx+r+1>R['n']: continue
        if fm[cy-r:cy+r+1,cx-r:cx+r+1].mean()<0.6: continue
        sub=M[cy-r:cy+r+1,cx-r:cx+r+1]+model
        try:
            Yg,Xg=np.mgrid[cy-r:cy+r+1,cx-r:cx+r+1]
            q,_=curve_fit(gauss2d,(Xg.ravel(),Yg.ravel()),sub.ravel(),p0=[p[0],cx,cy,p[3],p[4],p[5]],bounds=([0,cx-r,cy-r,0.35*R['beam'][0]/R['ps'],0.35*R['beam'][0]/R['ps'],-np.pi],[10,cx+r,cy+r,4*R['beam'][0]/R['ps'],4*R['beam'][0]/R['ps'],np.pi]),maxfev=5000)
            g=gfit_summary(q,R['ps'],R['beam'][0],(cx-(p[1]-S_pix[0]),cy-(p[2]-S_pix[1]))); draws.append([g['A'],g['dE'],g['dN'],g['flux']])
        except Exception: pass
    d=np.array(draws); base['err']=dict(A=float(d[:,0].std()),dE=float(d[:,1].std()),dN=float(d[:,2].std()),flux=float(d[:,3].std()),n=len(d)); return base
SET={'Xb4d':prep('Xb4d',1.05,24.0),'Xb4b':prep('Xb4b',1.05,10.0)}
# --- (1) wing maps, 2D Gaussian fits (Xb4d)
guess={'blue':(0.0,-1.9),'red':(0.15,+3.3)}
for ous in ('Xb4d','Xb4b'):
  for mol in ('HCN','HCO+'):
    R=SET[ous][mol]; Spix=np.array(R['pxy'](*S0))
    for nm,(a,b),gp in (('blue',BLUE,guess['blue']),('red',RED,guess['red']),('redC',REDC,guess['red'])):
        M=mom(R,a,b); pp=(Spix[0]-gp[0]/R['ps'],Spix[1]+gp[1]/R['ps'])
        try: g=fit_with_errors(R,M,pp,Spix,nmc=30 if ous=='Xb4d' else 15)
        except Exception as ex: P(f'{ous} {mol} {nm}: fit failed {ex}'); continue
        e=g['err']; res[f'fit|{ous}|{mol}|{nm}']=g
        P(f'FIT {ous}@1.05 {mol} {nm}({a}..{b}): peak {g["A"]:.3f}+-{e["A"]:.3f} Jy/b km/s ({g["A"]/e["A"]:.1f} sigma) at dE={g["dE"]:+.2f}+-{e["dE"]:.2f}" dN={g["dN"]:+.2f}+-{e["dN"]:.2f}" from S; FWHM {g["fwhm_a"]:.2f}x{g["fwhm_b"]:.2f}" (convolved); integrated {g["flux"]:.3f}+-{e["flux"]:.3f} Jy km/s')
# --- (2) aperture fluxes r=1" at blue lobe and red candidate, noise from 60 random apertures; Xb4d, Xb4b, recovery
apos={'Slobe':offset_pos(S0,0.15,-1.85),'Slobe_r1.5':offset_pos(S0,0.15,-1.85),'Red':offset_pos(S0,0.15,3.3),'S':S0,'N':N0}
for ous in ('Xb4d','Xb4b'):
    for mol in ('HCN','HCO+'):
        R=SET[ous][mol]; fm=field_mask(R); ys,xs=np.where(fm)
        for an,(a,b) in (('blue',BLUE),('red',RED),('redC',REDC),('core',(-205,95))):
            for pn,rad in (('Slobe',1.0),('Slobe_r1.5',1.5),('Red',1.0),('S',1.0),('N',1.0)):
                if an not in ('core',) and pn in ('S','N'): continue
                sp=aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],R['pxy'](*apos[pn]),rad)); f=W(R,sp,a,b)
                # noise: random apertures in the field
                nn=[]; Y,X=np.mgrid[:R['n'],:R['n']]
                for _ in range(80):
                    i=rng.integers(len(ys)); m=np.hypot(X-xs[i],Y-ys[i])*R['ps']<=rad
                    if not fm[m].all(): continue
                    nn.append(W(R,aspec(R['c'],R['ba'],m),a,b))
                s=np.std(nn); res[f'ap|{ous}|{mol}|{an}|{pn}']=dict(F=float(f),sig=float(s),n=len(nn))
                P(f'AP {ous}@1.05 {mol} {an}({a}..{b}) {pn} r={rad}": {f:+.3f}+-{s:.3f} Jy km/s ({f/s:+.1f} sigma, {len(nn)} random apertures)')
# --- (3) velocity extent of the blue lobe spectrum (r=1" aperture, rebinned 5 channels = 22 km/s), Xb4d
def edge(v,sp,e,nb=5,thr=2.0,side='blue'):
    o=np.argsort(v); v=v[o]; sp=sp[o]
    k=len(v)//nb*nb; vb=v[:k].reshape(-1,nb).mean(1); sb=sp[:k].reshape(-1,nb).mean(1); eb=e/np.sqrt(nb)
    ok=sb>thr*eb
    idx=np.where(ok&(vb<-150))[0]
    # contiguous run from v=-150 blueward
    j=np.where((vb<-100))[0][-1]; ext=None
    while j>=0 and sb[j]>thr*eb: ext=vb[j]; j-=1
    return ext
for mol in ('HCN','HCO+'):
    R=SET['Xb4d'][mol]; sp=aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],R['pxy'](*apos['Slobe']),1.0)); lf=(R['v']>-900)&(R['v']<-450); e=sp[lf].std()
    for thr in (2.0,3.0): P(f'EDGE Xb4d {mol} Slobe r=1": contiguous blue emission >{thr} sigma (22 km/s bins, rms {e*1e3:.2f} mJy/ch) reaches v-VSYS = {edge(R["v"],sp,e,thr=thr)} km/s')
    # cumulative-flux weighted velocity of the wing (v<-150): flux-weighted mean and 95% blue limit
    k=(R['v']>=-500)&(R['v']<=-150); w=np.clip(sp[k],0,None); vm=(w*R['v'][k]).sum()/w.sum()
    P(f'MEANV Xb4d {mol} Slobe r=1": flux-weighted mean velocity of the -500..-150 wing = {vm:.0f} km/s; total flux -500..-150 = {sp[k].sum()*R["dv"]:.3f} Jy km/s')
    res[f'meanv|{mol}']=float(vm)
# --- (4) wing kinematics vs offset along PA185 (first moment of -420..-100 above 2 sigma) from Xb4d PV
d=np.load(OUT+'/s75a_pv.npz'); offs=d['offs']
for mol in ('HCN','HCO+'):
    M=d[f'Xb4d_1.05_{mol}']; v=d['v_'+mol]; lf=(v>-900)&(v<-450); s=M[lf].std(); dv=abs(np.median(np.diff(v))); k=(v>=-420)&(v<=-100)
    row=[]
    for j,o in enumerate(offs):
        if o>0.1 or o<-4.1: continue
        F=M[k,j].sum()*dv; eF=s*np.sqrt(k.sum())*dv; w=np.clip(M[k,j],0,None); vm=(w*v[k]).sum()/w.sum() if w.sum()>0 else np.nan
        tot=M[(v>-420)&(v<420),j].sum()*dv; row.append((o,F,eF,vm,F/tot if tot>0 else np.nan))
    P(f'KIN {mol} Xb4d 1.05" slit 1" along PA185 (offset N of S; wing -420..-100): '+' | '.join(f'{o:+.2f}": F={F:.3f}+-{eF:.3f} v={vm:.0f} wing/tot={fr:.2f}' for o,F,eF,vm,fr in row))
    res[f'kin|{mol}']=row
json.dump(res,open(OUT+'/s75b.json','w'),indent=1,default=float); open(OUT+'/s75b.txt','w').write('\n'.join(txt))
