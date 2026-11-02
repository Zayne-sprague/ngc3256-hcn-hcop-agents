# Q1 (round 2): HCN/HCO+ 3-2 ratio maps (Xb4d 0.94" with both lines >5 sigma; Xb4b 0.25" with >4 sigma), radial profiles of the ratio in concentric apertures about the N and S
# continuum peaks (cumulative r<=0.25,0.5,1,2" and annuli), cross-array comparison (Xb4b vs Xb4d vs Xb4f) of aperture fluxes for beam-dilution vs filtering, Sblue / outflow region.
# Input: out/mom_<ous>.<spw>_<line>32.npz (s5; |v-VSYS|<350 km/s moment-0, pbcor, Jy/beam km/s, per-pixel noise 'core_err'). Output out/s39_ratio_maps.png/.txt/.json
import sys, json; sys.path.insert(0,'/workspace/scripts')
from common import *; from regions import REGIONS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
K=(267.55763/265.88618)**2
D={}
for ous in ('Xb4b','Xb4d','Xb4f'):
    fn=files(f'*{ous}.ngc_3256_sci.spw25.cube.I.pbcor.fits')[0] if ous!='Xb4f' else files('*Xb4f.ngc_3256_sci.spw16.cube.I.pbcor.fits')[0]
    hd=fits.getheader(fn); wc=WCS(hd).celestial; ps=abs(hd['CDELT1'])*3600
    nh=files(f'*{ous}.ngc_3256_sci.spw27.cube.I.pbcor.fits')[0] if ous!='Xb4f' else files('*Xb4f.ngc_3256_sci.spw18.cube.I.pbcor.fits')[0]
    hh=fits.getheader(nh)
    if ous=='Xb4f': zO=np.load('/workspace/out/mom_Xb4f.spw16_HCO+32.npz'); zH=np.load('/workspace/out/mom_Xb4f.spw18_HCN32.npz')
    else: zO=np.load(f'/workspace/out/mom_{ous}.spw25_HCO+32.npz'); zH=np.load(f'/workspace/out/mom_{ous}.spw27_HCN32.npz')
    D[ous]=dict(wc=wc,ps=ps,x0=int(zO['x0']),y0=int(zO['y0']),O=zO['core'],eO=zO['core_err'],H=zH['core'],eH=zH['core_err'],
                bO=hd['BMAJ']*hd['BMIN']*3600**2,bH=hh['BMAJ']*hh['BMIN']*3600**2)
def px(d,pos):
    x,y=d['wc'].world_to_pixel_values(*pos); return x-d['x0'],y-d['y0']
def apsum(d,mol,pos,r0,r1,coords=None):
    m=d[mol]; n=m.shape[0]; Y,X=np.mgrid[:n,:n]; x,y=px(d,pos); rr=np.hypot(X-x,Y-y)*d['ps']; mk=(rr>=r0)&(rr<r1)&np.isfinite(m)
    ba=1.1331*d['b'+mol]/d['ps']**2; f=m[mk].sum()/ba; e=np.sqrt((d['e'+mol][mk]**2).sum())/ np.sqrt(ba)  # noise: pixels correlated within a beam -> sum err ~ sigma*sqrt(Npix/ba)... conservative sqrt(N_ind)*sigma
    e=d['e'+mol][mk].mean()*np.sqrt(max(mk.sum()/ba,1.0)) if mk.sum()>0 else np.nan
    return f,e
for nm in ('N','S'):
    pos=REGIONS[nm][:2]; P(f'## {nm} nucleus: HCN/HCO+ 3-2 in cumulative apertures (flux Jy km/s; R_Jy; x1.0126 for K)')
    for ous in ('Xb4b','Xb4d','Xb4f'):
        d=D[ous]; row=[]
        for r in (0.25,0.5,1.0,2.0,3.0):
            if ous=='Xb4f' and r<1: continue
            h,eh=apsum(d,'H',pos,0,r); o,eo=apsum(d,'O',pos,0,r); R=h/o; eR=R*np.hypot(eh/h,eo/o) if h>0 and o>0 else np.nan
            row.append((r,h,eh,o,eo,R,eR)); res[f'{nm}_{ous}_r{r}']=dict(HCN=h,eHCN=eh,HCOp=o,eHCOp=eo,R=R,eR=eR)
        P(f'  {ous}: '+' | '.join(f'r<{r}": HCN {h:.2f}+-{eh:.2f} HCO+ {o:.2f}+-{eo:.2f} R={R:.2f}+-{eR:.2f}' for r,h,eh,o,eo,R,eR in row))
    P(f'## {nm}: annuli (Xb4b and Xb4d)')
    for ous in ('Xb4b','Xb4d'):
        d=D[ous]; s=[]
        for r0,r1 in ((0,0.25),(0.25,0.5),(0.5,1.0),(1.0,2.0),(2.0,3.0)):
            h,eh=apsum(d,'H',pos,r0,r1); o,eo=apsum(d,'O',pos,r0,r1); R=h/o if (h>0 and o>0) else np.nan; eR=R*np.hypot(eh/h,eo/o) if R==R else np.nan
            s.append(f'{r0}-{r1}": R={R:.2f}+-{eR:.2f} (HCN {h:.2f}, HCO+ {o:.2f})'); res[f'{nm}_{ous}_ann{r0}-{r1}']=dict(R=R,eR=eR,HCN=h,HCOp=o)
        P(f'  {ous}: '+' | '.join(s))
# centroid offsets of HCN vs HCO+ (Xb4b) near N and S
P('## HCN vs HCO+ peak positions (Xb4b 0.25", within 0.6" of continuum peak) and offset from the continuum peak position (from s1, 0.24" mfs)')
for nm in ('N','S'):
    d=D['Xb4b']; x,y=px(d,REGIONS[nm][:2]); n=d['H'].shape[0]; Y,X=np.mgrid[:n,:n]; mk=np.hypot(X-x,Y-y)*d['ps']<0.6; out=[]
    for mol in ('H','O'):
        m=np.where(mk,d[mol],-1e9); iy,ix=np.unravel_index(np.argmax(m),m.shape)
        # intensity-weighted centroid within 0.3" of peak
        mm=np.hypot(X-ix,Y-iy)*d['ps']<0.3; w=np.clip(d[mol],0,None)*mm; cx_=(w*X).sum()/w.sum(); cy_=(w*Y).sum()/w.sum()
        out.append((mol,(cx_-x)*d['ps'],(cy_-y)*d['ps'],d[mol][iy,ix]))
    P(f'  {nm}: '+'; '.join(f'{m}: centroid offset (pix-x,y)=({-dx:+.2f}" E,{dy:+.2f}" N) peak {pk:.2f} Jy/b km/s' for m,dx,dy,pk in out))
# outflow region vs S aperture
dd=D['Xb4d']; S0=REGIONS['S'][:2]; SB=(S0[0]-0.2/3600/np.cos(np.deg2rad(S0[1])),S0[1]-1.85/3600)
P('## S aperture vs outflow-side region (Xb4d native, |v|<350 km/s integrated line, Jy km/s)')
for nm,pos,r in (('S r<1.0"',S0,1.0),('S r<1.5"',S0,1.5),('Sblue r<1.0"',SB,1.0),('Sblue r<1.5"',SB,1.5),('Sblue r<2.5"',SB,2.5)):
    h,eh=apsum(dd,'H',pos,0,r); o,eo=apsum(dd,'O',pos,0,r); P(f'  {nm:14s} HCN {h:.2f}+-{eh:.2f} HCO+ {o:.2f}+-{eo:.2f} R={h/o:.2f}+-{h/o*np.hypot(eh/h,eo/o):.2f}')
# maps
fig,ax=plt.subplots(1,3,figsize=(16,5.2))
for a,(ous,half,nsig,pos) in zip(ax,(('Xb4d',6.0,5,CEN:=[0.5*(NUC['N'][i]+NUC['S'][i]) for i in (0,1)]),('Xb4b',1.6,4,NUC['N']),('Xb4b',1.6,4,NUC['S']))):
    d=D[ous]; x,y=px(d,pos); n=d['H'].shape[0]; ok=(d['H']>nsig*d['eH'])&(d['O']>nsig*d['eO']); Rm=np.where(ok,d['H']/d['O'],np.nan)
    w=half/d['ps']; x0,x1,y0,y1=int(x-w),int(x+w),int(y-w),int(y+w)
    im=a.imshow(Rm[y0:y1,x0:x1],origin='lower',vmin=0.2,vmax=1.4,cmap='viridis',extent=[-half,half,-half,half]); a.invert_xaxis()
    a.contour(np.linspace(half,-half,x1-x0),np.linspace(-half,half,y1-y0),d['O'][y0:y1,x0:x1],levels=np.array([0.1,0.3,0.6,0.9])*np.nanmax(d['O'][y0:y1,x0:x1]),colors='w',linewidths=.6)
    a.set_title(f'{ous} HCN/HCO+ 3-2 (Jy ratio), both >{nsig}$\\sigma$'); a.set_xlabel('dE [arcsec]'); plt.colorbar(im,ax=a)
plt.tight_layout(); plt.savefig('/workspace/out/s39_ratio_maps.png',dpi=85)
# radial/azimuthal sectors around S in Xb4d: ratio vs azimuth in annulus 0.5-2"
P('## Xb4d ratio vs azimuth in the annulus 0.5-2" about S and N (4 quadrants; E of N=90)')
for nm in ('N','S'):
    d=dd; x,y=px(d,REGIONS[nm][:2]); n=d['H'].shape[0]; Y,X=np.mgrid[:n,:n]; rr=np.hypot(X-x,Y-y)*d['ps']; ang=np.degrees(np.arctan2(-(X-x),(Y-y)))%360; s=[]
    for lab,(a0,a1) in {'N':(315,45),'E':(45,135),'S':(135,225),'W':(225,315)}.items():
        sec=((ang>=a0)|(ang<a1)) if a0>a1 else ((ang>=a0)&(ang<a1)); mk=sec&(rr>=0.5)&(rr<2.0)
        ba=1.1331*d['bH']/d['ps']**2; h=d['H'][mk].sum()/ba; o=d['O'][mk].sum()/ba; s.append(f'{lab}: R={h/o:.2f} (HCN {h:.2f},HCO+ {o:.2f})')
    P(f'  {nm}: '+' | '.join(s))
json.dump(res,open('/workspace/out/s39_ratio_maps.json','w'),indent=1,default=float); open('/workspace/out/s39_ratio_maps.txt','w').write('\n'.join(txt))
