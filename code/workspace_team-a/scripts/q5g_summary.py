# Q5g: (i) sizes of the matched sub-clumps (HCN, HCO+), (ii) line-continuum offsets with Condon errors, (iii) virial masses, N/S ratio, crossing times.
# python3 -I scripts/q5g_summary.py -> out/q5/q5g.log, q5g_clumps.csv, q5g_offsets.csv, q5g_virial.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
import pandas as pd
from scipy import ndimage as ndi
d=load(); dx=d['dx']; dy=d['dy']; log=open(O+'q5g.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
mH,nH=win_map(d,'H',*VW); mO,nO=win_map(d,'O',*VW)
Bc={'HCN':beamcov(*d['beamH']),'HCO+':beamcov(*d['beamO'])}; BM={'HCN':d['beamH'],'HCO+':d['beamO']}; IMG={'HCN':mH,'HCO+':mO}
def condon(F,fa,fb,bm,pk,rms):
    bmaj,bmin=bm[0],bm[1]; cm,cn=np.sqrt(fa**2+bmaj**2),np.sqrt(fb**2+bmin**2); k=cm*cn/(4*bmaj*bmin); snr=pk/rms
    rho=lambda a,b: np.sqrt(k*(1+bmaj**2/cm**2)**a*(1+bmin**2/cn**2)**b)*snr
    return cm/2.3548*np.sqrt(2)/rho(1.5,1.5), (cm/fa*np.sqrt(2)/rho(2.5,0.5)*cm if fa>0.03 else np.nan),(cn/fb*np.sqrt(2)/rho(0.5,2.5)*cn if fb>0.03 else np.nan)
def fit1(nm,x,y,radius):
    img=IMG[nm]; mask=np.hypot(dx-x,dy-y)<radius; Bm=Bc[nm]
    z=np.where(mask,ndi.gaussian_filter(img,2),-1e9); iy,ix=np.unravel_index(np.argmax(z),z.shape)
    p=[img[iy,ix]*1.3,dx[iy,ix],dy[iy,ix],0.1,0.1,0.,0,0,0]; lo=[0,x-radius,y-radius,0.005,0.005,-400,-9,-9,-9]; hi=[50,x+radius,y+radius,1.0,1.0,400,9,9,9]
    r=fit(img,dx,dy,mask,np.array(p),np.array(lo),np.array(hi),Bm_ if False else Bm,1,True); par=r.x
    res=img-model(par,dx,dy,Bm,1,True); rms=np.std(res[mask]); fa,fb=decon(par); pk=model(par[:6],dx,dy,Bm,1,False).max()
    ep,ea,eb=condon(par[0],fa,fb,BM[nm],pk,rms)
    return dict(flux=par[0],x=par[1],y=par[2],fwhm_maj=fa,fwhm_min=fb,pa=(par[5] if par[3]>=par[4] else par[5]+90)%180,peak=pk,snr=pk/rms,e_pos=ep,e_maj=ea,e_min=eb)
# (i) clumps: matched hp maxima from q5c (positions rounded) + the three main sources
CL={'N core':(0.0,2.63),'N-SW':(-0.55,1.9),'N-E':(1.24,2.13),'Br-E':(1.12,1.29),'Br-S':(0.37,0.14),'S core':(-0.15,-2.60),'S-E':(1.20,-2.52),'S-W2':(0.60,-2.75),'Sblue':(-0.2,-4.5),'NNW':(-1.44,5.07)}
rows=[]
for cn_,(x,y) in CL.items():
    for nm in ('HCN','HCO+'):
        for rad in (0.35,0.5):
            r=fit1(nm,x,y,rad); rows.append(dict(clump=cn_,line=nm,radius=rad,**r))
cd=pd.DataFrame(rows); cd.to_csv(O+'q5g_clumps.csv',index=False); pd.set_option('display.width',220)
pr('# sub-clump Gaussian fits (single Gaussian + plane in r<radius; map = +-350 km/s moment 0; FWHM intrinsic (beam-deconvolved) in arcsec, 1"=194 pc)')
pr(cd.round(3).to_string())
# (ii) offsets (from q5d_fits.csv)
q=pd.read_csv(O+'q5d_fits.csv'); rows=[]
def pick(src,mp,rad,ncomp=1,comp=0):
    r=q[(q.src==src)&(q['map']==mp)&(q.radius==rad)&(q.ncomp==ncomp)&(q.comp==comp)]; return r.iloc[0] if len(r) else None
for src,rad,nc in (('S',0.8,1),('S',1.4,2),('N',2.0,1),('N',0.6,1),('Bridge',1.0,1)):
    c=pick(src,'cont',rad,nc)
    if c is None: continue
    for mp in ('HCN','HCO+'):
        l=pick(src,mp,rad,nc)
        if l is None: continue
        ddx,ddy=l.x-c.x,l.y-c.y; sep=np.hypot(ddx,ddy); e=np.hypot(l.e_pos,c.e_pos); bm=0.24
        rows.append(dict(src=src,radius=rad,ncomp=nc,line=mp,d_east=ddx,d_north=ddy,sep=sep,e_sep=e,sep_beam=sep/bm,e_beam=e/bm,sep_pc=sep*PC,sig=sep/e))
od=pd.DataFrame(rows); od.to_csv(O+'q5g_offsets.csv',index=False)
pr('\n# line minus 1.1 mm continuum centroid offsets (same-array, same phase centre; Condon errors with fit-residual noise, conservative)')
pr(od.round(3).to_string())
# (iii) virial masses: M = 1040 sigma^2 R, R=1.91 sigma_r (Rosolowsky & Leroy 2006), sigma_r = geometric-mean FWHM/2.3548 of the deconvolved Gaussian
w=pd.read_csv(O+'q5e_width.csv')
rows=[]
for src,var in (('N',[('HCN',0.6,1,0),('HCO+',0.6,1,0),('HCN',2.0,1,0),('HCN',2.0,2,0)]),('S',[('HCN',0.8,1,0),('HCO+',0.8,1,0),('HCO+',1.4,2,0),('HCN',1.4,1,0)]),('Bridge',[('HCN',0.8,1,0),('HCO+',0.8,1,0)])):
    for nm,rad,nc,ic in var:
        r=pick(src,nm,rad,nc,ic)
        if r is None: continue
        fg=np.sqrt(r.fwhm_maj*max(r.fwhm_min,0.02)); sr=fg/2.3548*PC; R=1.91*sr
        for ap in (0.0,0.12):
            ww=w[(w.line==nm)&(w.src==src)&(w.aperture_r==ap)]
            if not len(ww): continue
            sg=ww.sigma.iloc[0]
            # beam smearing of the planar gradient: sigma_smear = G*theta/2.355 with G from q5e_grad.csv (R=0.6, same line)
            gdf=pd.read_csv(O+'q5e_grad.csv'); gg=gdf[(gdf.line==nm)&(gdf.src==src)&(gdf.R==0.6)].grad.iloc[0]
            sm=gg*0.24/2.3548; sint=np.sqrt(max(sg**2-sm**2,1.))
            rows.append(dict(src=src,line=nm,fit=f'r<{rad}"{"" if nc==1 else " 2comp"}',fwhm_geo_arcsec=fg,R_pc=R,aperture=('peak pixel' if ap==0 else f'r={ap}"'),sigma_obs=sg,sigma_smear=sm,sigma_int=sint,Mvir=1040*sint**2*R,Mvir_obs=1040*sg**2*R,Sigma_vir=1040*sint**2*R/(np.pi*R**2),tcross_Myr=R/sint*0.978))
vd=pd.DataFrame(rows); vd.to_csv(O+'q5g_virial.csv',index=False)
pr('\n# virial masses M=1040 sigma^2 R_pc (R=1.91*sigma_r, sigma_r from the deconvolved geometric FWHM); tcross=R/sigma [Myr]')
pr(vd.round(3).to_string(float_format=lambda x:f'{x:.3g}'))
