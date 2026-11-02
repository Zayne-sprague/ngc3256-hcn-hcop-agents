# Q1/Q6: integrated flux in circular apertures (centred on the field centre) for Xb4b/Xb4d/Xb4f (native resolution; total flux is beam independent) -> resolved-out flux.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
rows=[]
def apflux(img,d,R,sig):
    r=np.hypot(d['dx'],d['dy']); m=(r<R)&np.isfinite(img); npb=OM(d['beam'])/d['pix']**2
    return np.nansum(img[m])/npb, sig*np.sqrt(m.sum()/npb)     # correlated noise: sigma_pix*sqrt(N_ap/N_pix_per_beam)
for sp in ('HCN','HCO+'):
    for ous in ('Xb4f','Xb4d','Xb4b'):
        k=f'{sp} 3-2@{ous}'; d=load(k,crop=11.0 if ous=='Xb4b' else None); pr=noise_profile(d,None)
        sig=np.nanmedian([s for _,_,s,_,_ in pr[:2]])
        for R in (3,5,8,10):
            F,e=apflux(d['m0'],d,R,sig); rows.append(dict(what=k.split('@')[0]+' (Jy km/s)',ous=ous,R=R,F=F,eF=e))
for ous,f_ in (('Xb4f','member.uid___A001_X1273_Xb4f.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'),('Xb4d','member.uid___A001_X1273_Xb4d.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits'),('Xb4b','member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')):
    h=fits.open('/data/'+f_); hd=h[0].header; im=np.squeeze(h[0].data).astype(float); pix=abs(hd['CDELT2'])*3600
    w=WCS(hd).celestial; ny,nx=im.shape; Y,X=np.mgrid[:ny,:nx]; ra,de=w.pixel_to_world_values(X,Y); dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600
    d=dict(dx=dx,dy=dy,beam=(hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']),pix=pix); r=np.hypot(dx,dy)
    sig=1.4826*np.nanmedian(np.abs(im[(r>12)&(r<20)]-np.nanmedian(im[(r>12)&(r<20)])))
    for R in (3,5,8,10):
        F,e=apflux(im,d,R,sig); rows.append(dict(what='1.1mm continuum (Jy)',ous=ous,R=R,F=F,eF=e))
T=pd.DataFrame(rows); P=T.pivot_table(index=['what','R'],columns='ous',values=['F','eF'])
out=[]
for (w,R),r in P.iterrows():
    for a in ('Xb4d','Xb4b'):
        q=r['F'][a]/r['F']['Xb4f']; e=abs(q)*np.hypot(r['eF'][a]/r['F'][a],r['eF']['Xb4f']/r['F']['Xb4f'])
        out.append(dict(what=w,R=R,F_Xb4f=r['F']['Xb4f'],F_Xb4d=r['F']['Xb4d'],F_Xb4b=r['F']['Xb4b'],comp=f'{a}/Xb4f',frac=q,e_frac=e))
O=pd.DataFrame(out); O.to_csv('/workspace/out/q1/q1_flux_recovery.csv',index=False); T.to_csv('/workspace/out/q1/q1_aperture_fluxes.csv',index=False)
pd.set_option('display.width',250); print(O.round(3).to_string())
