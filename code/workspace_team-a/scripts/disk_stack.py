# Stacked spectra over the objectively defined extended-disk region (Xb55 HCO+ 2-1 island 1 at r>8", S/N>5 after bowl removal; see disk.py) in the Xb55, Xb47, Xb4f cubes.
# Integrated intensities (+-350 km/s, linear baseline) with empirical noise from the same-size random masks in emission-free pixels (r>16"), CH3OH Q-branch by template fit.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from scipy import ndimage as ndi
import pandas as pd
from ch3oh_matched import fit_one
from ch3oh import toK
RA0,DE0=156.96346,-43.90461
hB=fits.getheader(f('Xb55','spw18')); wB=WCS(hB).celestial; ny,nx=hB['NAXIS2'],hB['NAXIS1']; Y,X=np.mgrid[:ny,:nx]
ra,de=wB.pixel_to_world_values(X,Y); dx=(ra-RA0)*np.cos(np.deg2rad(DE0))*3600; dy=(de-DE0)*3600; r=np.hypot(dx,dy)
z=np.load('/workspace/out/m0_Xb55_spw18.npz'); m0,e0=z['m0'],z['e0']
bowl=np.nanmedian(m0[(r>14)&(r<22)]); e_ring=1.4826*np.nanmedian(np.abs(m0[(r>14)&(r<22)]-bowl))
mask=((m0-bowl)/max(e_ring,np.nanmedian(e0[r<14]))>5)&(r>8)&(r<20)&np.isfinite(m0)
lab,n=ndi.label(mask); sizes=ndi.sum(mask,lab,range(1,n+1)); isl=lab==(1+int(np.argmax(sizes)))
print('disk island pixels',isl.sum(),'centroid dRA,dDec=',dx[isl].mean().round(2),dy[isl].mean().round(2),'r range',r[isl].min().round(1),r[isl].max().round(1))
rad=ra[isl],de[isl]
np.savez('/workspace/out/disk_mask_Xb55.npz',ra=ra[isl],dec=de[isl],dRA=dx[isl],dDec=dy[isl])   # sky positions of the Xb55 pixels (1.2") defining the disk region
rows=[]
def stack_spec(cube,hd,coords,mask_pix=None):
    w=WCS(hd).celestial; x,y=w.world_to_pixel_values(*coords); x=np.round(x).astype(int); y=np.round(y).astype(int)
    ok=(x>=0)&(x<hd['NAXIS1'])&(y>=0)&(y<hd['NAXIS2']); pix=set(zip(y[ok],x[ok]))
    yy,xx=zip(*pix); return np.nanmean(np.array(cube[:,list(yy),list(xx)],float),axis=1),len(pix)
def nullpix(hd,npix,rng,nrand=40):
    w=WCS(hd).celestial; Yn,Xn=np.mgrid[:hd['NAXIS2'],:hd['NAXIS1']]; a,b=w.pixel_to_world_values(Xn,Yn); rr=np.hypot((a-RA0)*np.cos(np.deg2rad(DE0))*3600,(b-DE0)*3600)
    cand=np.argwhere((rr>16)&(rr<26)); out=[]
    for _ in range(nrand):
        c=cand[rng.integers(len(cand))]; sel=cand[(np.abs(cand[:,0]-c[0])<=3)&(np.abs(cand[:,1]-c[1])<=3)]
        out.append(sel[:npix] if len(sel)>=npix else sel)
    return out
rng=np.random.default_rng(3)
LINES=[('Xb55','spw18','HCO+ 2-1',178.375056),('Xb55','spw16','HCN 2-1',177.261111),('Xb47','spw18','HCO+ 4-3',356.734223),('Xb47','spw16','HCN 4-3',354.505478),
       ('Xb4f','spw16','HCO+ 3-2',267.557626),('Xb4f','spw18','HCN 3-2',265.886434)]
for ous,spw,nm,f0 in LINES:
    h=fits.open(f(ous,spw),memmap=True); hd=h[0].header; cube=h[0].data[0]; fr=freqaxis(hd)/1e9; v=(1-fr/f0)*C; dv=abs(np.median(np.diff(v)))
    s,npx=stack_spec(cube,hd,rad)
    def integ(sp):
        rel=v-2780.; bl=(np.abs(rel)>500)&(np.abs(rel)<1500)&np.isfinite(sp); p=np.polyfit(rel[bl],sp[bl],1); sc=sp-np.polyval(p,rel); w=np.abs(rel)<=350
        return np.nansum(sc[w])*dv
    I=integ(s); nulls=[]
    for sel in nullpix(hd,npx,rng):
        if len(sel)<max(3,npx//3): continue
        sp=np.nanmean(np.array(cube[:,sel[:,0],sel[:,1]],float),axis=1)
        if np.isfinite(sp).sum()<100: continue
        nulls.append(integ(sp))
    nulls=np.array(nulls); e=1.4826*np.median(np.abs(nulls-np.median(nulls)))
    rows.append(dict(ous=ous,line=nm,npix=npx,I_Jykms=I,eI_emp=e,snr=(I-np.median(nulls))/e,null_med=np.median(nulls)))
# CH3OH via template fits on the stacked spectra (K km/s)
for ous,spw,nm,fr0s in [('Xb55','spw22','CH3OH Q J=1-3',(165.050175,165.061130,165.099240)),('Xb4f','spw20','CH3OH 2(0)-1(-1) E',(254.015377,))]:
    h=fits.open(f(ous,spw),memmap=True); hd=h[0].header; cube=h[0].data[0]; fr=freqaxis(hd)/1e9; bm=(hd['BMAJ']*3600,hd['BMIN']*3600)
    s,npx=stack_spec(cube,hd,rad); yK=toK(s,fr,bm); A,eA,sg=fit_one(fr,yK,fr0s,v0=2780.,fw=150.)
    k2j=bm[0]*bm[1]*fr0s[0]**2*(1-2780/C)**2/1.222e6
    rows.append(dict(ous=ous,line=nm,npix=npx,I_Jykms=A*k2j,eI_emp=eA*k2j,snr=A/eA,I_Kkms=A,eI_Kkms=eA))
df=pd.DataFrame(rows); df.round(3).to_csv('/workspace/out/disk_stack.csv',index=False); print(df.round(3).to_string())
