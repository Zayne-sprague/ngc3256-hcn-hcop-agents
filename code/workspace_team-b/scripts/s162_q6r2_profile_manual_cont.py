"""Q6 r2 fixes 2-3: (a) model-independent azimuthal-profile half-max size and curve of growth; (b) imaging cross-check: regular Xb4b cube vs
manual 26.6 MHz cube (0.22x0.20", spw25); (c) native (unsmoothed) continuum peaks, T_b, dust Sigma. python3 -I s162_q6r2_profile_manual_cont.py"""
import sys, json, numpy as np, warnings
sys.path.insert(0,'/workspace/scripts'); warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter, label
from regions import REGIONS
C=299792.458; V0=2780.; HW=350.; PC=194.
REG={k:REGIONS[k] for k in ('N','S','Bridge','E1','E2')}
D='/data/member.uid___A001_X1273_Xb4b.ngc_3256_sci.'
def m0map(fn,f0,ra,de,half,pbcorr=True):
    h=fits.open(fn,memmap=True); hd=h[0].header; d=h[0].data
    while d.ndim>3: d=d[0]
    nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']; v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
    pix=abs(hd['CDELT2'])*3600; w=WCS(hd).celestial; cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]; n=int(half/pix); ix,iy=int(round(cx)),int(round(cy))
    win=np.abs(v-V0)<HW; lf=np.abs(v-V0)>600
    cut=np.array(d[:,iy-n:iy+n+1,ix-n:ix+n+1],dtype=np.float64); base=np.nanmedian(cut[lf],axis=0) if lf.sum()>2 else 0; sig=np.nanmedian(np.nanstd(cut[lf],axis=0)) if lf.sum()>2 else np.nan
    # edge channels partly inside the window: weight by overlap
    wgt=np.clip((HW-np.abs(v-V0))/dv+0.5,0,1)
    m0=np.nansum((cut-base)*wgt[:,None,None],axis=0)*dv
    return m0,pix,(hd['BMAJ']*3600,hd['BMIN']*3600),(cx-(ix-n),cy-(iy-n)),sig*np.sqrt(wgt.sum())*dv,int(win.sum()),dv
def halfmax(m,pix,cx,cy,rmax):
    ys,xs=np.indices(m.shape); sm=gaussian_filter(m,0.0376/pix); r=np.hypot(xs-cx,ys-cy)*pix
    k=np.argmax(np.where(r<0.4,sm,-np.inf)); py,px=np.unravel_index(k,m.shape); pk=sm[py,px]
    lab,_=label((sm>0.5*pk)&(np.hypot(xs-px,ys-py)*pix<rmax)); a=(lab==lab[py,px]).sum()*pix**2; return a,pk,(py,px)
def radial(m,pix,py,px,rmax=2.0):
    ys,xs=np.indices(m.shape); r=np.hypot(xs-px,ys-py)*pix; ed=np.arange(0,rmax+1e-6,pix*1.0)
    # first positive-going crossing of half of the central (r<beam/2) mean
    prof=np.array([np.nanmean(m[(r>=a)&(r<b)]) for a,b in zip(ed[:-1],ed[1:])]); rc=0.5*(ed[:-1]+ed[1:]); return rc,prof
def rhalf(rc,prof,pk):
    i=np.where(prof<0.5*pk)[0][0]; return float(np.interp(0.5*pk,[prof[i],prof[i-1]],[rc[i],rc[i-1]]))
def main():
    out={'lines':{}}
    for ln,f0,spw in (('HCO+3-2',267.5576,'spw25'),('HCN3-2',265.8864,'spw27')):
        out['lines'][ln]={}
        for rn,(ra,de,_) in REG.items():
            res={}
            variants=[('regular',D+spw+'.cube.I.pbcor.fits')]
            if spw=='spw25': variants.append(('manual',D+'spw25.cube.I.manual.image.pbcor.fits'))
            for vn,fn in variants:
                m,pix,bm,c,s0,nch,dv=m0map(fn,f0,ra,de,3.0)
                if min(m.shape)<int(5.5/pix) or not np.isfinite(m).all(): print(ln,rn,vn,'outside footprint/partial'); continue
                a,pk,(py,px)=halfmax(m,pix,c[0],c[1],1.2 if rn in('N','Bridge') else 1.5)
                dec=lambda A: np.sqrt(max(4*A/np.pi-bm[0]*bm[1],0))
                rc,prof=radial(m,pix,py,px); pk0=np.nanmean(m[np.hypot(*(np.indices(m.shape)-np.array([py,px])[:,None,None]))*pix<0.5*np.sqrt(bm[0]*bm[1])])  # mean within r<0.5*beam
                # azimuthal half-power: profile normalised to peak of profile's first bin (pixel peak)
                rh=rhalf(rc,prof,prof[0] if prof[0]>0 else pk0); Dobs=2*rh
                # curve of growth of flux (beam-normalised sum within r) and radius with half of the flux inside r<1.5"
                ys,xs=np.indices(m.shape); r=np.hypot(xs-px,ys-py)*pix
                cg=[(float(rr),float(np.nansum(m[r<rr])*pix**2/(1.1331*bm[0]*bm[1]))) for rr in (0.25,0.5,0.75,1.0,1.5,2.0)]
                res[vn]=dict(beam=bm,pix=pix,nchan=nch,dv=dv,peak_m0=float(np.nanmax(m)),halfmax_Ddec=float(dec(a)),radial_halfpower_Dobs=float(Dobs),radial_Ddec=float(np.sqrt(max(Dobs**2-bm[0]*bm[1],0))),cog=cg,sig0=float(s0))
                print(ln,rn,vn,'beam %.3fx%.3f'%bm,'halfmax Ddec %.2f'%dec(a),'radial Dobs %.2f -> Ddec %.2f'%(Dobs,np.sqrt(max(Dobs**2-bm[0]*bm[1],0))),'peak %.3f sig0 %.3f'%(np.nanmax(m),s0),flush=True)
            out['lines'][ln][rn]=res
    # continuum native peaks
    fn=D+'spw25_27_29_31.cont.I.pbcor.fits'; h=fits.open(fn); hd=h[0].header; img=np.squeeze(h[0].data).astype(float)*1e3
    pix=abs(hd['CDELT2'])*3600; bm=(hd['BMAJ']*3600,hd['BMIN']*3600); w=WCS(hd).celestial; yy,xx=np.mgrid[:img.shape[0],:img.shape[1]]
    rr=np.hypot(xx-img.shape[1]/2,yy-img.shape[0]/2)*pix; mk=(rr>8)&(rr<14); rms=1.4826*np.nanmedian(np.abs(img[mk]-np.nanmedian(img[mk])))
    nu=257.76; J=lambda T:(4.799e-11*nu*1e9)/np.expm1(4.799e-11*nu*1e9/T); Jbg=J(2.73); out['cont']={'rms_mJy':float(rms),'beam':bm}
    REGc=dict(REG); REGc['SE']=(156.96684,-43.90635,'SE')
    for n,(ra,de,_) in REGc.items():
        cx,cy=[float(t) for t in w.world_to_pixel_values(ra,de)]; k=int(0.4/pix); ix,iy=int(round(cx)),int(round(cy)); sub=img[iy-k:iy+k+1,ix-k:ix+k+1]
        S=float(np.nanmax(sub)); Tb=1.222e3*S/(nu**2*bm[0]*bm[1]); e={'peak_native_mJy':S,'Tb_RJ':Tb,'Tb_err':1.222e3*rms/(nu**2*bm[0]*bm[1])}
        for Td in (50,100):
            tau=-np.log(1-Tb/(J(Td)-Jbg)); kap=0.1*(nu/1000.)**1.5; e['Sigma_Td%d'%Td]=float(tau/kap*4785.); e['tau_Td%d'%Td]=float(tau)
        out['cont'][n]=e; print('cont',n,'peak %.3f mJy/beam (rms %.3f) Tb %.3f K Sigma50 %.0f Sigma100 %.0f'%(S,rms,Tb,e['Sigma_Td50'],e['Sigma_Td100']))
    json.dump(out,open('/workspace/out/q6r2/s162.json','w'),indent=1)
main()
