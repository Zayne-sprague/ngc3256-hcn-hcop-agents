# Q4c: line profiles / opacity & absorption diagnostics at 0.25" (Xb4b): HCN 3-2 (spw27) and HCO+ 3-2 (spw25), native beam, pbcor, r<0.25..2" annuli about the N and S continuum peaks.
# (i) Gaussian FWHM/centroid/peak of each line per annulus (opacity broadening / self-absorption dips would make HCN or HCO+ broader/flatter/asymmetric);
# (ii) velocity-resolved ratio R(v) in 100-km/s bins vs constant R (chi2) = test for dips in either line; (iii) peak brightness temperatures (Rayleigh-Jeans, native beam) and line/continuum
# brightness (continuum absorption is only possible if T_cont > T_ex; here T_cont,b << 3 K); (iv) negative-channel test inside the line window (a dip below the contsub zero would require absorption).
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import json
from scipy.optimize import curve_fit
FHCN,FHCO=265.886434,267.557626
def load(spw):
    fn=D+f'member.uid___A001_X1273_Xb4b.ngc_3256_sci.{spw}.cube.I.pbcor.fits'
    h=fits.open(fn,memmap=True); hd=h[0].header; nu=hd['CRVAL3']+(np.arange(hd['NAXIS3'])+1-hd['CRPIX3'])*hd['CDELT3']
    return h,hd,nu
def region_spectra(h,hd,ra,dec,radii,half=70):
    wc=WCS(hd).celestial; x,y=wc.world_to_pixel_values(ra,dec); ps=abs(hd['CDELT1'])*3600
    x0,y0=int(round(float(x))),int(round(float(y)))
    cube=np.array(h[0].data[0,:,y0-half:y0+half+1,x0-half:x0+half+1],dtype=np.float32)
    yy,xx=np.mgrid[-half:half+1,-half:half+1]; r=np.hypot(xx-(float(x)-x0),yy-(float(y)-y0))*ps
    out={}
    for lo,hi in radii:
        m=(r>=lo)&(r<hi) if lo>0 else (r<hi); out[(lo,hi)]=(np.nanmean(cube[:,m],axis=1),int(m.sum()))
    return out
gauss=lambda v,A,v0,s,c: A*np.exp(-0.5*((v-v0)/s)**2)+c
def fitg(v,y,sig):
    try:
        p,cv=curve_fit(gauss,v,y,p0=[y.max(),0,60,0],sigma=np.full_like(y,sig),absolute_sigma=True,maxfev=5000)
        r=(y-gauss(v,*p))/sig; return p,np.sqrt(np.diag(cv)),np.sum(r**2)/(len(v)-4)
    except Exception: return None,None,None
if __name__=='__main__':
    radii=[(0,.25),(.25,.5),(.5,1.),(1.,2.)]
    hA,hdA,nuA=load('spw27'); hB,hdB,nuB=load('spw25')
    bm=np.sqrt(hdA['BMAJ']*hdA['BMIN'])*3600; print('beam HCN',hdA['BMAJ']*3600,hdA['BMIN']*3600,'HCO+',hdB['BMAJ']*3600,hdB['BMIN']*3600)
    cont=load2d(D+'member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')[0]
    res={}; lines=[]
    for reg in ('N','S'):
        ra,dec,_=REGIONS[reg]
        sA=region_spectra(hA,hdA,ra,dec,radii); sB=region_spectra(hB,hdB,ra,dec,radii)
        vA=C*(1-nuA/(FHCN*1e9)*1.0)-VSYS; vB=C*(1-nuB/(FHCO*1e9))-VSYS
        for rr in radii:
            yA,n=sA[rr]; yB,_=sB[rr]
            fl=(abs(vA)>500)&(abs(vA)<1000); sgA=np.std(yA[fl]); flB=(abs(vB)>500)&(abs(vB)<1000); sgB=np.std(yB[flB])
            w=abs(vA)<500; wB=abs(vB)<500
            pA,eA,cA=fitg(vA[w],yA[w],sgA); pB,eB,cB=fitg(vB[wB],yB[wB],sgB)
            dv=abs(np.median(np.diff(vA)))
            IA=np.sum(yA[abs(vA)<350])*dv; IB=np.sum(yB[abs(vB)<350])*dv
            # velocity-resolved ratio in 100 km/s bins (relative to each line's rest at VSYS)
            edges=np.arange(-300,301,100); Rv=[];Re=[]
            for a,b in zip(edges[:-1],edges[1:]):
                ma=(vA>=a)&(vA<b); mb=(vB>=a)&(vB<b)
                fa,fb=yA[ma].sum()*dv,yB[mb].sum()*dv; ea=sgA*np.sqrt(ma.sum())*dv; eb=sgB*np.sqrt(mb.sum())*dv
                Rv.append(fa/fb if fb>0 else np.nan); Re.append(abs(fa/fb)*np.hypot(ea/max(fa,1e-9),eb/max(fb,1e-9)) if fb>0 and fa>0 else np.nan)
            Rv=np.array(Rv);Re=np.array(Re); Rm=IA/IB; g=np.isfinite(Rv)&np.isfinite(Re)&(Re>0)
            chi=np.sum(((Rv[g]-Rm)/Re[g])**2); dof=g.sum()-1
            # peak Tb [K] RJ: T=1.222e3*S[mJy/beam]/(nu^2 th1 th2)
            def tb(S,nu,hd): return 1.222e3*S*1e3/(nu**2*hd['BMAJ']*3600*hd['BMIN']*3600)
            TA=tb(gauss(*pA[:1],0,0,0)[()] if False else pA[0]+pA[3],FHCN*(1-VSYS/C),hdA); TB=tb(pB[0]+pB[3],FHCO*(1-VSYS/C),hdB)
            # continuum brightness at the same position (pixel nearest)
            nlow=int(np.sum((yA[abs(vA)<400]<-3*sgA)))+int(np.sum((yB[abs(vB)<400]<-3*sgB)))
            row=dict(region=reg,r_in=rr[0],r_out=rr[1],npix=n,sig_HCN_mJy=sgA*1e3,sig_HCO_mJy=sgB*1e3,I_HCN=IA,I_HCO=IB,R=Rm,
                     HCN_peak_mJy=1e3*pA[0],HCN_v0=pA[1],HCN_FWHM=2.355*pA[2],eFWHM_HCN=2.355*eA[2],HCO_peak_mJy=1e3*pB[0],HCO_v0=pB[1],HCO_FWHM=2.355*pB[2],eFWHM_HCO=2.355*eB[2],
                     FWHMratio=pA[2]/pB[2],chi2red_gHCN=cA,chi2red_gHCO=cB,Tpk_HCN_K=TA,Tpk_HCO_K=TB,Rv=list(np.round(Rv,3)),Rv_err=list(np.round(Re,3)),Rv_chi2=chi,Rv_dof=int(dof),n_neg3sig_channels=nlow)
            lines.append(row)
    import pandas as pd
    df=pd.DataFrame(lines); df.to_csv('/workspace/out/s43_profiles.csv',index=False)
    pd.set_option('display.width',280,'display.max_columns',50); print(df.drop(columns=['Rv','Rv_err']).round(3).to_string())
    for r in lines: print(r['region'],r['r_in'],r['r_out'],'R(v)=',r['Rv'],'+-',r['Rv_err'])
    # continuum at the peak pixels (0.24" beam) and Tb_cont
    wc=WCS(fits.getheader(D+'member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')).celestial
    hc=fits.getheader(D+'member.uid___A001_X1273_Xb4b.ngc_3256_sci.spw25_27_29_31.cont.I.pbcor.fits')
    for reg in ('N','S'):
        ra,dec,_=REGIONS[reg]; x,y=wc.world_to_pixel_values(ra,dec); x,y=int(round(float(x))),int(round(float(y)))
        pk=np.nanmax(cont[y-6:y+7,x-6:x+7]); nu=257.76
        T=1.222e3*pk*1e3/(nu**2*hc['BMAJ']*3600*hc['BMIN']*3600); print(reg,'cont peak',pk*1e3,'mJy/beam  Tb_cont =',round(T,3),'K (beam',hc['BMAJ']*3600,hc['BMIN']*3600,')')
