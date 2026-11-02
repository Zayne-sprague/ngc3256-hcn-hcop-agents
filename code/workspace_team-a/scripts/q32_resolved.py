# Q32: is the compact CO(2-1) core resolved at 0.18-0.26"? (a) beam-dilution curve fit T(theta)=Te+T0*s^2/(s^2+theta^2)
# from the convolved-beam T_b of Q23 (out/q23/co_hires.json); (b) 2D elliptical Gaussian+const fit to native 0.18" maps
# (peak channel, HCN-FWHM-window mean); (c) cross-check with the independent 0.82x0.68" cube (2015.1.00902.S) native T_b.
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import curve_fit, least_squares
C=299792.458; H=6.62607e-27; KB=1.380649e-16; NU=230.538e9; TBG=2.73
J=lambda T: H*NU/KB/(np.exp(H*NU/KB/T)-1); Jinv=lambda j: H*NU/KB/np.log1p(H*NU/KB/j)
d=json.load(open('out/q23/co_hires.json')); out={}
# (a) curve fit
th=np.array([0.18,0.26,0.38,0.57,0.75,1.0,1.7])
def model(t,Te,T0,s): return Te+T0*s**2/(s**2+t**2)
for c in 'SN':
    b=d['cores'][c]['beams']; T=np.array([b[k]['Tb_RJ'] for k in ['native','0.26','0.38','0.57','0.75','1.0','1.7']])
    th_eff=np.array([np.sqrt(0.1805*0.1739)]+list(th[1:]))
    res={}
    for lab,sel in [('all',slice(None)),('<=1.0',slice(0,6))]:
        # 8% flux-scale error + stat; a 10% missing-flux bias at >=0.75" is treated as sigma
        sig=0.08*T[sel]+0.3
        for Te0 in [None]:
            p,cv=curve_fit(model,th_eff[sel],T[sel],p0=[3,30,0.3],sigma=sig[...] ,bounds=([0,1,0.01],[30,300,3]),absolute_sigma=True)
            e=np.sqrt(np.diag(cv)); res[lab]={'Te':p[0],'T0':p[1],'s_fwhm_arcsec':p[2],'e':e.tolist(),'chi2':float((((model(th_eff[sel],*p)-T[sel])/sig)**2).sum()),'npts':int(np.size(T[sel]))}
    # no-pedestal pure single Gaussian
    p,cv=curve_fit(lambda t,T0,s:model(t,0,T0,s),th_eff,T,p0=[30,0.3],sigma=0.08*T+0.3,bounds=([1,0.01],[300,3]))
    res['nopedestal']={'T0':p[0],'s_fwhm_arcsec':p[1],'e':np.sqrt(np.diag(cv)).tolist(),'chi2':float((((model(th_eff,0,*p)-T)/(0.08*T+0.3))**2).sum())}
    # power-law alternative: T ~ theta^-a at theta>=0.38 (log slope between beams)
    sl=np.polyfit(np.log(th_eff[2:6]),np.log(T[2:6]),1)[0]; res['slope_0.38-1.0']=sl
    res['T']=T.tolist(); out[c]=res
    print(c,'T',np.round(T,1)); print(json.dumps(res,indent=0,default=float)[:900])
# MC for ff(0.26") and (0.18") with fit uncertainty
rng=np.random.default_rng(1)
for c in 'SN':
    b=d['cores'][c]['beams']; T=np.array([b[k]['Tb_RJ'] for k in ['native','0.26','0.38','0.57','0.75','1.0','1.7']])
    th_eff=np.array([np.sqrt(0.1805*0.1739)]+list(th[1:])); ff26=[];ff18=[];T0s=[];ss=[]
    for i in range(3000):
        Tn=T*(1+0.08*rng.standard_normal())*(1+0.04*rng.standard_normal(7))
        try: p,_=curve_fit(model,th_eff,Tn,p0=[3,30,0.3],bounds=([0,1,0.01],[30,300,3]),maxfev=2000)
        except Exception: continue
        s=p[2]; ff26.append(s**2/(s**2+0.26**2)); ff18.append(s**2/(s**2+th_eff[0]**2)); T0s.append(p[1]); ss.append(s)
    q=lambda x:[float(np.percentile(x,k)) for k in (16,50,84)]
    out[c]['MC']={'ff_0.26':q(ff26),'ff_0.18':q(ff18),'T0_RJ':q(T0s),'s_fwhm':q(ss)}
    print(c,'MC',out[c]['MC'])
json.dump(out,open('out/q32/curve_fit.json','w'),indent=1,default=float)
