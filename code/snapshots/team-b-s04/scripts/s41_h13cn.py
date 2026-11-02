# Q4a: H13CN 4-3 (rest 345.339760 GHz; sky ~342.14 GHz) in Xb47 spw20, which sits on the blue wing of CO 3-2 (rest 345.79599; band edge 342.396 GHz sky).
# Template = HCN 4-3 profile of the same region/aperture (Xb47 spw16, same velocity frame), amplitude alpha=I(H13CN)/I(HCN) (integrated); CO wing = smooth nuisance model, fitted
# simultaneously; profile-likelihood one-sided 3 sigma limit on alpha>=0 for every wing model and window; noise per channel from line-free 340.7-341.6 GHz sky, rescaled by sqrt(chi2_red) of each fit.
# Also: HC15N 4-3 (344.200) and the template-width variant (0.7x FWHM). Spectra: beam-sized aperture and peak pixel (out/regspec from s7a). Units mJy/beam, pbcor.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, json
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
C=299792.458; F13=345.339760; FCO=345.795990; FHCN=354.505478; VSYS=2780.
z20=np.load('/workspace/out/regspec/Xb47.spw20.cube.pbcor.npz'); z16=np.load('/workspace/out/regspec/Xb47.spw16.cube.pbcor.npz')
nu20=z20['nu']/1e9; nu16=z16['nu']/1e9
def template(reg,kind,wscale=1.0):
    T=z16[f'{reg}_{kind}']*1e3; v=(1-nu16/FHCN)*C-VSYS; o=np.argsort(v); v,T=v[o],T[o]
    base=np.median(T[(abs(v)>500)&(abs(v)<900)]); T=T-base
    I=np.trapezoid(T[abs(v)<600],v[abs(v)<600])      # mJy km/s (beam-aperture HCN 4-3)
    f=interp1d(v/wscale,T,bounds_error=False,fill_value=0.)
    return f,I
def fit_region(reg,kind,fmin,fmax,wscale=1.0,verbose=False):
    y_all=z20[f'{reg}_{kind}']*1e3; lf=(nu20>340.7)&(nu20<341.6); sg0=np.std(y_all[lf])
    Ti,I=template(reg,kind,wscale); m=(nu20>fmin)&(nu20<fmax); f=nu20[m]; y=y_all[m]; vh=(1-f/F13)*C-VSYS; t=Ti(vh)*wscale  # same integrated flux per unit alpha: scale amplitude by wscale^-1... (area conserved)
    t=t/wscale*wscale   # template shape compressed in v keeps peak; area = I*wscale  -> normalise so alpha is still integrated ratio
    area=np.trapezoid(Ti(np.linspace(-600,600,2401)*1.0),np.linspace(-600,600,2401)); t=Ti(vh)*(I/area)
    x=(f-FCO*(1-VSYS/C))/0.1; u=(f-342.0)/0.1
    M={'exp-quad (Gaussian wing)':(lambda p:np.exp(p[0]+p[1]*u+p[2]*u**2),[[2,.5,0],[0,1,.05],[3,.5,.02],[1,.3,-.02]]),
       'exp-cubic':(lambda p:np.exp(p[0]+p[1]*u+p[2]*u**2+p[3]*u**3),[[2,.5,0,0],[0,1,.05,0],[3,.5,.02,.01]]),
       'exp-linear':(lambda p:np.exp(p[0]+p[1]*u),[[2,.5],[0,1],[3,.8]]),
       'power-law @CO':(lambda p:p[0]*np.abs(x)**(-p[1]),[[200,2],[100,3],[500,1.5]]),
       'CO Gaussian (fixed centre)+offset':(lambda p:p[0]*np.exp(-0.5*((f-FCO*(1-VSYS/C))/(p[1]/C*FCO))**2)+p[2],[[150,150,0],[300,200,0],[100,100,0]])}
    rows=[]
    for nm,(fn,starts) in M.items():
        def prof(al):
            best=None
            for p0 in starts:
                try: r=least_squares(lambda p:(y-al*t-fn(p))/sg0,p0,max_nfev=600)
                except Exception: continue
                if best is None or r.cost<best[0]: best=(r.cost,r.x)
            return (2*best[0],best[1]) if best else (np.inf,None)
        als=np.linspace(0,0.6,121); ch=np.array([prof(a)[0] for a in als]); i=ch.argmin()
        npar=len(starts[0]); dof=len(y)-npar-1; s2=max(ch[i]/dof,1.0)
        d=(ch-ch[i])/s2
        ul=als[np.where((d>=9)&(als>=als[i]))[0][0]] if np.any((d>=9)&(als>=als[i])) else np.nan
        sel=als[d<1]; rows.append(dict(region=reg,kind=kind,wscale=wscale,window=f'{fmin:.2f}-{fmax:.2f}',model=nm,alpha_best=als[i],a1_lo=sel.min(),a1_hi=sel.max(),alpha_3s_ul=ul,
                    chi2red=ch[i]/dof,ndata=int(m.sum()),HCN_H13CN_ll=(1/ul if ul and ul>0 else np.nan),I_HCN_mJykms=I,sig_chan_mJy=sg0*np.sqrt(s2)))
    return rows
if __name__=='__main__':
    rows=[]
    for reg in ('N','S','Bridge'):
        for kind in ('ap','pix'):
            for fmin,fmax in ((341.80,342.20),(341.85,342.12),(341.80,342.28)):
                for w in (1.0,0.7):
                    if w!=1.0 and (kind=='pix' or (fmin,fmax)!=(341.80,342.20)): continue
                    rows+=fit_region(reg,kind,fmin,fmax,w)
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/s41_h13cn.csv',index=False)
    pd.set_option('display.width',250,'display.max_rows',500); print(df.round(3).drop(columns=['ndata']).to_string())
    # acceptable models: chi2red<1.3
    ok=df[df.chi2red<1.3]
    print('\nACCEPTED (chi2red<1.3), 3-sigma alpha limit range per region/kind/window:')
    print(ok.groupby(['region','kind','window','wscale']).alpha_3s_ul.agg(['min','max','count']).round(3).to_string())
