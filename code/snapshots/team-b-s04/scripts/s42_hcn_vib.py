# Q4b: HCN v2=1f (l-type doublet 'f') J=3-2 (267.199283 GHz; Eu=1050 K) and J=4-3 (356.256 GHz) at the nuclei: IR-pumping indicator. The 1f line lies at +401 km/s from HCO+ (red wing of HCO+, i.e. 0.36 GHz lower in frequency);
# the 1e partner (265.8527 / 354.4604 GHz, same S*mu^2 within 0.1%) lies +38 km/s (red) from HCN v=0, so alpha=I(1f)/I(HCN v=0) also bounds the 1e contamination of HCN J=3-2/4-3 (relevant to Q1 HCN/HCO+).
# Method as s41: template = HCN v=0 profile in the same aperture (spw with HCN; integral = I_HCN), amplitude alpha (integrated ratio), HCO+ wing = nuisance (exp-linear, exp-quad, power-law about
# HCO+ line centre, constant), profile-likelihood 3-sigma limit (alpha>=0). Noise from line-free channels (>700 km/s from HCO+ and >500 km/s from the vib line), rescaled by sqrt(chi2red).
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
C=299792.458; VSYS=2780.
R='/workspace/out/regspec/'
SETUPS={ # label: (vib file, HCN file, f_vib, f_HCN, f_HCO+)
 'J=3-2 Xb4b':('Xb4b.spw25','Xb4b.spw27',267.199283,265.886434,267.557626),
 'J=3-2 Xb4d':('Xb4d.spw25','Xb4d.spw27',267.199283,265.886434,267.557626),
 'J=3-2 Xb4f':('Xb4f.spw16','Xb4f.spw18',267.199283,265.886434,267.557626),
 'J=4-3 Xb47':('Xb47.spw18','Xb47.spw16',356.256,354.505478,356.734223)}
def vel(nu,f0): return (1-nu/f0)*C-VSYS
def run(label,reg,kind,wlo=-250.,whi=250.):
    fv,fh,Fv,Fh,Fo=SETUPS[label]
    zv=np.load(R+fv+'.cube.pbcor.npz'); zh=np.load(R+fh+'.cube.pbcor.npz')
    if f'{reg}_{kind}' not in zv.files: return []
    nuv=zv['nu']/1e9; nuh=zh['nu']/1e9
    yv=zv[f'{reg}_{kind}']*1e3; yh=zh[f'{reg}_{kind}']*1e3
    vh=vel(nuh,Fh); o=np.argsort(vh); vh,yhs=vh[o],yh[o]
    yhs=yhs-np.median(yhs[(abs(vh)>500)&(abs(vh)<1000)])
    msk=abs(vh)<600; I=np.trapezoid(yhs[msk],vh[msk])       # mJy km/s
    Ti=interp1d(vh,yhs,bounds_error=False,fill_value=0.)
    vv=vel(nuv,Fv); vo=vel(nuv,Fo)
    m=(vv>wlo)&(vv<whi)&np.isfinite(yv)
    if m.sum()<15: return []
    lf=(abs(vo)>600)&(abs(vv)>450)&np.isfinite(yv); sg0=np.std(yv[lf])
    t=Ti(vv[m]); y=yv[m]; u=vo[m]/100.            # t integrates to I for alpha=1
    M={'exp-linear':(lambda p:np.exp(p[0]+p[1]*u),[[0,-.5],[-1,-1],[1,-.2]]),
       'exp-quad':(lambda p:np.exp(p[0]+p[1]*u+p[2]*u**2),[[0,-.5,0],[-1,-1,.05],[1,-.2,-.05]]),
       'power-law @HCO+':(lambda p:p[0]*np.abs(u)**(-p[1]),[[1,1],[.3,2],[3,.5]]),
       'constant':(lambda p:p[0]+0*u,[[0.]])}
    rows=[]
    for nm,(fn,starts) in M.items():
        def prof(al):
            best=None
            for p0 in starts:
                try: r=least_squares(lambda p:(y-al*t-fn(p))/sg0,p0,max_nfev=500)
                except Exception: continue
                if best is None or r.cost<best[0]: best=(r.cost,r.x)
            return 2*best[0] if best else np.inf
        als=np.linspace(0,0.5,101); ch=np.array([prof(a) for a in als]); i=ch.argmin(); dof=m.sum()-len(starts[0])-1; s2=max(ch[i]/dof,1.)
        d=(ch-ch[i])/s2; ok=np.where((d>=9)&(als>=als[i]))[0]; ul=als[ok[0]] if len(ok) else np.nan
        sel=als[d<1] if np.any(d<1) else np.array([np.nan])
        rows.append(dict(setup=label,region=reg,kind=kind,model=nm,alpha_best=als[i],a1_lo=np.nanmin(sel),a1_hi=np.nanmax(sel),alpha_3s_ul=ul,chi2red=ch[i]/dof,n=int(m.sum()),
                         sig_mJy=sg0*np.sqrt(s2),I_HCN_mJykms=I,vmin_covered=float(vv[m].min()),vmax=float(vv[m].max())))
    return rows
if __name__=='__main__':
    rows=[]
    for lab in SETUPS:
        for reg in ('N','S','Bridge'):
            for kind in ('pix','ap'):
                rows+=run(lab,reg,kind)
    df=pd.DataFrame(rows); df.to_csv('/workspace/out/s42_hcn_vib.csv',index=False)
    pd.set_option('display.width',250,'display.max_rows',500); print(df.round(4).to_string())
    ok=df[df.chi2red<1.5]
    print('\nACCEPTED chi2red<1.5; range of 3-sigma alpha(1f/HCN) limit'); print(ok.groupby(['setup','region','kind']).alpha_3s_ul.agg(['min','max','count']).round(3).to_string())
