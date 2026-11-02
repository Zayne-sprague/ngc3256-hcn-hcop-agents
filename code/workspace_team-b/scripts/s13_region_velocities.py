# Per-region systemic (line-centre) velocities from HCO+/HCN 3-2 (1" and 0.25" beams): flux-weighted centroid of >3 sigma channels, half-maximum midpoint, and
# 1G centroid with its chi2r (flag invalid where chi2r>2). Bootstrap-free: error = scatter between the 4 line/array estimates and formal rms term.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import pandas as pd
fits_=pd.read_csv('/workspace/out/s7_line_fits.csv')
SETS=[('Xb4b','spw25',267.55763,'HCO+ 3-2'),('Xb4b','spw27',265.88618,'HCN 3-2'),('Xb4d','spw25',267.55763,'HCO+ 3-2'),('Xb4d','spw27',265.88618,'HCN 3-2')]
rows=[]
for ous,spw,f0,lab in SETS:
    z=np.load(f'/workspace/out/regspec/{ous}.{spw}.cube.pbcor.npz'); v=vrad(z['nu'],f0*1e9); o=np.argsort(v)
    for reg in REGIONS:
        y=z[reg+'_ap'][o]; vv=v[o]; ok=np.isfinite(y)
        fl=ok&(abs(vv-VSYS)>450)&(abs(vv-VSYS)<900); 
        if fl.sum()<20: continue
        yb=y[fl]; base=np.median(yb); rms=np.std(yb); yy=y-base
        sm=np.convolve(yy,np.ones(3)/3,'same'); win=abs(vv-VSYS)<300
        ip=np.argmax(np.where(win,sm,-1e9)); pk=sm[ip]
        if pk<6*rms/np.sqrt(3): continue
        # half-max crossings (outermost within +-300 of peak)
        hm=pk/2; sel=np.where(win&(sm>=hm))[0]; lo,hi=sel.min(),sel.max()
        def cross(i0,i1):
            return vv[i0]+(hm-sm[i0])*(vv[i1]-vv[i0])/(sm[i1]-sm[i0]) if sm[i1]!=sm[i0] else vv[i0]
        vlo=cross(lo-1,lo) if lo>0 else vv[lo]; vhi=cross(hi,hi+1) if hi<len(vv)-1 else vv[hi]
        cl=win&(yy>3*rms); cen=np.sum(yy[cl]*vv[cl])/np.sum(yy[cl]); 
        cen_err=rms*np.sqrt(np.sum(((vv[cl]-cen))**2))/np.sum(yy[cl])
        f=fits_[(fits_.ous==ous)&(fits_.region==reg)&(fits_.spec=='ap')&(fits_.line==lab)]
        g1=f.g1_v0.values[0] if len(f) else np.nan; chi=f.g1_chi2r.values[0] if len(f) else np.nan
        rows.append(dict(ous=ous,line=lab,region=reg,v_centroid=cen,v_centroid_err=cen_err,v_halfmax_mid=0.5*(vlo+vhi),fwhm_halfmax=vhi-vlo,v_peak=vv[ip],v_g1=g1,g1_chi2r=chi,
                         g1_valid=bool(chi<2) if np.isfinite(chi) else False,peak_mJy=pk*1e3,snr_peak=pk/rms))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s13_region_velocities.csv',index=False)
pd.set_option('display.width',220); print(df.round(1).to_string())
g=df.groupby('region').agg(v_mid_mean=('v_halfmax_mid','mean'),v_mid_std=('v_halfmax_mid','std'),v_cen_mean=('v_centroid','mean'),v_cen_std=('v_centroid','std'),n=('v_centroid','size'),fwhm_mean=('fwhm_halfmax','mean'))
g['dv_vs_2780']=g.v_mid_mean-VSYS; g.to_csv('/workspace/out/s13_region_vsys.csv'); print(g.round(1).to_string())
