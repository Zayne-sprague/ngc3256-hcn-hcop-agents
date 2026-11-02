# Q1 fix 2: is the disk 3-2 deficit (disk 3-2/2-1 Tb ratio below the nuclear one) a missing-short-spacing effect?  Image products only -> parametric uv-hole model.
# Assumptions (stated): both the Band 5 (2-1) and Band 6 Xb4f (3-2) maps are compact-array maps whose longest baseline is ~ lambda/theta_maj (-> B_max ~ 40 m, i.e. 7-m-array-like),
# shortest baseline B_min = 8.7 m (7-m dish shadowing limit) or 15 m (12-m shadowing limit) [ALMA Technical Handbook values, external assumption].
# A smooth Gaussian disk of FWHM th retains a fraction exp[-(pi th B_min/lambda)^2/(4 ln2)] of its central brightness after a uv hole of radius B_min (exact for a Gaussian).
# Observed: disk/compact 3-2/2-1 deficit d_obs (Tb ratios).  Model: d(th)= retained_3-2(th)/retained_2-1(th).  Solve th; compare to the mask extent.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.optimize import brentq
C=299792.458; ARC=206265.
nu={'2-1 (Band 5)':177.26/(1+2780/C),'2-1 HCO+':178.375/(1+2780/C),'3-2 (Band 6 Xb4f)':267.56/(1+2780/C)}
lam=lambda g:2.998e8/(g*1e9)       # m
def ret(th,B,g): x=(np.pi*th/ARC*B/lam(g))**2/(4*np.log(2)); return np.exp(-x)
rows=[]
for B in (8.7,15.0):
    for th in (6,8,10,12,14,16,20,25,30,40):
        r6=ret(th,B,nu['3-2 (Band 6 Xb4f)']); r5=ret(th,B,nu['2-1 HCO+'])
        rows.append(dict(B_min_m=B,disk_FWHM_arcsec=th,retained_3_2=r6,retained_2_1=r5,deficit_32_over_21=r6/r5,MRS_6_arcsec=0.983*lam(nu['3-2 (Band 6 Xb4f)'])/B*ARC,MRS_5_arcsec=0.983*lam(nu['2-1 HCO+'])/B*ARC))
M=pd.DataFrame(rows); M.to_csv('/workspace/out/q1/q1_disk_recovery_model.csv',index=False); print(M.round(3).to_string())
# observed deficit: disk-mask 3-2(Xb4f)/2-1 Tb ratio over the 'compact' 3-2(Xb4f)/2-1 ratio (deblended N,S and 10" matched N,S,Bridge)
D=pd.read_csv('/workspace/out/q1/q1_disk_ratios.csv',keep_default_na=False); L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv')
comp={'HCO+':[0.428,0.401,0.497,0.478],'HCN':[0.350,0.333,0.383,0.365]}      # deblend N,S ([Xb4f 3-2]/2-1) ; 10" matched N,S point values
out=[]
for sp,nm in (('HCO+','HCO+ 3-2/2-1'),('HCN','HCN 3-2/2-1')):
    ref=np.mean(comp[sp]); eref=np.std(comp[sp])
    for zl in ('none','ring','null'):
        q=D[(D.zero_level==zl)&(D.ratio==nm)].iloc[0]; d=q.R/ref; 
        row=dict(species=sp,zero_level=zl,disk_ratio=q.R,e_stat=q.eR_stat,compact_ratio=ref,compact_spread=eref,deficit_obs=d)
        for B in (8.7,15.0):
            f=lambda th: ret(th,B,nu['3-2 (Band 6 Xb4f)'])/ret(th,B,nu['2-1 HCO+'])-d
            try: row[f'theta_FWHM_for_deficit_B{B}']=brentq(f,2,200)
            except Exception: row[f'theta_FWHM_for_deficit_B{B}']=np.nan
        out.append(row)
O=pd.DataFrame(out); O.to_csv('/workspace/out/q1/q1_disk_deficit_interpretation.csv',index=False); print(O.round(3).to_string())
# recovery-independent diagnostics: (HCN/HCO+ 3-2)/(HCN/HCO+ 2-1) double ratio (disk) vs the same at the nuclei (10" matched) -> equals 1 if HCN and HCO+ share the spatial distribution and recovery
dd=[]
for zl in ('none','ring','null'):
    a=D[(D.zero_level==zl)&(D.ratio=='HCN/HCO+ 3-2')].iloc[0]; b=D[(D.zero_level==zl)&(D.ratio=='HCN/HCO+ 2-1')].iloc[0]
    dd.append(dict(zero_level=zl,disk_HCN_HCOp_3_2=a.R,disk_HCN_HCOp_2_1=b.R,double_ratio_disk=a.R/b.R,double_ratio_N=0.478/0.620,double_ratio_S=0.429/0.563,e_stat=a.R/b.R*np.hypot(a.eR_stat/a.R,b.eR_stat/b.R)))
DD=pd.DataFrame(dd); DD.to_csv('/workspace/out/q1/q1_disk_double_ratio.csv',index=False); print(DD.round(3).to_string())
