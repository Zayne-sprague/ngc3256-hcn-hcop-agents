# Free Gaussian fit (amplitude, centroid, FWHM, linear baseline) of the 251.66-251.68 GHz feature in all regions where it is seen; centroid -> velocity for CH3OH 2(0)-1(-1)E (254.015377) and NS 11/2-9/2 f (253.968393/253.970633) vs the HCO+ 3-2 velocity of the same region/OUS.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
from ch3oh import d,toK,beam,C
from ch3oh_lines import ref
def g(f,a,fc,fw,b0,b1): return a*np.exp(-4*np.log(2)*((f-fc)/fw)**2)+b0+b1*(f-251.67)
rows=[]
for ous,spw in (('Xb4f','spw20'),('Xb4d','spw29'),('Xb4b','spw29')):
    fr=d[f'{ous}__{spw}__freq']
    for reg in ['N','S','Bridge','Soutflow','Sext','SWarm','NEclump','Eclump']:
        k=f'{ous}__{spw}__{reg}__pk'
        if k not in d.files: continue
        s=d[k]; m=(fr>251.40)&(fr<251.95)&np.isfinite(s)
        v0,fw0=ref(ous,reg)
        try:
            p,cv=curve_fit(g,fr[m],s[m]*1e3,p0=[5,251.67,0.1,0,0],sigma=None,bounds=([0,251.60,0.03,-20,-100],[100,251.74,0.35,20,100]),maxfev=20000)
        except Exception as e: continue
        e=np.sqrt(np.diag(cv)); res=s[m]*1e3-g(fr[m],*p); sg=np.std(res); e=e*np.sqrt(1.0)  # curve_fit scales by residuals
        fc=p[1]; vC=(1-fc/254.015377)*C; vN=(1-fc/253.968393)*C; dv=e[1]/fc*C
        I=p[0]*p[2]*1.0645/fc*C*1e-3  # Jy/beam km/s
        rows.append(dict(ous=ous,region=reg,peak_mJy=p[0],snr_peak=p[0]/sg,fc_GHz=fc,e_fc_MHz=e[1]*1e3,fwhm_kms=p[2]/fc*C,I_Jykms=I,v_if_CH3OH=vC,v_if_NSf=vN,e_v=dv,v_HCO=v0,dv_CH3OH=vC-v0,dv_NS=vN-v0))
df=pd.DataFrame(rows); df.round(3).to_csv('/workspace/out/line251_centroids.csv',index=False)
pd.set_option('display.width',250); print(df.round(2).to_string())
