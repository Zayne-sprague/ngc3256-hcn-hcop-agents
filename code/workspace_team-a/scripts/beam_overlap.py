# N-S separation vs beam: response of each cube beam at the other nucleus
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
sep=np.array([-0.15-0.0,-2.55-2.55]); d=np.hypot(*sep); pa_sep=np.degrees(np.arctan2(sep[0],sep[1]))%180  # PA east of north of N->S line (mod 180)
rows=[]
for ous,spw in [('Xb55','spw18'),('Xb47','spw18'),('Xb4f','spw16'),('Xb4d','spw25'),('Xb4b','spw25')]:
    hd=fits.getheader(f(ous,spw)); bmaj,bmin,bpa=hd['BMAJ']*3600,hd['BMIN']*3600,hd['BPA']
    dpa=np.radians(pa_sep-bpa)
    fw=1/np.sqrt((np.cos(dpa)/bmaj)**2+(np.sin(dpa)/bmin)**2)   # FWHM of beam along N-S separation
    resp=np.exp(-4*np.log(2)*(d/fw)**2)
    rows.append(dict(ous=ous,bmaj=bmaj,bmin=bmin,bpa=bpa,sep_arcsec=d,beam_fwhm_along_sep=fw,response_at_other_nucleus=resp,sep_over_fwhm=d/fw))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/beam_overlap.csv',index=False); print(df.round(4).to_string())
