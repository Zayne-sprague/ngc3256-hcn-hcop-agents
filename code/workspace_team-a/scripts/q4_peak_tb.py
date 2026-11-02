# Q4h: peak brightness temperature (per channel, native beam, Rayleigh-Jeans) of HCN/HCO+ at the region peak pixels: filling-factor / opacity context. From out/line_table_v2.csv (peakpix). Output out/q4/peak_tb.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
LT=pd.read_csv('/workspace/out/line_table_v2.csv'); LT=LT[(LT.kind=='peakpix')&LT.line.isin(['HCN 2-1','HCO+ 2-1','HCN 3-2','HCO+ 3-2','HCN 4-3','HCO+ 4-3'])&LT.region.isin(['N','S','Bridge'])]
NU={'HCN 2-1':177.261111,'HCO+ 2-1':178.375056,'HCN 3-2':265.886434,'HCO+ 3-2':267.557626,'HCN 4-3':354.505478,'HCO+ 4-3':356.734223}
rows=[]
for r in LT.itertuples():
    k=1.222e6/(NU[r.line]**2*r.bmaj*r.bmin)*1e-3      # K per mJy/beam (rest frequency; sky/rest differ by 0.9%)
    rows.append(dict(ous=r.ous,region=r.region,line=r.line,beam=f'{r.bmaj:.2f}x{r.bmin:.2f}',peak_mJy_beam=r.peak_mJy,rms_mJy=r.rms_mJy,Tb_peak_K=r.peak_mJy*k,eTb_K=r.rms_mJy*k,dv=r.dv))
R=pd.DataFrame(rows); R.round(3).to_csv('/workspace/out/q4/peak_tb.csv',index=False); pd.set_option('display.width',200); print(R.round(2).to_string())
