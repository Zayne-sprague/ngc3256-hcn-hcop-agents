# Q8-2: T43/T32 (Tb) at N/S/Bridge from out/q8/matched_m0.csv for each matched beam and each 3-2 array; N/S contrast; 3-2 array-to-array recovery.
# Statistical errors only from null-window rms; the 10% flux scale per OUS cancels in the N/S contrast (same maps) but not in absolute ratios.
import pandas as pd, numpy as np
R=pd.read_csv('/workspace/out/q8/matched_m0.csv'); rows=[]; rec=[]
def get(sp,tr,beam,reg,col): 
    r=R[(R.line==sp)&(R.transition==tr)&(R.beam==beam)&(R.region==reg)].iloc[0]; return r[col],r['eTb']*(1 if col.startswith('Tb') else 1)
for sp in ('HCN','HCO+'):
  for beam in R.beam.unique():
    for a32 in ('3-2 Xb4b','3-2 Xb4d','3-2 Xb4f'):
      if R[(R.line==sp)&(R.transition==a32)&(R.beam==beam)].empty: continue
      for col in ('Tb_pos_Kkms','Tb_peak_Kkms'):
        d={}
        for reg in ('N','S','Bridge'):
          t43,e43=get(sp,'4-3 Xb47',beam,reg,col); t32,e32=get(sp,a32,beam,reg,col)
          r=t43/t32; e=r*np.hypot(e43/t43,e32/t32); d[reg]=(r,e)
          rows.append(dict(line=sp,beam=beam,arr32=a32,measure=col[3:7],region=reg,T43_T32=r,err=e))
        for a,b in (('N','S'),('N','Bridge'),('S','Bridge')):
          c=d[a][0]/d[b][0]; ec=c*np.hypot(d[a][1]/d[a][0],d[b][1]/d[b][0])
          rows.append(dict(line=sp,beam=beam,arr32=a32,measure=col[3:7],region=f'{a}/{b}',T43_T32=c,err=ec))
    for reg in ('N','S','Bridge'):
      for col in ('Tb_pos_Kkms',):
        f=get(sp,'3-2 Xb4f',beam,reg,col) if not R[(R.line==sp)&(R.transition=='3-2 Xb4f')&(R.beam==beam)].empty else None
        if f is None: continue
        for a in ('3-2 Xb4b','3-2 Xb4d'):
          x=get(sp,a,beam,reg,col); rec.append(dict(line=sp,beam=beam,region=reg,array=a,ratio_to_Xb4f=x[0]/f[0],err=x[0]/f[0]*np.hypot(x[1]/x[0],f[1]/f[0])))
C=pd.DataFrame(rows); C.to_csv('/workspace/out/q8/T43_T32_contrast.csv',index=False)
Q=pd.DataFrame(rec); Q.to_csv('/workspace/out/q8/recovery_32_arrays.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500)
print(C[(C.measure=='pos_')].round(3).to_string()); print(Q.round(3).to_string())
