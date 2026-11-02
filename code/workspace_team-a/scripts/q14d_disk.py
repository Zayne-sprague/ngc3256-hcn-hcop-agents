# Q14d: disk-mask CH3OH limit relative to HCO+ 2-1 (stack from out/disk_stack.csv, Xb55 beam 9.36x6.28, empirical noise): K km/s ratio
import pandas as pd, numpy as np
d=pd.read_csv('/workspace/out/disk_stack.csv').set_index(['ous','line'])
k=1.222e6/(178.375**2*9.4*6.3)   # K per Jy/beam at HCO+ 2-1, native Xb55 beam
H=d.loc[('Xb55','HCO+ 2-1'),'I_Jykms']*k; eH=d.loc[('Xb55','HCO+ 2-1'),'eI_emp']*k
C=d.loc[('Xb55','CH3OH Q J=1-3'),'I_Kkms']; eC=d.loc[('Xb55','CH3OH Q J=1-3'),'eI_Kkms']
r=C/H; e=eC/H
out=pd.DataFrame([dict(region='Disk (Xb55 HCO+2-1 S/N>5, r>8", 136 px)',HCOp21_Kkms=H,eHCOp21=eH,CH3OH_blend_Kkms=C,eCH3OH=eC,ratio=r,eratio=e,lim3=3*e+max(r,0))])
out.to_csv('/workspace/out/q14/disk_limits.csv',index=False); print(out.round(4).to_string())
