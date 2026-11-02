# Q13f: prediction of the lobe HCN 4-3 blue-window flux (Xb47 at 4.5") if the 3-2 excess is real HCN vs if it is an HDO blend (no 4-3 counterpart). Reads out/q13/q13c_table.csv
import pandas as pd, numpy as np
T=pd.read_csv('/workspace/out/q13/q13c_table.csv'); a=T[(T.pos=='Soutflow')&(T.J=='3-2')].iloc[0]; b=T[(T.pos=='Soutflow')&(T.J=='4-3')].iloc[0]
nu=(354.505477/265.886434)**2     # S43/S32 for equal Tb
r43=(b.HCN_core/a.HCN_core)/nu; er43=r43*np.hypot(a.sHCN_core/a.HCN_core,b.sHCN_core/b.HCN_core)
S32=a.HCN_blue; eS=a.sHCN_blue; E32=a.E_blue
real=S32*nu*r43; ereal=real*np.hypot(eS/S32,er43/r43)
hdo=(S32-E32)*nu*r43; ehdo=abs(S32-E32)*nu*np.hypot(eS/abs(S32-E32)*r43,er43) if S32!=E32 else np.nan
obs=b.HCN_blue; eo=b.sHCN_blue
print(f'4.5" beam, lobe aperture r=1": HCN 3-2 blue {S32:.3f}+-{eS:.3f}; E32 (k_core={a.k_core:.3f}) {E32:.3f}; core Tb r43/32 = {r43:.2f}+-{er43:.2f}; (nu43/nu32)^2={nu:.3f}')
print(f'predicted HCN 4-3 blue-window (Jy/beam km/s): real HCN wing {real:.3f}+-{ereal:.3f}; HDO-only excess (wing = k*HCO+ only) {hdo:.3f}+-{ehdo:.3f}; observed {obs:.3f}+-{eo:.3f}')
print(f'observed - real = {(obs-real)/np.hypot(eo,ereal):+.2f} sigma ; observed - HDO-only = {(obs-hdo)/np.hypot(eo,ehdo):+.2f} sigma; Tb r43/32 of the wing = {obs/S32/nu:.2f}+-{obs/S32/nu*np.hypot(eo/obs,eS/S32):.2f} vs core {r43:.2f}')
