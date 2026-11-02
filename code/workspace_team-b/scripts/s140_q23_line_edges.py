# Q23: channel offset of each ladder line from the nearest edge of its spw, per OUS
import numpy as np
C=299792.458
L={'HCN2-1':177.2611,'HCO+2-1':178.3751,'HCN3-2':265.8864,'HCO+3-2':267.5576,'HCN4-3':354.5055,'HCO+4-3':356.7342,'H13CN4-3':345.3398*0+345.3398}
H13=345.33977  # H13CN 4-3 rest GHz
L['H13CN4-3']=H13
for ous in ['Xb55','Xb47','Xb4f','Xb4d','Xb4b']:
  for s in [16,18,20,22,25,27,29,31]:
    try: z=np.load(f'/workspace/out/regspec/{ous}.spw{s}.cube.pbcor.npz')
    except Exception: continue
    nu=z['nu']/1e9; n=len(nu); d=abs(nu[1]-nu[0])
    for k,f0 in L.items():
      for v in (2780.,):
        f=f0*(1-v/C)
        if nu.min()-0.01<f<nu.max()+0.01:
          i=np.argmin(abs(nu-f)); 
          lo=i if nu[1]<nu[0]==False else None
          # low-frequency edge channel index
          ilo=int(np.argmin(nu)); 
          off=(f-nu.min())/d; offh=(nu.max()-f)/d
          print(f'{ous} spw{s} {k:9s} obs {f:.3f} GHz  ch-from-lowfreq-edge {off:6.0f}  from-high-edge {offh:6.0f} (dnu {d*1e3:.2f} MHz, {n} ch)  nu0={nu[0]:.3f}')
