import numpy as np
def J(nu,T):
    x=4.799e-11*nu*1e9; return x/np.expm1(x/T)
for nu,nm,vals in [(265.886,'HCN',{'N':(2.64,4.7),'S':(1.34,2.6)}),(267.558,'HCO+',{'N':(2.82,5.1),'S':(2.94,5.7)})]:
  for Tk in (20,25,30,50,80,150):
    d=J(nu,Tk)-J(nu,2.73)
    print(nm,Tk,{k:tuple(round(v/d,3) for v in vv) for k,vv in vals.items()})
# 10" core share check: S HCO+ 3-2 10"=9.95 K km/s, core Tint*dil*1.0645*FWHMv
for nm,T,dil,fw,F in [('S HCO+ gauss',5.7,0.52,117,9.95),('N HCN gauss',4.7,0.56,150,4.81)]:
  print(nm,T*dil*1.0645*fw/F)
