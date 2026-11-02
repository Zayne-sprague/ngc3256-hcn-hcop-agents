import numpy as np
from scipy.interpolate import RegularGridInterpolator as R
d=np.load('/workspace/out/q3/grid_lvg.npz'); T0=np.log10(d['TK']);N0=d['LN'];C0=d['LNC']
print(d['hcn'][0,0,0])
for key,(T,n,N,ref) in {'a':(150,10**4.62,10**13.5,[.03602,.01178,.004214]),'b':(150,10**4.62,10**14.3,[.2293,.07581,.0266]),'c':(60,10**4.3,1e14,[.03304,.007259,.001878]),'d':(35,10**4.8,10**13.8,[.04231,.008499,.001522])}.items():
    out=[]
    for j in range(3):
        f=R((T0,N0,C0),np.log(np.clip(d['hcn'][...,j,0],1e-30,None)))
        out.append(np.exp(f([[np.log10(T),np.log10(n),np.log10(N)]])[0]))
    print(key,np.round(np.array(out)/ref,2))
