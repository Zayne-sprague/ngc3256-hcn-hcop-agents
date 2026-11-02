# Q7 injection: flat intrinsic ratio R0 (HCN = R0*HCO+ moment map + real null-window HCN noise, N and S apertures) must return a flat annulus profile (no pipeline-induced gradient).
import sys; sys.path.insert(0,'/workspace/scripts')
from q7_lib import *
sys.argv=['x']; 
exec(open('/workspace/scripts/q7_profiles.py').read().split("rows=[]")[0])   # reuse loaders, null_m0
tb=(0.30,0.30,0.); H,O=cb[('Xb4b','HCN')],cb[('Xb4b','HCO+')]; w=(-350,350)
mO=smooth(O.m0(*w),O,tb); nH=smooth(null_m0(H,w),H,tb)
for R0 in (0.40,0.60,1.00):
  for reg,(x0,y0) in {'N':(0.0,2.55),'S':(-0.15,-2.55)}.items():
    x,y=peak_near(O,mO,x0,y0); out=[]
    for a,b in zip(EDGES['Xb4b@0.30'][:-1],EDGES['Xb4b@0.30'][1:]):
        fH=annulus_flux(H,R0*mO+nH,x,y,a,b,tb)[0]; fO=annulus_flux(O,mO,x,y,a,b,tb)[0]; out.append(fH/fO)
    print(f'R0={R0:.2f} {reg}: recovered annulus R (0-.15,.15-.25,.25-.5,.5-1,1-1.5,1.5-2,2-3") =',np.round(out,3))
