# Q38 step 2: Gaussian dust disk (FWHM F, total flux from FIR SED), seen through the 12m uv hole (u_min = 0.983/MRS) and a 10.5" taper.
# Recovered beam-centre Tb vs total; F constrained by observed dust Tb; gives missing Tb and J_bg bound.
import numpy as np
H=6.62607e-27;KB=1.380649e-16;C=2.998e10; arc=np.pi/180/3600
def B(nu,T): return 2*H*nu**3/C**2/np.expm1(H*nu/(KB*T))
nu160=C/160e-4
S=lambda nu,Td,beta: 93.48*(nu/nu160)**beta*B(nu,Td)/B(nu160,Td)
def disk(nu,F,Td,beta,mrs,b=10.5):
    om=np.pi*F**2/(4*np.log(2))*arc**2; ct=C**2/(2*KB*nu**2); Ipk=S(nu,Td,beta)*1e-23/om*ct
    sm=Ipk*F**2/(F**2+b**2)
    # smoothed peak = integral of the (Gaussian) visibility distribution; fraction at u>umin is exp(-a umin^2)
    a=np.pi**2*(F**2+b**2)*arc**2/(4*np.log(2)); umin=0.983/(mrs*arc)
    rf=np.exp(-a*umin**2)
    return Ipk,sm,sm*rf,rf
L=['Obs 10.5" total Tb (N,S) at 258 GHz: 6.3,5.3 mK; dust share ~0.5 -> dust Tb ~2.6-3.1 mK (incl. any compact dust); disk recovered Tb must be <= 3.5 mK']
for mrs in (8,11,16):
 for Td,beta in ((30,1.8),(30,1.5),(34,2.0)):
  for F in (10,15,20,30,45,60):
    Ipk,sm,rec,rf=disk(257.7e9,F,Td,beta,mrs)
    L.append(f'MRS={mrs}" Td={Td} beta={beta} F={F}": smoothed total {sm*1e3:.1f} mK, recovered {rec*1e3:.2f} mK ({rf*100:.0f}%), missing {(sm-rec)*1e3:.1f} mK; unsmoothed pk {Ipk*1e3:.1f} mK'+('  EXCLUDED' if rec*1e3>3.5 else ''))
print('\n'.join(L)); open('/workspace/out/q38/q38_filter.log','w').write('\n'.join(L))
