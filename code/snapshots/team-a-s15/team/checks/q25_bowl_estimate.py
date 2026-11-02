# back-of-envelope: uniform negative bowl implied by cumulative-flux decline between 16 and 30" (from q25_growth.log)
import numpy as np
cum={'Mich HCN':(26.16,14.42),'Mich HCO+':(44.35,30.63),'Xb55 HCN':(49.86,31.37),'Xb55 HCO+':(94.13,63.65)}
A=np.pi*(30**2-16**2); a48=np.pi*(64-16); a812=np.pi*(144-64)
F={'Mich HCN':(11.79,6.54),'Mich HCO+':(19.03,11.69),'Xb55 HCN':(23.41,11.39),'Xb55 HCO+':(42.26,23.13)}
b={k:(c[0]-c[1])/A for k,c in cum.items()}   # upper bound on bowl depth (Jy km/s/arcsec2), assumes no emission 16-30"
for s in ('Mich','Xb55'):
    r=[]
    for ap,ar in ((0,a48),(1,a812)):
        h=F[s+' HCN'][ap]+b[s+' HCN']*ar; o=F[s+' HCO+'][ap]+b[s+' HCO+']*ar; r.append(h/o)
    print(s,'bowl-corrected (max) 4-8, 8-12:',np.round(r,3),'factor',round(r[1]/r[0],3))
