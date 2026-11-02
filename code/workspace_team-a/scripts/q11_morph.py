# fit sky morphology (cores+exp disk) and transfer-function scale f to the 2-1 HCO+ and 3-2 HCO+ (Xb4f) annulus profiles jointly (T from physical scaling u_a=f*Bmin/lambda)
import sys; sys.path.insert(0,'/workspace/scripts')
from q11_forward import *
from scipy.optimize import differential_evolution
obsH={}
for tag,k in {'21':'HCO+ 2-1','32f':'HCO+ 3-2@Xb4f','32d':'HCO+ 3-2@Xb4d'}.items(): obsH[tag]=obs_annuli(k)[0]
print(obsH)
SETS={'21':(177.3,55.),'32f':(265.9,37.),'32d':(265.9,22.)}
def model(p,tag,Bm=8.7,ramp=2.):
    h,sc,fc,f=p[:4]; fr,pbw=SETS[tag]; return annuli(observe(sky_hco(h,sc,fc),Tf(fr,Bm,f,ramp),pbw))
for ramp in (2.,4.):
  for tag in ('21','32f'):
    y=obsH[tag]; s=0.05*np.abs(y)+1.
    def chi(p):
        m=model(p,tag,ramp=ramp); A=(m*y/s**2).sum()/(m*m/s**2).sum(); return (((A*m-y)/s)**2).sum()
    r=differential_evolution(chi,[(3,30),(2,10),(0.02,0.9),(0.2,3)],seed=1,maxiter=60,popsize=12,tol=1e-6)
    m=model(r.x,tag,ramp=ramp); A=(m*y/(0.05*np.abs(y)+1)**2).sum()/(m*m/(0.05*np.abs(y)+1)**2).sum()
    print('ramp',ramp,tag,'chi2 %.1f'%r.fun,'h,sc,fc,f=',r.x.round(2),'model',(A*m).round(1),'obs',y.round(1))
