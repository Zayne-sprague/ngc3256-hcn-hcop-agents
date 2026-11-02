import numpy as np
from scipy.optimize import curve_fit
from q4_lib import robust_sigma
def gauss(v,A,v0,fw,a,b): return A/(fw/2.3548*2.5066)*np.exp(-0.5*((v-v0)/(fw/2.3548))**2)+a+b*(v-v0)/100.
def bright(v,y,v00=2780.,fw0=150.):
    m=np.abs(v-v00)<700
    p,cv=curve_fit(gauss,v[m],y[m],p0=[y[m].max()*fw0*2.5/2.3548*0.9,v00,fw0,0,0],maxfev=20000)
    r=y[m]-gauss(v[m],*p); far=np.abs(v[m]-p[1])>1.5*p[2]; sg=robust_sigma(r[far]) if far.sum()>8 else robust_sigma(r)
    dv=abs(np.median(np.diff(v))); return p,np.sqrt(abs(cv[0,0])),sg
