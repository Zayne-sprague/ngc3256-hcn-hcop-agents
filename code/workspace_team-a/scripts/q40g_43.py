# Q40g: HCN/HCO+ 4-3 (Xb47, Band 7, native 4.1x2.7") at the clumps with the same decomposition at the common 8.64x6.0 beam (native Xb47 smoothed): amplitudes, noise from null windows, 3-sigma HCN upper limits, R upper limits.
import sys; sys.path.insert(0,'/workspace/scripts')
import q40b_fit as B
import numpy as np
from q1_lib import *
from scipy.optimize import nnls
B.NU[4]=(354.505478,356.734223); B.KF[4]=(356.734223/354.505478)**2
KEY=['HCN 4-3','HCO+ 4-3']
def getm(name,ln):
    d=load(KEY[ln]); return d,d['m0'],[n[1] for n in d['nulls']]
B.getmaps=lambda name,ln: getm(name,ln) if name=='Xb47' else B.__dict__['_gm'](name,ln)
log=[]
def P(s=''): print(s); log.append(s)
d,m,ns=getm('Xb47',0); P(f'Xb47 HCN 4-3 beam {d["beam"]} pix {d["pix"]:.2f} map shape {m.shape} n_null {len(ns)} null-window rms of smoothed maps: ')
Z={}
for ln,lab in ((0,'HCN'),(1,'HCO+')):
    d,_,_=getm('Xb47',ln); pr=noise_profile(d,7.2); Z[lab]=pr
    P(f'{lab} 4-3 empirical sigma of the 7.2"-smoothed moment-0 (Jy/beam km/s, pbcor, 25th-percentile null window) by radius: '+'; '.join(f'{a}-{b}": {s:.3f}' for a,b,s,_,_ in pr))
sg=lambda lab,r: sigma_at(Z[lab],r)
for variant,sz in (('unresolved',{}),('nuc3',{'N':3,'S':3}),('nuc3+clump4',{'N':3,'S':3,'E1':4,'E2':4,'W':4})):
    o={lab:B.run('Xb47',ln,'b86',{},sz,noise=False) for ln,lab in ((0,'HCN'),(1,'HCO+'))}
    P(f'## {variant}: amplitudes Jy km/s @8.64x6.0 (value +- sigma_pix at the clump radius)')
    for reg in ('N','S','E1','E2','W'):
        i=o['HCN']['names'].index(reg); ex,ey=B.EXY[reg]; r=np.hypot(ex,ey)
        h=(o['HCN']['amp'][i],sg('HCN',r)); q=(o['HCO+']['amp'][i],sg('HCO+',r)); ul=max(h[0],0)+3*h[1]
        P(f'  {reg} (r={r:.1f}"): HCN {h[0]:.2f}+-{h[1]:.2f}  HCO+ {q[0]:.2f}+-{q[1]:.2f}  HCN 3sig UL {ul:.2f}  -> R43 (T_B) '+(f'= {h[0]/q[0]*B.KF[4]:.2f}+-{B.KF[4]*h[0]/q[0]*np.hypot(h[1]/h[0],q[1]/q[0]):.2f}' if h[0]>3*h[1] else f'< {ul/(max(q[0]-3*q[1],1e-3))*B.KF[4]:.2f} (3 sigma HCN UL / (HCO+ - 3 sigma))'))
open('/workspace/out/q40b/q40g_43.log','w').write('\n'.join(log))
