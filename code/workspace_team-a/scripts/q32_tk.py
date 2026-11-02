# Q32d: ff range from the beam-dilution fit -> Tk lower limits on the intrinsic compact component; missing-flux variant
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import curve_fit
H=6.62607e-27; KB=1.380649e-16; NU=230.538e9; TBG=2.73
J=lambda T: H*NU/KB/(np.exp(H*NU/KB/T)-1); Jinv=lambda j: H*NU/KB/np.log1p(H*NU/KB/j)
d=json.load(open('out/q23/co_hires.json')); cf=json.load(open('out/q32/curve_fit.json')); tt=json.load(open('out/q30/tau_tk.json'))
model=lambda t,Te,T0,s: Te+T0*s**2/(s**2+t**2)
th=np.array([np.sqrt(0.1805*0.1739),.26,.38,.57,.75,1.0,1.7]); keys=['native','0.26','0.38','0.57','0.75','1.0','1.7']
tk=lambda T,tau,ff: float(Jinv(T/(ff*(1-np.exp(-tau)))+J(TBG)))
out={}
for c in 'SN':
    T=np.array([d['cores'][c]['beams'][k]['Tb_RJ'] for k in keys])
    # variant: 12m-only deficit of 11% at >=0.57" (Q23b) corrected
    Tc=T.copy(); Tc[3:]*=1.11
    r={}
    for lab,TT in [('as_measured',T),('missflux_corrected',Tc)]:
        p,cv=curve_fit(model,th,TT,p0=[3,30,.3],sigma=0.08*TT+0.3,bounds=([0,1,.01],[30,300,3]))
        s=p[2]; r[lab]={'Te':p[0],'T0':p[1],'s':s,'ff_0.26':s**2/(s**2+.26**2),'ff_0.18':s**2/(s**2+th[0]**2),'ff_1.0':s**2/(s**2+1.)}
    # Tk: intrinsic compact T0 (ff=1 at its own size), tau from Q30 (R=40,60,90 central), peak ratio and HCN-mean scaling
    ratio_mean_to_peak=0.644 if c=='S' else None
    tab={}
    for R in ('R40','R60','R90'):
        tau=tt[f'{c}_ratio_peak'][R]['tau'][1]
        tab[R]={'tau':tau,'Tk_T0':tk(r['as_measured']['T0'],tau,1.0),'Tk_T0_low':tk(r['as_measured']['T0']-3,tau,1.0),
                'Tk_native_ff_fit':tk(T[0],tau,r['as_measured']['ff_0.18']),'Tk_0.26_ff1':tk(T[1],tau,1.0),'Tk_0.26_fffit':tk(T[1],tau,r['as_measured']['ff_0.26'])}
    r['Tk']=tab; out[c]=r; print(c,json.dumps(r,default=float,indent=0)[:1800])
json.dump(out,open('out/q32/tk_ff.json','w'),indent=1,default=float)
