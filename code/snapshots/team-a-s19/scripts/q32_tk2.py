# Q32 round 2: Tk lower limits per model (peak-channel, tau from Q30 R=40/60/90 central values)
import json,numpy as np,sys
sys.argv=['x']; 
H=6.62607e-27; KB=1.380649e-16; NU=230.538e9
Jf=lambda T: H*NU/KB/(np.exp(H*NU/KB/T)-1); Jinv=lambda j: H*NU/KB/np.log1p(H*NU/KB/j)
tk=lambda T,tau,ff=1.0: float(Jinv(T/(ff*(1-np.exp(-tau)))+Jf(2.73)))
d=json.load(open('/workspace/out/q32/models.json')); tt=json.load(open('/workspace/out/q30/tau_tk.json'))
out={}
for c in 'SN':
    tau={R:tt[f'{c}_ratio_peak'][R]['tau'][1] for R in ('R40','R60','R90')}
    T=d['1d'][f'{c}_peak']['T']; r={'tau':tau}
    r['rigorous_ff1_native_0.18']={R:tk(T[0],t) for R,t in tau.items()}
    r['rigorous_ff1_obs_0.26']={R:tk(T[1],t) for R,t in tau.items()}
    r['smooth_core_T0']={}
    for mn,v in d['1d'][f'{c}_peak'].items():
        if mn=='T': continue
        r['smooth_core_T0'][mn]={'T0':v['T0'],'s':v['s_fwhm'],'ff26':v['ff_0.26'],'Tk':{R:tk(v['T0'],t) for R,t in tau.items()}}
    for k,v in d['2d'].items():
        if k.startswith(f'{c}_peak') and k.endswith('nc3'): r['smooth_core_T0']['2D '+k]={'T0':v['T0'],'s':v['s_core'],'ff26':v['ff26'],'Tk':{R:tk(v['T0'],t) for R,t in tau.items()}}
    # hypothetical sub-beam clumping: ff_sub applied to the deconvolved core T0 (M1)
    T0=d['1d'][f'{c}_peak']['M1 gauss+const']['T0']
    r['subclump_grid_on_T0_M1']={f'ff{f}':{R:tk(T0,t,f) for R,t in tau.items()} for f in (1.0,0.7,0.5,0.3)}
    r['observed_0.26_ff1_over_ffgrid']={f'ff{f}':{R:tk(T[1],t,f) for R,t in tau.items()} for f in (1.0,0.9,0.77,0.6,0.4)}
    out[c]=r
    print(c,json.dumps({k:v for k,v in r.items() if k!='smooth_core_T0'},indent=0,default=lambda x:round(x,1))[:1500])
    for k,v in r['smooth_core_T0'].items(): print('  ',k,'T0=%.1f s=%.2f ff26=%.2f Tk(R40,60,90)='%(v['T0'],v['s'],v['ff26']),[round(x,1) for x in v['Tk'].values()])
json.dump(out,open('/workspace/out/q32/tk2.json','w'),indent=1)
