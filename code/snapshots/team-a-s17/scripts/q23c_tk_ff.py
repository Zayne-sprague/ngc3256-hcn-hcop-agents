# Q23c: conditional Tk lower limits from CO(2-1) T_b: Tk_min vs filling factor and tau, plus a beam-dilution (Gaussian source) fit of T_b(beam).
# Inputs: out/q23b/vmatch_pos.json (0.18" cube convolved to circular beams), out/q23/co_hires.json native values. Output out/q23b/tk_ff.json
import numpy as np, json
from scipy.optimize import curve_fit
H=6.62607e-27;KB=1.380649e-16;NU=230.538e9;TBG=2.73
J=lambda T: H*NU/KB/(np.exp(H*NU/KB/T)-1)
Jinv=lambda j: H*NU/KB/np.log1p(H*NU/KB/j)
vp=json.load(open('out/q23b/vmatch_pos.json')); out={}
rj2pl=lambda t: t  # T_b here is RJ (K); J(T) is the RJ-equivalent of the Planck temperature in units of K: J = hnu/k /(e^x-1)
for nm in 'SN':
    bm=vp['cores'][nm]['beams']; b=np.array([float(k) for k in bm])
    pk=np.array([bm[k]['at_core_interp']['peak_RJ'] for k in bm]); hm=np.array([bm[k]['at_core_interp']['HCN']['T_mean_FWHM'] for k in bm])
    r={'beams':b.tolist(),'peak_RJ':pk.tolist(),'HCN_meanFWHM_RJ':hm.tolist()}
    f=lambda b,Ts,th: Ts*th**2/(th**2+b**2)   # Gaussian source (FWHM th) in circular Gaussian beam (FWHM b), peak-on-source
    for lab,y in [('peak',pk),('HCNmean',hm)]:
        try:
            p,c=curve_fit(f,b,y,p0=[50,0.5],bounds=([y.max(),0.05],[500,5]));e=np.sqrt(np.diag(c))
            r['fit_'+lab]={'Ts_RJ':float(p[0]),'e_Ts':float(e[0]),'theta_FWHM_arcsec':float(p[1]),'e_theta':float(e[1]),'ff_at_0.26':float(f(0.26,*p)/p[0]),'ff_at_0.57':float(f(0.57,*p)/p[0]),'ff_at_1.7':float(f(1.7,*p)/p[0]),'chi_rms':float(np.std(y-f(b,*p)))}
        except Exception as ex: r['fit_'+lab]=str(ex)
    # Tk_min table: observed RJ T_b at 0.26" (peak and HCN-velocity-mean); T_b = ff*(J(Tex)-J(Tbg))*(1-exp(-tau))
    tab={}
    for lab,T in [('peak',pk[0]),('HCNmean',hm[0])]:
        t={}
        for tau in [0.5,1,3,10,30]:
            for ff in [0.1,0.2,0.3,0.5,0.7,1.0]:
                j=T/(ff*(1-np.exp(-tau)))+J(TBG); t[f'tau{tau}_ff{ff}']=float(Jinv(j))
        tab[lab]=t
    r['Tk_min_0.26']=tab; r['T_RJ_0.26']={'peak':float(pk[0]),'HCNmean':float(hm[0])}
    out[nm]=r
    print(nm,'T_RJ(0.26) peak/HCNmean',round(pk[0],1),round(hm[0],1)); print(' fits',{k:{a:round(b,2) for a,b in v.items()} for k,v in r.items() if k.startswith('fit') and isinstance(v,dict)})
    for lab in tab: print(' ',lab,' Tk_min [rows tau 1,3,10 ; cols ff .1 .2 .3 .5 1]:',[[round(tab[lab][f'tau{tau}_ff{ff}'],0) for ff in [0.1,0.2,0.3,0.5,1.0]] for tau in [1,3,10]])
json.dump(out,open('out/q23b/tk_ff.json','w'),indent=1)
