# (1) H13CN 4-3 vs CO 3-2 wing: template fit on the CO-free (red, low-frequency) side; (2) CS 7-6 / H2CS 10(0,10)-9(0,9) two-line template fit
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
C=299792.458; V=2780.
d=np.load('/workspace/out/region_spectra.npz')
def vel(fr,f0): return (1-fr/f0)*C
out={}
for reg in ['N','S','Bridge']:
    f16=d['Xb47__spw16__freq']; fx=d['Xb47__spw18__freq']
    T_hcn=d[f'Xb47__spw16__{reg}__pk']; v_hcn=vel(f16,354.505478)
    T_hcop=d[f'Xb47__spw18__{reg}__pk']; v_hcop=vel(fx,356.734223)
    # ---- H13CN
    f20=d['Xb47__spw20__freq']; s20=d[f'Xb47__spw20__{reg}__pk']; v20=vel(f20,345.339760)
    res={}
    for tn,(vt,T) in {'HCN4-3':(v_hcn,T_hcn),'HCO+4-3':(v_hcop,T_hcop)}.items():
        Ti=interp1d(vt-V,T,bounds_error=False,fill_value=0.)
        m=(v20-V>-20)&(v20-V<300)   # red half (low-freq side) of H13CN: CO 3-2 wing is on the blue (high-freq) side
        x=Ti(v20[m]-V); y=s20[m]
        # fit y = alpha*T + b with noise from first 150 channels of the spectrum (line-free)
        sg=np.std(s20[:150])
        A=np.vstack([x,np.ones_like(x)]).T; coef,_,_,_=np.linalg.lstsq(A,y,rcond=None)
        cov=sg**2*np.linalg.inv(A.T@A); al,e=coef[0],np.sqrt(cov[0,0])
        # same fit, b fixed 0
        al0=np.sum(x*y)/np.sum(x*x); e0=sg/np.sqrt(np.sum(x*x))
        # integrals
        IT=np.sum(Ti(v20-V))*abs(np.median(np.diff(v20)))
        res[tn]=dict(alpha=float(al),e_alpha=float(e),alpha_b0=float(al0),e_alpha_b0=float(e0),I_template_JyBkms=float(IT),I_H13CN_JyBkms=float(al*IT),e_I=float(e*IT),
                     ratio_H13CN_over_template_3sig_upper=float(max(al,0)+3*e),sigma_chan=float(sg))
    out[f'H13CN_{reg}']=res
    # ---- CS / H2CS
    f22=d['Xb47__spw22__freq']; s22=d[f'Xb47__spw22__{reg}__pk']; v22=vel(f22,342.883000)
    shift=C*(1-342.944369*(1-V/C)/342.883000)-V     # apparent velocity offset of H2CS in CS frame (km/s), negative=bluer
    sg=np.std(np.r_[s22[:100],s22[-100:]]) 
    m=(np.abs(v22-V)<500)
    for tn,(vt,T) in {'HCN4-3':(v_hcn,T_hcn),'HCO+4-3':(v_hcop,T_hcop)}.items():
        Ti=interp1d(vt-V,T,bounds_error=False,fill_value=0.)
        def model(p,withH):
            a,h,dv,b=p; return a*Ti(v22[m]-V-dv)+(h*Ti(v22[m]-V-shift-dv) if withH else 0)+b
        r0=least_squares(lambda p:(model([p[0],0,p[1],p[2]],False)-s22[m])/sg,[0.1,0,0],bounds=([0,-60,-0.05],[5,60,0.05]))
        r1=least_squares(lambda p:(model(p,True)-s22[m])/sg,[0.1,0.01,0,0],bounds=([0,0,-60,-0.05],[5,5,60,0.05]))
        J=r1.jac; cov=np.linalg.inv(J.T@J) if np.linalg.cond(J.T@J)<1e12 else np.full((4,4),np.nan)
        IT=np.sum(Ti(v22-V))*abs(np.median(np.diff(v22)))
        out[f'CS_H2CS_{reg}_{tn}']=dict(shift_H2CS_kms=float(shift),chi2_CSonly=float(2*r0.cost),chi2_CS_H2CS=float(2*r1.cost),ndata=int(m.sum()),
            a_CS=float(r1.x[0]),h_H2CS=float(r1.x[1]),e_h=float(np.sqrt(cov[1,1])),vshift=float(r1.x[2]),I_CS=float(r1.x[0]*IT),I_H2CS=float(r1.x[1]*IT),e_I_H2CS=float(np.sqrt(cov[1,1])*IT),
            I_CS_if_noH2CS=float(r0.x[0]*IT),sigma=float(sg))
json.dump(out,open('/workspace/out/blend_fits.json','w'),indent=1)
for k,v in out.items():
    print(k)
    for kk,vv in (v.items() if 'CS' in k else [(t,w) for t,w in v.items()]): print('   ',kk,vv if not isinstance(vv,float) else round(vv,4))
