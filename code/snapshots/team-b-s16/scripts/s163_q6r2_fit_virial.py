"""Q6 r2: T_b(scale) power-law and core+pedestal fits (from s160), virial-mass variants (from s151, s162). python3 -I s163_q6r2_fit_virial.py"""
import json, numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
t=json.load(open('/workspace/out/q6r2/s160.json')); s151=json.load(open('/workspace/out/q6/s151.json')); s162=json.load(open('/workspace/out/q6r2/s162.json'))
PC=194.; out={'tb':{},'vir':{}}; fig,ax=plt.subplots(1,2,figsize=(10,4),sharey=True)
keys=['native 0.24"','55pc(0.28")','120pc(0.62")','250pc(1.29")','Xb4d beam 0.93x0.74"']
W=1.8
print('T_b(theta) (K RJ): native | 55 | 120 | 250 pc | slope(0.24->1.29") | core+ped fit: Tc, d, Tp (w=%.1f") | Tb(Xb4d native)/Tb(Xb4b smoothed to Xb4d beam)'%W)
for li,ln in enumerate(t):
    for rn,res in t[ln].items():
        th=np.array([res[k]['theta'] for k in keys]); T=np.array([res[k]['Tb'] for k in keys]); e=np.array([max(res[k]['sig'],0.03*res[k]['Tb']) for k in keys])
        sl=np.log(res[keys[3]]['Tb']/res[keys[0]]['Tb'])/np.log(res[keys[3]]['theta']/res[keys[0]]['theta'])
        sl2=np.log(res[keys[2]]['Tb']/res[keys[0]]['Tb'])/np.log(res[keys[2]]['theta']/res[keys[0]]['theta'])
        f=lambda p:(p[0]*p[2]**2/(p[2]**2+th**2)+p[1]*W**2/(W**2+th**2)-T)/e
        best=None
        for d0 in (0.15,0.3,0.5,0.8):
            r=least_squares(f,[T[0],0.5,d0],bounds=([0,0,0.02],[100,100,2.0]))
            if best is None or r.cost<best.cost: best=r
        Tc,Tp,d=best.x; J=best.jac; cv=np.linalg.pinv(J.T@J); chi2=2*best.cost
        dil=lambda th0: d**2/(d**2+th0**2)
        res_=dict(slope_native_to_250pc=float(sl),slope_native_to_120pc=float(sl2),Tc=float(Tc),Tp=float(Tp),d_core=float(d),chi2=float(chi2),eTc=float(np.sqrt(cv[0,0])),
                  filling_core_in_native_beam=float(dil(0.239)),ratio_Xb4d_native_over_smoothed=float(res['Xb4d native']['Tb']/res[keys[4]]['Tb']),T=T.tolist(),theta=th.tolist())
        out['tb'].setdefault(ln,{})[rn]=res_
        print(ln,rn,' | '.join('%.2f'%x for x in T[:4]),'| slope %.2f (to 120pc %.2f)'%(sl,sl2),'| Tc %.1f d %.2f Tp %.2f chi2 %.1f'%(Tc,d,Tp,chi2),'| %.2f'%res_['ratio_Xb4d_native_over_smoothed'])
        ax[li].errorbar(th,T,e,marker='o',label=rn); ax[li].set_title(ln)
for a in ax: a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('beam geometric-mean FWHM [arcsec] (1" = 194 pc)'); a.axhline(0,c='k')
ax[0].set_ylabel('peak T_b (RJ) [K]'); ax[1].legend(); plt.tight_layout(); plt.savefig('/workspace/out/q6r2/tb_vs_scale.png',dpi=110)
# virial variants
print('\nvirial: region line | A: R=1.91 sig_r (D_hm), sigma_1G | B: R=D_hm/2 | C: sigma_mom2 | C2: sigma_mom2 grad-corrected | D: R from radial D | A@35Mpc')
Mf=lambda sig,R: 1040*sig**2*R
for ln in ('HCO+3-2','HCN3-2'):
    for rn in ('N','S','Bridge','E1','E2'):
        L=s151['lines'][ln][rn]; D=L['size']['Ddec']; sp=L['spec']['peakpix']; g=sp['gauss1']; sg=(g['fwhm'] if (g and g['chi2r']<2) else sp['fwhm_halfmax'])/2.355
        Dr=s162['lines'][ln][rn]['regular']['radial_Ddec']
        R191=lambda D: 1.91*D/2.355*PC
        v=dict(sigma_1G=sg,sigma_mom2=sp['sigma_mom2'],sigma_mom2_int=sp.get('sigma_int_mom2'),D_hm=D,D_radial=Dr,
               A=Mf(sg,R191(D)),B=Mf(sg,D/2*PC),C=Mf(sp['sigma_mom2'],R191(D)),C2=Mf(sp.get('sigma_int_mom2') or np.nan,R191(D)),Dv=Mf(sg,R191(Dr)),A35=Mf(sg,R191(D)*170/194))
        out['vir'].setdefault(ln,{})[rn]=v
        print(ln,rn,'sig1G %.0f mom2 %.0f int %.0f | A %.1e B %.1e C %.1e C2 %.1e D %.1e A35 %.1e'%(sg,sp['sigma_mom2'],v['sigma_mom2_int'] or np.nan,v['A'],v['B'],v['C'],v['C2'],v['Dv'],v['A35']))
json.dump(out,open('/workspace/out/q6r2/s163.json','w'),indent=1)
