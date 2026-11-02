# Q1: unresolved-nuclei deblending. The Xb4d (0.94x0.74") 3-2 moment-0 map of each molecule is split into zones N (north of the N-S midpoint, r<8"),
# S (south, r<8") and Ext (r>=8", S/N>3 only), convolved to the beam of the J map being fitted, and fitted (weighted linear LS, coefficients >=0 enforced
# by lsq_linear) to the matched-beam moment-0 map of J=2-1 (Xb55), 3-2 (Xb4f, check) or 4-3 (Xb47) at 10" and 6.5" circular beams.
# Coefficient f_zone,J = (zone flux at J)/(zone flux at 3-2 Xb4d, Jy km/s). R_J(zone) = R_3(zone, Xb4d) * f_HCN/f_HCO+  (Jy km/s ratios; x1.0126 for K km/s).
# Assumes each zone's morphology at J equals the Xb4d 3-2 morphology (excitation gradients inside a zone are not modelled).
# Also: Sblue wing/core HCN/HCO+ in matched apertures (Xb4d native). Output out/s31_zones.json/.txt, out/s31_wing.txt
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s30_matchedJ_ratios import *
from scipy.optimize import lsq_linear
def templates(zone_split_dec_off=0.0):
    o=run('Xb4d',3,None,windows=(350,)); T={}
    for mol in ('HCN','HCO+'):
        R=o[mol]; m0=R['m0_350'].copy(); sig=R['err_350']; m0=np.where(m0>3*sig,m0,0.0)
        n=m0.shape[0]; Y,X=np.mgrid[:n,:n]; wc_xy=R['pxy']
        cx,cy=R['pxy'](*CEN); rr=R['rr']; dy=(Y-cy)*R['ps']-zone_split_dec_off   # +y = north
        zones={'N':(rr<8)&(dy>0),'S':(rr<8)&(dy<=0),'Ext':rr>=8}
        T[mol]={z:m0*mk for z,mk in zones.items()}; T[mol]['_meta']=(R['native'],R['ps'],R['beamarea'])
    return T
def smooth_map(m,native,ps,tgt):
    c=m[None].astype(np.float32); return conv(c,ps,native,tgt)[0]
def fit(ous,J,fwhm,T):
    o=run(ous,J,fwhm,windows=(350,)); res={}
    for mol in ('HCN','HCO+'):
        R=o[mol]; tgt=R['beam']; nat,ps_t,ba_t=T[mol]['_meta']
        # template maps are on the Xb4d grid; resample onto the J-map grid via world coordinates is avoided by working in arcsec offsets from the midpoint
        # (both grids are centred on the N-S midpoint to <1 pixel; interpolate with map_coordinates)
        zs=[z for z in ('N','S','Ext')]
        tm={z:smooth_map(T[mol][z],nat,ps_t,tgt) for z in zs}
        n=R['m0_350'].shape[0]; Y,X=np.mgrid[:n,:n]; cx,cy=R['pxy'](*CEN)
        # sample template at J-map pixel centres: offsets in arcsec
        n_t=tm['N'].shape[0]; cxt,cyt=T['cen']
        xs=(X-cx)*R['ps']/ps_t+cxt; ys=(Y-cy)*R['ps']/ps_t+cyt
        A=np.stack([map_coordinates(tm[z],[ys.ravel(),xs.ravel()],order=1,mode='constant') for z in zs],1)
        y=R['m0_350'].ravel(); e=R['err_350'].ravel(); mk=(R['rr'].ravel()<18)&(e>0)
        # convert template (Jy/beam of target beam, Xb4d flux units) -> coefficients are dimensionless f
        W=1/e[mk]; sol=lsq_linear(A[mk]*W[:,None],y[mk]*W,bounds=(0,np.inf))
        f=sol.x; chi2=np.sum((A[mk]@f-y[mk])**2*W**2)/(mk.sum()-len(f))
        infl=max(1.0,chi2)*R['beamarea']   # correlated pixels + model misfit
        cov_=np.linalg.inv((A[mk]*W[:,None]).T@(A[mk]*W[:,None]))*infl
        res[mol]=dict(f=f.tolist(),ef=np.sqrt(np.diag(cov_)).tolist(),chi2r=float(chi2),corr_NS=float(cov_[0,1]/np.sqrt(cov_[0,0]*cov_[1,1])),
                      resid_rms_over_peak=float(np.std((A[mk]@f-y[mk]))/y[mk].max()))
    return res
if __name__=='__main__':
    T=templates(); o=run('Xb4d',3,None,windows=(350,)); T['cen']=o['HCN']['pxy'](*CEN)
    # native zone fluxes and ratios from the 3-2 template itself
    Z={}; ba=T['HCN']['_meta'][2]
    for z in ('N','S','Ext'):
        Z[z]=(T['HCN'][z].sum()/ba,T['HCO+'][z].sum()/ba)
    txt=['## Xb4d 3-2 zone fluxes (Jy km/s, S/N>3 pixels, native beam, no bowl correction): zone HCN HCO+ R_Jy']
    for z,(a,b) in Z.items(): txt.append(f'  {z:4s} {a:8.2f} {b:8.2f} {a/b:.3f}')
    out={'zones_3-2_Xb4d':{z:dict(HCN=a,HCOp=b,R=a/b) for z,(a,b) in Z.items()}}
    for ous,J,fwhm in (('Xb4f',3,10.0),('Xb4f',3,6.5),('Xb55',2,10.0),('Xb47',4,10.0),('Xb47',4,6.5)):
        r=fit(ous,J,fwhm,T); tag=f'{ous}_J{J}_{fwhm:g}'
        txt.append(f'## {tag}: chi2r HCN {r["HCN"]["chi2r"]:.1f} HCO+ {r["HCO+"]["chi2r"]:.1f}; corr(N,S) of coefficients HCN {r["HCN"]["corr_NS"]:.2f} HCO+ {r["HCO+"]["corr_NS"]:.2f}')
        out[tag]=dict(fit=r,R={})
        for i,z in enumerate(('N','S','Ext')):
            fH,eH=r['HCN']['f'][i],r['HCN']['ef'][i]; fO,eO=r['HCO+']['f'][i],r['HCO+']['ef'][i]
            R3=Z[z][0]/Z[z][1]
            if fH>3*eH and fO>3*eO:
                rr_=R3*fH/fO; er=rr_*np.hypot(eH/fH,eO/fO); rel=fH/fO
            else: rr_=er=rel=float('nan')
            txt.append(f'  {z:4s} f_HCN {fH:.3f}+-{eH:.3f}  f_HCO+ {fO:.3f}+-{eO:.3f}  (J/3-2 flux ratio incl. recovery)  R_J(Jy)={rr_:.3f}+-{er:.3f}  [R_3(Xb4d)={R3:.3f}]')
            out[tag]['R'][z]=dict(f_HCN=fH,ef_HCN=eH,f_HCOp=fO,ef_HCOp=eO,R=rr_,eR=er,R3=R3)
        print('\n'.join(txt[-4:]),flush=True)
    json.dump(out,open('/workspace/out/s31_zones.json','w'),indent=1,default=float)
    open('/workspace/out/s31_zones.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
