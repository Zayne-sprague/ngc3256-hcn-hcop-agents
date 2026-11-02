# Q1 (round 2): N/S/Ext deblend of the 10" (and 6.5") HCN and HCO+ J=2-1 (Xb55), 3-2 (Xb4f, same-beam reference) and 4-3 (Xb47) moment-0 maps.
# Templates: LINE-SPECIFIC (HCN template from the Xb4d HCN 3-2 cube, HCO+ from HCO+), zones N (dec>mid), S (dec<=mid) within r<8", Ext (r>=8", S/N>3), each convolved to the target beam.
# Model: y = fN*TN + fS*TS + fE*TE + c0 + c1*G20 ; ALL coefficients signed (WLS, noise = per-pixel line-free rms, pixel-correlation variance inflation = beam area in pixels,
# further x chi2_red when chi2_red>1: misfit inflation). c0 = flat zero-level (bowl) term, G20 = broad Gaussian (FWHM 20" about the N-S midpoint, a smooth extended/bowl term).
# Template-shape systematics: 150 bootstrap draws: N,S templates shifted by U(+-0.3") each, N/S split line shifted U(+-1") in Dec, Ext scaled U(0.9,1.1), fit radius +-2".
# The SAME perturbed templates are used for all arrays in a draw so that template errors cancel in ladder ratios f(J)/f(Xb4f 3-2).
# Injection test: synthetic maps built from unclipped templates with a different split/shift (truth), known f, a bowl, and realistic correlated noise (random-sign sums of
# the array's own line-free channels); refit with the nominal templates.   Output: out/s37_deblend2.json/.txt
import sys, json; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
from s30_matchedJ_ratios import *
from scipy.ndimage import shift as ndshift
rng=np.random.default_rng(7)
S.HALF=24.0
JOBS=[('Xb4f',3,10.0,12),('Xb55',2,10.0,16),('Xb47',4,10.0,10),('Xb4f',3,6.5,12),('Xb47',4,6.5,10)]
FRAD={f'{a}_J{b}_{c:g}':r for a,b,c,r in JOBS}
o4=S.run('Xb4d',3,None,windows=(350,)); cen_t=o4['HCN']['pxy'](*CEN)
def zones(mol,split=0.0,sh=((0,0),(0,0)),clip=True,escale=1.0):
    R=o4[mol]; m0=R['m0_350'].copy(); sig=R['err_350']
    if clip: m0=np.where(m0>3*sig,m0,0.0)
    n=m0.shape[0]; Y,X=np.mgrid[:n,:n]; dy=(Y-cen_t[1])*R['ps']-split; rr=R['rr']
    Z={'N':m0*((rr<8)&(dy>0)),'S':m0*((rr<8)&(dy<=0)),'Ext':m0*(rr>=8)*escale}
    for z,(dx,dyy) in zip(('N','S'),sh): Z[z]=ndshift(Z[z],(dyy/R['ps'],dx/R['ps']),order=1)
    return Z
def smooth(m,tgt,mol):
    R=o4[mol]; return S.conv(m[None].astype(np.float32),R['ps'],R['native'],tgt)[0]
ARR={}
for ous,J,fw,_ in JOBS:
    o=S.run(ous,J,fw,windows=(350,)); ARR[f'{ous}_J{J}_{fw:g}']=o
def design(R,Zs,mol,tgt):
    tm={z:smooth(Zs[z],tgt,mol) for z in Zs}; n=R['m0_350'].shape[0]; Y,X=np.mgrid[:n,:n]; cx,cy=R['pxy'](*CEN)
    ps_t=o4[mol]['ps']; xs=(X-cx)*R['ps']/ps_t+cen_t[0]; ys=(Y-cy)*R['ps']/ps_t+cen_t[1]
    cols=[map_coordinates(tm[z],[ys.ravel(),xs.ravel()],order=1,mode='constant') for z in ('N','S','Ext')]
    cols.append(np.ones(n*n)); G=np.exp(-4*np.log(2)*((X-cx)**2+(Y-cy)**2)*R['ps']**2/20.0**2).ravel(); cols.append(G)
    return np.stack(cols,1)
def wls(A,y,e,mk,ba,infl_chi=True):
    W=1/e[mk]; Aw=A[mk]*W[:,None]; yw=y[mk]*W; f,*_=np.linalg.lstsq(Aw,yw,rcond=None)
    chi=np.sum((Aw@f-yw)**2)/(mk.sum()-A.shape[1]); cov=np.linalg.inv(Aw.T@Aw)*ba*(max(1,chi) if infl_chi else 1)
    return f,cov,chi
def fit_one(key,mol,Zs,rfit):
    R=ARR[key][mol]; tgt=R['beam']; A=design(R,Zs,mol,tgt)
    y=R['m0_350'].ravel(); e=R['err_350'].ravel(); mk=(R['rr'].ravel()<rfit)&(e>0)
    return wls(A,y,e,mk,R['beamarea']),A,y,mk
txt=['## deblend2: coefficients f (flux at J / Xb4d 3-2 template flux, Jy km/s units) per zone; err = noise x pixel-corr x chi2r inflation; boot = template-shape scatter']
out={}; nom={}
for key in ARR:
    nom[key]={}
    for mol in ('HCN','HCO+'):
        (f,cov,chi),A,y,mk=fit_one(key,mol,zones(mol),FRAD[key])
        nom[key][mol]=dict(f=f,ef=np.sqrt(np.diag(cov)),chi=chi,corr=cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov))))
        txt.append(f'{key:14s} {mol:5s} chi2r {chi:6.1f} rfit {FRAD[key]}"  f_N {f[0]:.3f}+-{np.sqrt(cov[0,0]):.3f} f_S {f[1]:.3f}+-{np.sqrt(cov[1,1]):.3f} f_Ext {f[2]:.3f}+-{np.sqrt(cov[2,2]):.3f} c0 {f[3]:.3f} c1 {f[4]:.3f}  corr(N,S) {nom[key][mol]["corr"][0,1]:.2f} corr(N,Ext) {nom[key][mol]["corr"][0,2]:.2f} corr(S,Ext) {nom[key][mol]["corr"][1,2]:.2f}')
print('\n'.join(txt),flush=True)
# bootstrap
NB=150; boot={k:{m:[] for m in ('HCN','HCO+')} for k in ARR}
for b in range(NB):
    split=rng.uniform(-1,1); sh=[tuple(rng.uniform(-.3,.3,2)) for _ in range(2)]; esc=rng.uniform(.9,1.1); dr=rng.uniform(-2,2)
    for mol in ('HCN','HCO+'):
        Zs=zones(mol,split,sh,True,esc)
        for key in ARR:
            (f,cov,chi),*_=fit_one(key,mol,Zs,FRAD[key]+dr); boot[key][mol].append(f[:3])
    if b%30==0: print('boot',b,flush=True)
boot={k:{m:np.array(v) for m,v in d.items()} for k,d in boot.items()}
ZN=['N','S','Ext']
# nominal 3-2 zone fluxes (template totals) for R_J
ba=o4['HCN']['beamarea']; T3={m:{z:zones(m)[z].sum()/ba for z in ZN} for m in ('HCN','HCO+')}
txt.append('## Xb4d 3-2 template zone fluxes (Jy km/s): '+'; '.join(f'{z} HCN {T3["HCN"][z]:.1f} HCO+ {T3["HCO+"][z]:.1f} R {T3["HCN"][z]/T3["HCO+"][z]:.3f}' for z in ZN))
res={}
def kfac(J): return (265.886/ {2:177.261,3:265.886,4:354.505}[J])**2   # T_b ratio factor (nu3/nuJ)^2 for HCN
for key in ARR:
    if key.startswith('Xb4f'): continue
    ref='Xb4f_J3_'+key.split('_')[2]; J=int(key.split('_')[1][1:]); res[key]={}
    txt.append(f'## {key} (ref {ref}): per zone')
    for i,z in enumerate(ZN):
        d={}
        for mol in ('HCN','HCO+'):
            f=nom[key][mol]['f'][i]; ef=nom[key][mol]['ef'][i]; f3=nom[ref][mol]['f'][i]; e3=nom[ref][mol]['ef'][i]
            L=f/f3 if (f>0 and f3>0) else np.nan
            bl=boot[key][mol][:,i]/boot[ref][mol][:,i]; bl=bl[np.isfinite(bl)&(boot[ref][mol][:,i]>0)]
            eL_st=abs(L)*np.hypot(ef/f,e3/f3) if np.isfinite(L) else np.nan
            d[mol]=dict(f=f,ef=ef,f3=f3,ef3=e3,ladder=L,ladder_err_stat=eL_st,ladder_boot_sd=float(np.std(bl)) if len(bl)>5 else np.nan,boot_f_sd=float(np.std(boot[key][mol][:,i])))
        # R_J = R_3(template) * fH/fO  and also ratio of ladder factors
        fH,fO=d['HCN']['f'],d['HCO+']['f']
        bR=(boot[key]['HCN'][:,i]/boot[key]['HCO+'][:,i])*(T3['HCN'][z]/T3['HCO+'][z]); bR=bR[np.isfinite(bR)&(bR>0)&(bR<5)]
        RJ=T3['HCN'][z]/T3['HCO+'][z]*fH/fO if fH>0 and fO>0 else np.nan
        eRst=abs(RJ)*np.hypot(d['HCN']['ef']/fH,d['HCO+']['ef']/fO) if np.isfinite(RJ) else np.nan
        eRb=float(np.std(bR)) if len(bR)>5 else np.nan
        # N-S anticorrelation in the bootstrap
        cNS=float(np.corrcoef(boot[key]['HCN'][:,0],boot[key]['HCN'][:,1])[0,1])
        d['R_J']=RJ; d['R_J_stat']=eRst; d['R_J_boot']=eRb; d['R_J_tot']=float(np.hypot(eRst,eRb)) if np.isfinite(eRst) and np.isfinite(eRb) else np.nan; d['boot_corr_NS_HCN']=cNS
        d['K']=kfac(J)
        res[key][z]=d
        txt.append(f'  {z:4s} f_HCN {d["HCN"]["f"]:.3f}+-{d["HCN"]["ef"]:.3f}(boot {d["HCN"]["boot_f_sd"]:.3f})  f_HCO+ {d["HCO+"]["f"]:.3f}+-{d["HCO+"]["ef"]:.3f}(boot {d["HCO+"]["boot_f_sd"]:.3f})  R_J(Jy)={RJ:.3f}+-{eRst:.3f}(stat+misfit) +-{eRb:.3f}(template)  '
                   f'ladder J/3 [Jy] HCN {d["HCN"]["ladder"]:.2f}+-{d["HCN"]["ladder_err_stat"]:.2f}/{d["HCN"]["ladder_boot_sd"]:.2f}  HCO+ {d["HCO+"]["ladder"]:.2f}+-{d["HCO+"]["ladder_err_stat"]:.2f}/{d["HCO+"]["ladder_boot_sd"]:.2f}  [K: x{kfac(J):.2f}]')
# injection
txt.append('## injection test (truth: unclipped templates, split +0.8", N shifted +0.3" E, true f below, flat bowl -0.3*peak/100, noise from own line-free channels), 40 draws each')
inj={}
for key in ('Xb55_J2_10','Xb47_J4_10'):
    mol='HCN'; R=ARR[key][mol]; ftrue=np.array([0.8,0.5,0.7]) if 'J2' in key else np.array([0.5,0.3,0.4])
    Zt=zones(mol,0.8,((0.3,0),(0,0)),clip=False)
    A=design(R,Zt,mol,R['beam']); sig_true=A[:,:3]@ftrue
    cube=R['cube_lf']; nlf=cube.shape[0]; ncore=int((abs(R['v'])<350).sum()); rec=[]
    for k in range(40):
        idx=rng.choice(nlf,ncore,replace=True); sg=rng.choice([-1,1],ncore)
        noise=(cube[idx]*sg[:,None,None]).sum(0).ravel()*R['dv']
        y=sig_true+noise-0.3; e=R['err_350'].ravel(); mk=(R['rr'].ravel()<FRAD[key])&(e>0)
        An=design(R,zones(mol),mol,R['beam']); f,cov,chi=wls(An,y,e,mk,R['beamarea']); rec.append(f[:3])
    rec=np.array(rec); inj[key]=dict(true=ftrue.tolist(),mean=rec.mean(0).tolist(),sd=rec.std(0).tolist(),corrNS=float(np.corrcoef(rec[:,0],rec[:,1])[0,1]))
    txt.append(f'  {key}: true N,S,Ext {ftrue}  recovered mean {np.round(rec.mean(0),3)} sd {np.round(rec.std(0),3)}  corr(N,S) {inj[key]["corrNS"]:.2f}')
json.dump(dict(nominal={k:{m:dict(f=v['f'],ef=v['ef'],chi=v['chi']) for m,v in d.items()} for k,d in nom.items()},res=res,inj=inj,T3=T3),open('/workspace/out/s37_deblend2.json','w'),indent=1,default=lambda x:x.tolist() if hasattr(x,'tolist') else float(x))
open('/workspace/out/s37_deblend2.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
