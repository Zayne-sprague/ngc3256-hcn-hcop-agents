# Q2 round 2, fix 1: non-LTE (RADEX LVG grid out/q3d/grid_lvg500.npz from s62_radex_grid.py: HCN, HCO+ J=1-0..4-3, Tk 10-500, n 1e3-1e8, N 1e11.5-1e19 cm-2 at dv=150 km/s) importance-sampled posterior for the S-lobe broad component.
# Parameters (N(H2) of the clumps capped at 3e23 cm-2 in the primary run, see NH2cap): log n(H2) [1e4,1e7 primary], Tk [10,50 primary], log N_HCN/dv, log N_HCO+/dv, filling factor eta of the r=1.0" aperture (log-uniform 0.02-1), abundances X(HCN), X(HCO+) log-uniform 1e-9..1e-8 (prior only, not constrained).
# Data (ln-Gaussian): integrated brightness of the broad component of HCN 3-2 and HCO+ 3-2 in the r=1.0" aperture (Omega=pi arcsec^2; s77c scenarios, flux-scale systematic 25%), and the 3-2/1-0 brightness-temperature ratios r31 in the Michiyama+2018 window
# (HCN 0.38+-0.12, HCO+ 0.19+-0.06 from s77d), one-sided limits T(4-3)/T(3-2) (HCN <1.4, HCO+ <3: s77d 10" aggregate).  Mass: M = N_mol/X * eta*Omega*D^2 * mu m_H, mu=2.72, D=40 Mpc.
# Output: out/q2r2/s77e.txt/json, s77e_post.npz.  usage: python3 -I s77e_radex_mass.py
import json, numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); rng=np.random.default_rng(42)
d=np.load('/workspace/out/q3d/grid_lvg500.npz'); T0=np.log10(d['TK']); N0=d['LN']; C0=d['LNC']
def mk(key):
    a=d[key]; TR=a[...,0]; ok=np.isfinite(TR)&(TR>0)&(a[...,2]>0); l=np.where(ok,np.log(np.where(ok,TR,1)),np.nan)
    return [RGI((T0,N0,C0),l[...,j],bounds_error=False,fill_value=np.nan) for j in range(4)]
IH=mk('hcn'); IO=mk('hcop')
nuH=[88.6316,177.2611,265.8864,354.5055]; nuO=[89.1885,178.3751,267.5576,356.7342]
D_cm=40*3.0857e24; OM=np.pi*(1/206265.)**2; MU=2.72; MH=1.6735e-24; MSUN=1.989e33
def Tb_from_S(S,nu_GHz,om_arcsec2): return S*1e-23*2.998e10**2/(2*1.380649e-16*(nu_GHz*1e9)**2*(om_arcsec2*(1/206265.)**2))   # S in Jy km/s -> K km/s (S*1e-23 cgs per Jy; velocity unit kept)
DG=150.0
SCEN={ # name: (F_HCN, F_HCO, dv_HCN, dv_HCO) broad-component total flux (Jy km/s, r=1.0") and FWHM (km/s), from s77c.json r=1.0
 'C0 shared kinematics':(0.834,1.157,2.355*110.4,2.355*110.4),
 'Cs own kinematics':(0.993,1.099,2.355*141.7,2.355*89.5),
 'HDO100 (HCN-HDO removed)':(0.760,1.124,2.355*101.1,2.355*101.1),
 'HDO195':(0.698,1.111,2.355*96.0,2.355*96.0)}
def run(name,F_H,F_O,dvH,dvO,logn=(4,7),Tk=(10,50),ns=4000000,NH2cap=3e23,use10=True,use43=True,sysf=0.25,r31=((0.38,0.12),(0.19,0.06)),fscale=1.0,cap=None,Xh=(-9,-8),Xo=(-9,-8)):
    lT=rng.uniform(np.log10(Tk[0]),np.log10(Tk[1]),ns); ln=rng.uniform(*logn,ns)
    lXh=rng.uniform(*Xh,ns); lXo=rng.uniform(*Xo,ns)
    # log N/dv (cm-2 per km/s) uniform from 1e9 up to the value that gives N(H2)=NH2cap for the drawn abundance (physical cap on the column of the clumps: Sakamoto+2014 peak nuclear Sigma_mol 3e3 X20 Msun/pc2 = 1.9e23 H2 cm-2)
    lNH=rng.uniform(9.0,np.minimum(17.0,np.log10(NH2cap*10**lXh/dvH)));lNO=rng.uniform(9.0,np.minimum(17.0,np.log10(NH2cap*10**lXo/dvO)))
    leta=rng.uniform(np.log10(0.02),0,ns)
    pH=np.stack([lT,ln,lNH+np.log10(DG)],1); pO=np.stack([lT,ln,lNO+np.log10(DG)],1)
    TH=[np.exp(IH[j](pH)) for j in range(4)]; TO=[np.exp(IO[j](pO)) for j in range(4)]
    ok=np.isfinite(TH[2])&np.isfinite(TO[2])&np.isfinite(TH[0])&np.isfinite(TO[0])&np.isfinite(TH[3])&np.isfinite(TO[3])
    eta=10**leta
    IHm=eta*TH[2]*1.0645*dvH; IOm=eta*TO[2]*1.0645*dvO                     # K km/s beam-averaged over aperture
    IHo=Tb_from_S(F_H*fscale,265.8864*(1-2780/2.998e5),OM/(1/206265.)**2); IOo=Tb_from_S(F_O*fscale,267.5576*(1-2780/2.998e5),OM/(1/206265.)**2)
    lnL=np.zeros(ns); s=np.hypot(sysf,0.07)
    with np.errstate(all='ignore'):
        lnL+=-0.5*((np.log(IHm)-np.log(IHo))/s)**2-0.5*((np.log(IOm)-np.log(IOo))/s)**2
        if use10:
            (mh,sh),(mo,so)=r31; rH=TH[2]/TH[0]; rO=TO[2]/TO[0]
            lnL+=-0.5*((rH-mh)/sh)**2-0.5*((rO-mo)/so)**2
        if use43:
            rH4=TH[3]/TH[2]; rO4=TO[3]/TO[2]
            lnL+=np.where(rH4>1.4,-0.5*((rH4-1.4)/0.3)**2,0)+np.where(rO4>3.0,-0.5*((rO4-3.0)/0.6)**2,0)
    lnL=np.where(ok,lnL,-np.inf)
    # abundance priors
    NHt=10**lNH*dvH; NOt=10**lNO*dvO
    area=eta*OM*D_cm**2
    MH_=NHt/10**lXh*area*MU*MH/MSUN; MO_=NOt/10**lXo*area*MU*MH/MSUN
    w=np.exp(lnL-np.nanmax(lnL[ok]))
    if cap is not None: w=w*((MH_<cap)&(MO_<cap))
    w=np.where(np.isfinite(w),w,0); w/=w.sum(); ess=1/np.sum(w**2)
    def q(x,p=(2.5,16,50,84,97.5)):
        o=np.argsort(x); c=np.cumsum(w[o]); return [float(x[o][np.searchsorted(c,pp/100)]) for pp in p]
    out=dict(ess=float(ess),chi2min=float(-2*np.nanmax(lnL[ok])))
    for k,x in (('logn',ln),('Tk',10**lT),('logNH_dv',lNH),('logNO_dv',lNO),('eta',eta),('MH',MH_),('MO',MO_),('NH2_from_HCN',NHt/10**lXh),('A_NH_NO',NHt/NOt),('tauH32',None),('rH',TH[2]/TH[0]),('rO',TO[2]/TO[0]),('rH43',TH[3]/TH[2]),('rO43',TO[3]/TO[2]),('rH21',TH[1]/TH[2]),('rO21',TO[1]/TO[2])):
        if x is None: continue
        out[k]=q(np.log10(x)) if k in ('MH','MO','Tk','eta','NH2_from_HCN','A_NH_NO') else q(x)
    # variance decomposition of log M_HCN: spread from X prior alone vs rest
    out['sd_logMH_total']=float(np.sqrt(np.sum(w*(np.log10(MH_)-np.sum(w*np.log10(MH_)))**2)))
    out['sd_logX_prior']=float(np.sqrt(np.sum(w*(lXh-np.sum(w*lXh))**2)))
    # M conditional on X fixed (remove X): M*X/X0
    out['MO_X1e-9']=q(np.log10(MO_*10**lXo/1e-9)); out['MO_X1e-8']=[x-1 for x in out['MO_X1e-9']]; out['MO_X3e-9']=[x+np.log10(1/3) for x in out['MO_X1e-9']]
    out['logNmol_O']=q(np.log10(NOt*area)); out['logNmol_H']=q(np.log10(NHt*area))
    MHfix=MH_*10**lXh/1e-9; out['MH_X1e-9']=q(np.log10(MHfix)); out['MH_X1e-8']=[x-1 for x in out['MH_X1e-9']]
    # fraction of weight with n>1e6 etc
    out['P_logn_gt_6']=float(np.sum(w*(ln>6))); out['P_logn_lt_5']=float(np.sum(w*(ln<5)))
    return out,(w,ln,lT,MH_,MO_,lNH,lNO,eta)
R={}
def show(tag,o):
    f=lambda k,i=None: o[k] if i is None else o[k][i]
    P(f'{tag}: ESS={o["ess"]:.0f} chi2min={o["chi2min"]:.1f}  log n(H2) median {o["logn"][2]:.2f} [16-84: {o["logn"][1]:.2f},{o["logn"][3]:.2f}]  Tk median {10**o["Tk"][2]:.0f} K [{10**o["Tk"][1]:.0f},{10**o["Tk"][3]:.0f}]  eta {10**o["eta"][2]:.2f} [{10**o["eta"][1]:.2f},{10**o["eta"][3]:.2f}]')
    P(f'   log M_HCN median {o["MH"][2]:.2f} (95%: {o["MH"][0]:.2f}..{o["MH"][4]:.2f}; 68%: {o["MH"][1]:.2f}..{o["MH"][3]:.2f}) Msun [X(HCN) 1e-9..1e-8 prior];  at X=1e-9 fixed: {o["MH_X1e-9"][2]:.2f} (68%: {o["MH_X1e-9"][1]:.2f}..{o["MH_X1e-9"][3]:.2f}); X=1e-8: {o["MH_X1e-8"][2]:.2f};  log M_HCO+ median {o["MO"][2]:.2f} (95%: {o["MO"][0]:.2f}..{o["MO"][4]:.2f}) [X(HCO+) 1e-9..1e-8]')
    P(f'   HCO+-based mass: X(HCO+)=1e-9: {o["MO_X1e-9"][2]:.2f} (68% {o["MO_X1e-9"][1]:.2f}..{o["MO_X1e-9"][3]:.2f}), 3e-9: {o["MO_X3e-9"][2]:.2f}, 1e-8: {o["MO_X1e-8"][2]:.2f} (68% {o["MO_X1e-8"][1]:.2f}..{o["MO_X1e-8"][3]:.2f}) [log Msun]; number of HCO+ molecules log {o["logNmol_O"][2]:.2f}, HCN log {o["logNmol_H"][2]:.2f}')
    P(f'   N_HCN/N_HCO+ (log) {o["A_NH_NO"][2]:.2f} [{o["A_NH_NO"][1]:.2f},{o["A_NH_NO"][3]:.2f}];  model 3-2/1-0 T ratio HCN {o["rH"][2]:.2f}, HCO+ {o["rO"][2]:.2f}; 4-3/3-2 HCN {o["rH43"][2]:.2f} HCO+ {o["rO43"][2]:.2f}; 2-1/3-2 HCN {o["rH21"][2]:.2f} HCO+ {o["rO21"][2]:.2f}; P(n>1e6)={o["P_logn_gt_6"]:.2f}, P(n<1e5)={o["P_logn_lt_5"]:.2f}; sd(log M_HCN)={o["sd_logMH_total"]:.2f} of which X prior alone {o["sd_logX_prior"]:.2f}')
for nm,(a,b,c,dd) in SCEN.items():
    o,_=run(nm,a,b,c,dd); R[nm]=o; show(f'[{nm}] primary prior n 1e4-1e7, Tk 10-50',o)
base=SCEN['Cs own kinematics']
for tag,kw in (('wide prior n 1e3-1e8, Tk 10-300',dict(logn=(3,8),Tk=(10,300))),('Tk 10-15 (JWST outer, <15 K)',dict(Tk=(10,15))),('Tk 40-100 (JWST inner >40 K)',dict(Tk=(40,100))),('N(H2) cap 1e22 (Michiyama N_CO/1e-4 scale)',dict(NH2cap=1e22)),('N(H2) cap 1e24',dict(NH2cap=1e24)),('N(H2) cap 1e25 (no effective cap)',dict(NH2cap=1e25)),('no 1-0 constraint',dict(use10=False)),
               ('r31 from the "total" aperture (HCN 0.40, HCO+ 0.19; same)',dict(r31=((0.40,0.12),(0.19,0.06)))),('flux scale x1.25',dict(fscale=1.25)),('flux scale x0.8',dict(fscale=0.8)),
               ('mass cap 5e7 (Sakamoto HV CO mass 2.5e7 X20, X20<=2)',dict(cap=5e7)),('X(HCN) 3e-10..1e-8',dict(Xh=(-9.5,-8)))):
    o,_=run(tag,*base,**kw); R['Cs|'+tag]=o; show(f'[Cs, {tag}]',o)
# prior-only (no data) mass range for reference
json.dump(R,open(OUT+'/s77e.json','w'),indent=1,default=float); open(OUT+'/s77e.txt','w').write('\n'.join(txt))
