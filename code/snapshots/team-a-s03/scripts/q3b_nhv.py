# Q3 round 2, fixes 4+5: (4) limits on HCN/HCO+ 3-2 at 225-750 km/s (LOS) blue/red of the northern nucleus (pole-on outflow, Sakamoto+2014 v_max 750 km/s) with empirical off-source noise incl. baseline ripple;
# (5) dense-gas mass of the Soutflow blue lobe and N limits as ranges over alpha_HCN and r31, vs CO outflow masses with X20=1 and the 2x lower alpha_CO of Pereira-Santaella+2024.  Outputs out/q3b/nhv.log, nhv_*.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
import pandas as pd
rng=np.random.default_rng(7)
log=open(OUT+'nhv.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
DL=40.; z=0.0093
Lp=lambda S,nu_obs: 3.25e7*S/nu_obs**2*DL**2*(1+z)**-3          # K km/s pc^2, S in Jy km/s, nu_obs GHz
NUOBS={'HCN':F0['HCN']/(1+z),'HCO+':F0['HCO+']/(1+z)}
C_={sp:vcube(sp,'Xb4d',1.05,9.) for sp in ('HCN','HCO+')}
BASEFREE={'HCN':[(1380,2000)],'HCO+':[(1280,1680),(1880,2000)]}
def wsum_ap(sp,x,y,rad,lo,hi,order=1):
    c=C_[sp]; s=aper(c,x,y,rad)*1e3; s,_=baseline(c['v'],s,BASEFREE[sp],order); k=(c['v']>=lo)&(c['v']<hi); return s[k].sum()*c['dv']/1e3,k.sum()
# off-source positions: r 4.5-8.5" from phase centre, no emission in +-300 km/s of 2775, >=2.5" from nuclei/S-lobe
c=C_['HCN']; cand=[]
r_=np.hypot(c['dx'],c['dy'])
for _ in range(4000):
    x,y=rng.uniform(-8.5,8.5,2)
    if not 4.5<np.hypot(x,y)<8.5: continue
    if min(np.hypot(x-a,y-b) for a,b in ((0,2.55),(-0.15,-2.55),(-0.2,-4.5),(0.3,0.15)))<3.0: continue
    cand.append((x,y))
def core_flux(sp,x,y,rad): return wsum_ap(sp,x,y,rad,2625,2925)[0]
def offsource(rad,lo,hi,sp,vN,nmax=80):
    pos=[];vals=[]
    for x,y in cand:
        if core_flux('HCN',x,y,rad)>0.03 or core_flux('HCO+',x,y,rad)>0.05: continue
        if any(np.hypot(x-a,y-b)<2*rad for a,b in pos): continue
        pos.append((x,y)); vals.append(wsum_ap(sp,x,y,rad,lo+vN-2777,hi+vN-2777)[0])    # same absolute-velocity window
        if len(pos)>=nmax: break
    return np.array(vals),pos
vN=2777.; rows=[]
APN={'N (r=1")':(0.0,2.55,1.0),'N (r=1.5")':(0.0,2.55,1.5),'M18 North-OF-Blue (r=1")':(POS['M18 North-OF-Blue'][0],POS['M18 North-OF-Blue'][1],1.0)}
WINS={'225-350':(225,350),'350-750':(350,750)}
res=[]
for an,(x,y,rad) in APN.items():
    for sp in ('HCN','HCO+'):
        core=wsum_ap(sp,x,y,rad,vN-150,vN+150)[0]
        for side,sgn in (('blue',-1),('red',+1)):
            for wn,(a,b) in WINS.items():
                lo,hi=(vN-b,vN-a) if sgn<0 else (vN+a,vN+b)
                # clip to the band
                vv=C_[sp]['v']; lo=max(lo,vv.min()); hi=min(hi,vv.max())
                if hi-lo<50: continue
                I,n=wsum_ap(sp,x,y,rad,lo,hi); off,pos=offsource(rad,lo,hi,sp,vN)
                mu=np.mean(off); sg=1.4826*np.median(np.abs(off-np.median(off))); I0=I-mu
                S_lim=3*max(sg,0)*np.pi*rad**2/OM((1.05,1.05,0)) if False else 3*sg*np.pi*rad**2/(1.133*1.05**2)
                res.append(dict(ap=an,line=sp,side=side,win=wn,v_lo=lo,v_hi=hi,I_JybeamKms=I,offsrc_mean=mu,offsrc_sigma=sg,n_off=len(off),I_corr=I0,snr=I0/sg,lim3_JybeamKms=3*sg,core_JybeamKms=core,lim3_over_core=3*sg/core,S_lim3_Jykms=S_lim,Lp_lim3=Lp(S_lim,NUOBS[sp])))
R=pd.DataFrame(res); R.round(4).to_csv(OUT+'nhv_N_limits.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
pr('== N high-velocity windows (LOS velocity relative to v_N=2777; baseline: linear fit to line-free channels; sigma = robust scatter of the same absolute-velocity window over off-source apertures)')
pr(R[['ap','line','side','win','v_lo','v_hi','I_JybeamKms','offsrc_mean','offsrc_sigma','n_off','snr','lim3_JybeamKms','lim3_over_core']].round(4).to_string(index=False))
# stack N (r=1") and M18 N-OF-Blue (r=1"): mean, noise from off-source mean of the two-aperture average ~ sigma/sqrt2 only if independent; use the average of the two sigmas /sqrt(2)
st=[]
for sp in ('HCN','HCO+'):
    for side in ('blue','red'):
        for wn in WINS:
            a=R[(R.ap=='N (r=1")')&(R.line==sp)&(R.side==side)&(R.win==wn)]; b=R[(R.ap=='M18 North-OF-Blue (r=1")')&(R.line==sp)&(R.side==side)&(R.win==wn)]
            if len(a) and len(b): st.append(dict(line=sp,side=side,win=wn,mean_corr=(a.I_corr.iloc[0]+b.I_corr.iloc[0])/2,sigma=np.hypot(a.offsrc_sigma.iloc[0],b.offsrc_sigma.iloc[0])/2,))
S_=pd.DataFrame(st); S_['snr']=S_.mean_corr/S_.sigma; S_['lim3']=3*S_.sigma; pr('== stack of N(r=1") and M18 North-OF-Blue (r=1")'); pr(S_.round(4).to_string(index=False)); S_.round(4).to_csv(OUT+'nhv_N_stack.csv',index=False)
# ---- dense-gas mass limits for N (HCN, 350-750 blue+red windows, r=1.5" aperture) vs Sakamoto CO HV mass (6e7 X20)
pr('== dense-gas mass limits at N (3 sigma, HCN 3-2, r=1.5", 350-750 km/s windows; blue and red separately)')
MCO={'N':6e7}
for side in ('blue','red'):
    r=R[(R.ap=='N (r=1.5")')&(R.line=='HCN')&(R.side==side)&(R.win=='350-750')].iloc[0]; L=r.Lp_lim3
    lo=3*L/1.0; hi=10*L/0.3
    pr(f'HCN {side} 350-750: S<{r.S_lim3_Jykms:.3f} Jy km/s, L\'(3-2)<{L:.2e} K km/s pc2; M_dense < {lo:.1e} (alpha=3, r31=1) .. {hi:.1e} Msun (alpha=10, r31=0.3); fraction of CO HV mass 6e7*X20 (X20=1): {lo/6e7:.2f}..{hi/6e7:.2f}; with X20=0.5: {lo/3e7:.2f}..{hi/3e7:.2f}; wing/core flux <{r.lim3_over_core:.3f}')
# ---- Soutflow lobe masses
pr('\n== Soutflow blue lobe dense-gas luminosity and mass (HCN, HCO+; box -1.7<dRA<1.3, -6.0<dDec<-3.3; blue half = v_loc-400..v_loc-100 with v_loc=2711; also the 200-450 km/s blue range of the earlier analysis)')
box=lambda c:(c['dx']>-1.7)&(c['dx']<1.3)&(c['dy']>-6.0)&(c['dy']<-3.3)
out=[]
for wn,(lo,hi) in {'blue half (vloc-400..-100)':(2311,2611),'2300-2575 (earlier)':(2300,2575)}.items():
    for sp in ('HCN','HCO+'):
        c=C_[sp]; ppb=1.133*1.05**2/c['pix']**2; m=box(c)
        # sum of spectrum over box pixels (per pixel baseline-subtracted)
        sp_pix=c['a'][:,m]; sp_sum=sp_pix.sum(1)/ppb*1e3     # mJy (flux in box)
        s,e=baseline(c['v'],sp_sum,BASEFREE[sp],1); k=(c['v']>=lo)&(c['v']<hi); S=s[k].sum()*c['dv']/1e3
        # empirical noise of box-summed window sums: sliding sums in line-free channels
        n=k.sum(); fm=np.zeros(len(c['v']),bool)
        for a,b in BASEFREE[sp]: fm|=(c['v']>=a)&(c['v']<=b)
        sm=np.convolve(s,np.ones(n),'valid')*c['dv']/1e3; mv=np.convolve(fm.astype(float),np.ones(n),'valid')>n-0.5; eS=robust_sigma(sm[mv])
        out.append(dict(win=wn,line=sp,S_Jykms=S,eS=eS,Lp32=Lp(S,NUOBS[sp])))
O=pd.DataFrame(out); pr(O.round(4).to_string(index=False)); O.round(4).to_csv(OUT+'nhv_lobe_flux.csv',index=False)
# 1-0 from Michiyama+2018 Table 3 (South-OF-Blue: HCN 0.07+-0.02, HCO+ 0.10+-0.02 Jy km/s); L'(1-0) directly
nu10={'HCN':88.6316/(1+z),'HCO+':89.1885/(1+z)}; M18={'HCN':(0.07,0.02),'HCO+':(0.10,0.02)}
pr("\nDense mass ranges for the blue lobe, alpha 3-10 Msun/(K km/s pc2) (Gao & Solomon 2004: 10 for virialised cores; 3 as outflow-like lower value, cf. Garcia-Burillo+2012), two routes:")
rowsm=[]
for sp in ('HCN','HCO+'):
    S,eS=M18[sp]; L10=Lp(S,nu10[sp]); e10=Lp(eS,nu10[sp])
    L32=O[(O.win=='blue half (vloc-400..-100)')&(O.line==sp)].Lp32.iloc[0]
    r31=L32/L10; 
    # route B: 3-2 only: L'(1-0) between L'(3-2) (r31=1, thermalised) and L'(3-2)/0.3
    rowsm.append(dict(line=sp,Lp10_M18=L10,e_Lp10=e10,Lp32=L32,r31_implied=r31,M_route_A_a3=3*L10,M_route_A_a10=10*L10,M_route_B_min=3*L32/1.0,M_route_B_max=10*L32/0.3))
    pr(f'{sp}: L\'(1-0)(M18)={L10:.2e}+-{e10:.1e}; L\'(3-2) blue half={L32:.2e}; implied r31={r31:.2f} (2.1" M18 beam vs our box: differing apertures, M18 zero velocity unknown -> +-50%); M_dense route A (M18 1-0 x alpha 3..10) = {3*L10:.1e}..{10*L10:.1e}; route B (3-2 only, r31 0.3..1) = {3*L32/1.0:.1e}..{10*L32/0.3:.1e} Msun')
pd.DataFrame(rowsm).to_csv(OUT+'nhv_mass.csv',index=False)
L10=rowsm[0]['Lp10_M18']
for nm,Mco in (('S outflow HV CO mass 2.5e7 X20, X20=1',2.5e7),('with alpha_CO halved (Pereira-Santaella+2024, X20=0.5)',1.25e7)):
    pr(f'dense/CO fraction vs {nm}: HCN route A {3*L10/Mco:.2f}..{10*L10/Mco:.2f}; route B {rowsm[0]["M_route_B_min"]/Mco:.2f}..{rowsm[0]["M_route_B_max"]/Mco:.2f}. NB the CO mass is for the whole S outflow (both lobes, all velocities, >~ 260 km/s deprojected) while the HCN lobe is only the blue half within the box, so these fractions are only indicative (velocity-window/aperture factor ~2).')
