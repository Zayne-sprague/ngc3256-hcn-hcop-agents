# Q2 round 2, fix 1 (dense fraction) and 4 (outflow rate budget): consolidates s77a (flux systematic), s77c (kinematic decomposition), s77e (non-LTE masses) into
#  (1) a flux-scale systematic table across imagings, (2) CO outflow luminosity/mass from Michiyama+2018 Table 3, dense fraction,
#  (3) a table of M, v, R, t, Mdot, Pdot, Edot for every assumption (inclination 70/80/85 deg from the line of sight - Sakamoto+2014 70-85, fiducial 80 -, R=365 / 700 pc projected, 3 mass options, 2 velocities, lobe M v/R and thin-shell 3 M v/R).
# External inputs (assumptions): D=40 Mpc (1"=194 pc), mu=2.72, Michiyama+2018 Table 3 CO(1-0) 0.8+-0.08 (inner 2.1") / 2.2 Jy km/s (whole [-360,-180] map), alpha_CO 0.4-2.1, Sakamoto+2014 M_mol(HV,S)=2.5e7 X20, Mdot_S=50 X20 (i=80) / 25 X20 (70),
#  SFR 49+-2 (Michiyama+2020), S-nucleus SFR 1-6.8 (as quoted by GOALS-JWST arXiv:2403.14751), JWST cold Mdot<=33, Michiyama+2018 CO Mdot 40-104 and dense 1-48 Msun/yr.
import json, numpy as np
from scipy.special import erf
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
A=json.load(open(OUT+'/s77a.json')); C=json.load(open(OUT+'/s77c.json')); E=json.load(open(OUT+'/s77e.json'))
# (1) flux systematic
P('## (1) lobe wing flux (-405..-150, r=1.0", baseline B field-subtracted) in independent imagings (Jy km/s), all at 1.05" circular beam')
for mol in ('HCN','HCO+'):
    rows=[(p,A[f'{p}|{mol}|r1.0|B|blue|lobe']) for p in ('std_Xb4d','tap_Xb4d','std_Xb4b','man_Xb4b','manflat_Xb4b') if f'{p}|{mol}|r1.0|B|blue|lobe' in A]
    for p,r in rows: P(f'  {mol:5s} {p:13s}: {r["F"]:.3f} +- {r["sig"]:.3f}')
    fx=[r['F'] for p,r in rows if 'Xb4d' in p]; fb=[r['F'] for p,r in rows if 'Xb4b' in p]
    P(f'  {mol}: Xb4d std/tapered mean {np.mean(fx):.3f} (spread +-{(max(fx)-min(fx))/2/np.mean(fx)*100:.0f}%); Xb4b products mean {np.mean(fb):.3f}; Xb4b/Xb4d = {np.mean(fb)/np.mean(fx):.2f}')
    res[f'flux|{mol}']=dict(xb4d=fx,xb4b=fb)
# (2) CO
Dl=40.0; z=2780/2.998e5
def Lp(S,nu): return 3.25e7*S*(nu/(1+z))**-2*Dl**2*(1+z)**-3
LCO={'inner2.1':Lp(0.8,115.2712),'map':Lp(2.2,115.2712)}; LHCN10={'inner2.1':Lp(0.07,88.6316),'map':Lp(0.12,88.6316)}
P(f"## (2) CO(1-0) L' of South-OF-Blue [-360,-180] (Michiyama+2018 Tab.3 / text): inner 2.1\": {LCO['inner2.1']:.2e}, whole wing map: {LCO['map']:.2e} K km/s pc^2; HCN(1-0) L': {LHCN10['inner2.1']:.2e} / {LHCN10['map']:.2e}")
for al in (0.4,0.8,2.1):
    P(f'   alpha_CO={al}: M_mol(lobe, window) = {al*LCO["inner2.1"]:.2e} (inner) .. {al*LCO["map"]:.2e} (map) Msun')
P(f"   alpha route for the dense gas (what p-b/p0006 used): alpha_HCN=10 x L'HCN(1-0)={10*LHCN10['map']:.2e} Msun (map) vs M_mol(alpha_CO=0.8)={0.8*LCO['map']:.2e}: ratio {10*LHCN10['map']/(0.8*LCO['map']):.2f} -> alpha_HCN=10 gives dense > total molecular (unphysical in an unbound outflow)")
res['co']=dict(LCO=LCO,LHCN10=LHCN10)
# kinematics of the broad component
def Ev(c,s): # mean |v| and mean v over v<0 part for a gaussian N(c,s) weighted by flux
    from scipy.stats import norm
    x=np.linspace(c-6*s,c+6*s,20001); w=norm.pdf(x,c,s); return float(np.sum(w*np.abs(x))/w.sum()), float(np.sum((w*x)[x<0])/w[x<0].sum())
kin={}
for nm,(cb,sb) in {'HCN Cs':(-177.5,141.7),'HCO+ Cs':(-89.3,89.5),'common C0':(-118.2,110.4)}.items():
    a,b=Ev(cb,sb); kin[nm]=(abs(cb),a,b); P(f'  broad component {nm}: centroid {cb:.0f} sigma {sb:.0f} -> mean|v| {a:.0f}, mean v of v<0 part {b:.0f} km/s')
VSET={'window-mean (-500..-150 HCN, s75b)':283.0,'broad |v|, common C0':kin['common C0'][1],'broad centroid HCO+ (low)':89.3}
# (3) mass options
G=lambda k: 10**E[k]['MO_X1e-8'][2]
_lo=10**E['Cs own kinematics']['MO_X1e-8'][2]; _cap=10**E['Cs own kinematics']['MO_X3e-9'][2]
MASS={'HCO+ X=1e-8 (low)':_lo,'central (geometric mean of low and CO cap, X(HCO+)~5e-9)':float(np.sqrt(_lo*_cap)),'HCO+ X=3e-9 (CO cap: ~ Sakamoto M_mol(HV,S)=2.5e7 X20=1)':_cap,'HCO+ X=1e-9 / HCN X=1e-8 (high)':10**E['Cs own kinematics']['MO_X1e-9'][2],'HCN X=1e-8, HDO195 scenario':10**E['HDO195']['MH_X1e-8'][2],'HCN X=1e-8, own kinematics':10**E['Cs own kinematics']['MH_X1e-8'][2]}
for k,m in MASS.items(): P(f'  mass option {k}: {m:.2e} Msun')
res['mass']=MASS
fw=lambda c,sg,a,b: 0.5*(erf((b-c)/(np.sqrt(2)*sg))-erf((a-c)/(np.sqrt(2)*sg)))
FW={'common C0 (-118,110)':fw(-118.2,110.4,-360,-180),'HCO+ Cs (-89,90)':fw(-89.3,89.5,-360,-180),'HCN Cs (-178,142)':fw(-177.5,141.7,-360,-180)}
P('## dense fraction in the Michiyama window: fraction of the broad component inside -360..-180 km/s: '+', '.join(f'{k}: {v:.2f}' for k,v in FW.items()))
for mk in ('HCO+ X=1e-8 (low)','HCO+ X=3e-9 (CO cap: ~ Sakamoto M_mol(HV,S)=2.5e7 X20=1)'):
    for ap in ('inner2.1','map'):
        for al in (0.4,0.8,2.1):
            Md=MASS[mk]*FW['common C0 (-118,110)']; Mc=al*LCO[ap]; P(f'  M_dense(window)={Md:.1e} [{mk[:22]}] / M_CO(window, {ap}, alpha_CO={al})={Mc:.1e} -> f_dense={Md/Mc:.2f}')
P("  luminosity ratio L'(HCO+ 1-0)/L'(CO 1-0) = 0.20+-0.04 at South-OF-Blue vs 0.08+-0.01 at the S nucleus (Michiyama+2018 Table 4)")
res['FW']=FW
SEC=3.156e7; PC=3.0857e13; MSUN=1.989e33
def budget(M,v,Rp,i,shell=1.0):
    ci,si=np.cos(np.deg2rad(i)),np.sin(np.deg2rad(i)); vt=v/ci; R=Rp/si; t=R*PC/vt/SEC/1e6      # Myr
    Md=shell*M/(t*1e6); Pd=Md*MSUN/SEC*vt*1e5; Ed=0.5*Md*MSUN/SEC*(vt*1e5)**2
    return dict(vt=vt,R=R,t=t,Mdot=Md,Pdot=Pd,Edot=Ed)
P('## (3) budget rows: mass option | v_los (km/s) | R_proj (pc) | i (deg from LOS) | formula -> v_true, R_true, t[Myr], Mdot[Msun/yr], Pdot[dyn], Edot[erg/s]')
rows=[]
for mk in ('HCO+ X=1e-8 (low)','central (geometric mean of low and CO cap, X(HCO+)~5e-9)','HCO+ X=3e-9 (CO cap: ~ Sakamoto M_mol(HV,S)=2.5e7 X20=1)','HCO+ X=1e-9 / HCN X=1e-8 (high)'):
    for vk,v in VSET.items():
        for Rp in (365.,700.):
            for i in (70.,80.,85.):
                for sh,shn in ((1.0,'lobe M v/R'),(3.0,'shell 3 M v/R')):
                    b=budget(MASS[mk],v,Rp,i,sh); rows.append((mk,vk,v,Rp,i,shn,b))
for mk,vk,v,Rp,i,shn,b in rows:
    if shn=='lobe M v/R' and ((vk.startswith('broad |v|') and Rp==365.) or (vk.startswith('window') and Rp==365.)):
        P(f'  {mk:38s}| v={v:4.0f} ({vk[:12]}) | R={Rp:3.0f} | i={i:2.0f} | {shn} -> vt={b["vt"]:5.0f} R={b["R"]:4.0f} t={b["t"]:.2f} Mdot={b["Mdot"]:6.1f} Pdot={b["Pdot"]:.1e} Edot={b["Edot"]:.1e}')
Md=np.array([r[6]['Mdot'] for r in rows]); P(f'  all {len(rows)} combinations: Mdot range {Md.min():.1f} .. {Md.max():.0f} Msun/yr, median {np.median(Md):.0f}')
cen=budget(MASS['central (geometric mean of low and CO cap, X(HCO+)~5e-9)'],kin['common C0'][1],365.,80.,1.0)
P(f'  CENTRAL (M={MASS["central (geometric mean of low and CO cap, X(HCO+)~5e-9)"]:.1e}, v_los={kin["common C0"][1]:.0f}, R_proj=365, i=80, M v/R): Mdot={cen["Mdot"]:.0f} Msun/yr, Pdot={cen["Pdot"]:.1e} dyn, Edot={cen["Edot"]:.1e} erg/s, t={cen["t"]:.2f} Myr, v_true={cen["vt"]:.0f} km/s')
# marginal dependences (log-derivatives): fix others at central
def ch(**kw):
    a=dict(M=MASS['central (geometric mean of low and CO cap, X(HCO+)~5e-9)'],v=kin['common C0'][1],Rp=365.,i=80.,shell=1.0); a.update(kw); return budget(a['M'],a['v'],a['Rp'],a['i'],a['shell'])['Mdot']
P('  one-at-a-time Mdot from the central case: '+'; '.join(f'{k}={v_:.0f}' for k,v_ in dict(M_low=ch(M=MASS['HCO+ X=1e-8 (low)']),M_capCO=ch(M=MASS['HCO+ X=3e-9 (CO cap: ~ Sakamoto M_mol(HV,S)=2.5e7 X20=1)']),M_high_uncapped=ch(M=MASS['HCO+ X=1e-9 / HCN X=1e-8 (high)']),v_283=ch(v=283.),R700=ch(Rp=700.),i70=ch(i=70.),i85=ch(i=85.),shell3=ch(shell=3.0)).items()))
res['central']=cen; res['rows_range']=[float(Md.min()),float(Md.max()),float(np.median(Md))]
P('  Literature: Sakamoto+2014 CO Mdot_S = 50 X20 (i=80) / 25 X20 (i=70); Michiyama+2018 Table 2: CO 40-104, HCN 1-48 Msun/yr; GOALS-JWST cold <=33 (warm 2.0+-0.8); SFR total 49+-2 (Michiyama+2020), S nucleus 1-6.8.')
json.dump(res,open(OUT+'/s77h.json','w'),indent=1,default=float); open(OUT+'/s77h.txt','w').write('\n'.join(txt))
