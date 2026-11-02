# Q2 step d: clean-mask coverage of the wing, dense-gas luminosity/mass and (projected) mass-outflow rate of the S blue lobe and the red candidate.
# Inputs: out/q2/s75b.json (Xb4d@1.05" fluxes). External assumptions are in ASSUME below and printed with every number.
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
OUT='/workspace/out/q2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
J=json.load(open(OUT+'/s75b.json'))
S0=REGIONS['S'][:2]
# ---- (1) clean-mask coverage of the blue wing at native resolution (Xb4d)
for mol,spw in (('HCN','spw27'),('HCO+','spw25')):
    c,v,ps,wc,x0,y0,b0=load_cube('Xb4d',spw,F0[mol][3],half=24.0)
    m,_,_,_,_,_,_=load_cube('Xb4d',spw,F0[mol][3],half=24.0,kind='mask')
    c=c-np.median(c[(v>-900)&(v<-450)],axis=0); pxy=pxy_fn(wc,x0,y0); n=c.shape[1]
    for nm,pos,(a,b) in (('Slobe',offset_pos(S0,0.15,-1.85),(-405,-150)),('Red',offset_pos(S0,0.15,3.3),(180,320))):
        ap_=apmask(n,ps,pxy(*pos),1.5); k=(v>=a)&(v<=b); dv=abs(np.median(np.diff(v)))
        tot=c[k][:,ap_].sum()*dv; inm=(c[k]*(m[k]>0.5))[:,ap_].sum()*dv; fr=(m[k][:,ap_]>0.5).mean()
        P(f'MASK Xb4d native {mol} {nm} r=1.5" {a}..{b}: wing sum {tot/ (1.1331*b0[0]*b0[1]/ps**2):.3f} Jy km/s total, fraction of the signal inside the clean mask {inm/tot if tot else float("nan"):.2f}; fraction of aperture-channel pixels in mask {fr:.2f}')
        J[f'mask|{mol}|{nm}']=dict(frac_signal=float(inm/tot) if tot else None,frac_pix=float(fr))
# ---- (2) luminosity and mass
D=40.0; Z=VSYS/C; DL=D                    # assumption: D_L = 40 Mpc (1" = 194 pc); z from 2780 km/s
ASSUME=dict(D_L_Mpc=D,z=Z,alpha_HCN10=10.0,r31=0.35,r31_range=(0.25,0.6),X_HCN=(1e-9,1e-8),X_HCO=(1e-9,1e-8),Tex=(10,50),mu=2.72)
rest={'HCN':265.88618,'HCO+':267.55763}
def Lp(S,mol): nu=rest[mol]/(1+Z); return 3.25e7*S*nu**-2*DL**2*(1+Z)**-3     # Solomon+1997 K km/s pc^2
# Einstein A and rotational constants (LAMDA/CDMS; B: HCN 44.316 GHz, HCO+ 44.594 GHz)
AUL={'HCN':8.36e-4,'HCO+':1.45e-3}; BROT={'HCN':44.316e9,'HCO+':44.5944e9}
h=6.62607e-34; k=1.380649e-23; c_=2.99792458e8; Mpc=3.0857e22; mH=1.6735e-27; Msun=1.989e30
def Nmol(S,mol,T):                 # optically thin LTE total number of molecules, S in Jy km/s
    nu=rest[mol]*1e9/(1+Z); F=S*1e-26*(nu/c_)*1e3             # W m^-2 : S[Jy km/s]*1e3 m/s -> Jy m/s ; F=int S dnu = S dv nu/c
    Nu=4*np.pi*(DL*Mpc)**2*F/(AUL[mol]*h*(rest[mol]*1e9))      # emitted at rest frequency
    Q=k*T/(h*BROT[mol])+1/3; Eu=BROT[mol]*h/k*12; gu=7
    return Nu*Q/gu*np.exp(Eu/T)
def mass_thin(S,mol,T,X): return Nmol(S,mol,T)/X*ASSUME['mu']*mH/Msun
cases={}
FL={ # Jy km/s in window -405..-150 (blue) / 180..320 (red); [low, nominal, high] low = r=1.5" aperture, nominal = r=1.5" aperture x recovery-free, high = 2D Gaussian integral
 'blue|HCN':(J['ap|Xb4d|HCN|blue|Slobe_r1.5']['F'],J['fit|Xb4d|HCN|blue']['flux']),
 'blue|HCO+':(J['ap|Xb4d|HCO+|blue|Slobe_r1.5']['F'],J['fit|Xb4d|HCO+|blue']['flux']),
 'redC|HCN':(J['ap|Xb4d|HCN|redC|Red']['F'],J['fit|Xb4d|HCN|redC']['flux']),
 'redC|HCO+':(J['ap|Xb4d|HCO+|redC|Red']['F'],J['fit|Xb4d|HCO+|redC']['flux'])}
v_out={'blue':abs(J['meanv|HCN']),'redC':250.0}; R_pc={'blue':1.88*194.0,'redC':3.4*194.0}
res={}
for key,(lo,hi) in FL.items():
    wing,mol=key.split('|'); Llo,Lhi=Lp(lo,mol),Lp(hi,mol)
    P(f'LUM {key}: S = {lo:.3f} (r~1-1.5" aperture) .. {hi:.3f} (2D Gaussian integral) Jy km/s -> L\'_{mol}(3-2) = {Llo:.2e} .. {Lhi:.2e} K km/s pc^2')
    # method A (HCN only; HCN-based alpha): M = alpha * L'(3-2)/r31
    if mol=='HCN':
        for r31 in (0.6,0.35,0.25):
            P(f'  MASS-A {key}: alpha(HCN 1-0)=10 Msun/(K km/s pc^2), r31={r31}: M_dense = {10*Llo/r31:.2e} .. {10*Lhi/r31:.2e} Msun')
        res[key+'|A']=dict(lo=10*Llo/0.35,hi=10*Lhi/0.35,lo_r31=(10*Llo/0.6,10*Lhi/0.25))
    # method B: optically thin LTE with abundance X
    for T in (10,20,50):
        for X in (1e-8,1e-9):
            P(f'  MASS-B {key}: thin LTE Tex={T} K, X({mol})={X:.0e}: M_gas = {mass_thin(lo,mol,T,X):.2e} .. {mass_thin(hi,mol,T,X):.2e} Msun  (N_mol={Nmol(lo,mol,T):.2e}..{Nmol(hi,mol,T):.2e})')
    res[key+'|B']={f'T{T}_X{X:.0e}':(mass_thin(lo,mol,T,X),mass_thin(hi,mol,T,X)) for T in (10,20,50) for X in (1e-8,1e-9)}
# mass-outflow rate (projected; time-average: Mdot = M v / R), nominal case = HCN method A r31=0.35 alpha=10
for wing in ('blue','redC'):
    M=res[f'{wing}|HCN|A']; v=v_out[wing]; R=R_pc[wing]; t=R*3.0857e13/v/3.156e7/1e6    # Myr (1 pc = 3.0857e13 km)
    for lab,MM in (('low',M['lo']),('high',M['hi'])):
        Md=MM*v/(R*3.0857e13/3.156e7)             # Msun/yr
        P(f'MDOT {wing} ({lab}): M_dense={MM:.2e} Msun, v={v:.0f} km/s (projected), R={R:.0f} pc (projected) -> t_dyn={t:.2f} Myr, Mdot_dense = {Md:.2f} Msun/yr (x tan i; = {Md/49:.3f} of SFR 49 Msun/yr [Michiyama+2020])')
        res[f'mdot|{wing}|{lab}']=Md
    # also lower M for r31 range
    P(f'   {wing}: r31 range (0.6..0.25) -> M_dense {M["lo_r31"][0]:.2e} .. {M["lo_r31"][1]:.2e} Msun')
json.dump(dict(res=res,assume={k:v for k,v in ASSUME.items()}),open(OUT+'/s75c.json','w'),indent=1,default=float)
open(OUT+'/s75c.txt','w').write('\n'.join(txt+['ASSUMPTIONS: '+json.dumps(ASSUME)]))
