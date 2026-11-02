# Q17 step c: (i) 3-2/1-0 brightness-temperature ratio of OF-Red (Michiyama+2018 Table 3: HCN 0.08+-0.02, HCO+ 0.08+-0.02 Jy km/s, [240,360], 'inner 2.1"') from our 3-2 flux in the same window at the M18 OF-Red position,
# compared with the same ratio at South-OF-Blue (s77d: HCN 0.38+-0.12, HCO+ 0.19+-0.06 for window -360..-180 at the lobe);
# (ii) ridge mass / Mdot / momentum bookkeeping scaled from the lobe by the matched-window flux ratio (same RADEX scaling as Q2r2: dense mass in lobe window -360..-180 = 2.1e6 (HCO+ X=1e-8) .. 7.0e6 (X=3e-9) Msun, s77e/s75c; Q2r2 report).
# T_b ratio: (S32/S10) * (nu10/nu32)^2 for equal beam and aperture (Rayleigh-Jeans, Jy km/s per same area).  1-0 data from the paper's inner-2.1" aperture; we try r=1.05", 1.5", 2.1" 3-2 apertures (1.05" beam).
# usage: python3 -I s104c_q17_r31_mass.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
OUT='/workspace/out/q17'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
def hms(h,m,s): return (h+m/60+s/3600)*15
OFR=(hms(10,27,51.20),-(43+54/60+15.64/3600)); SOB=(hms(10,27,51.17),-(43+54/60+21.43/3600))
nu32={'HCN':354.5055*0+265.8864,'HCO+':267.5576}; nu10={'HCN':88.6316,'HCO+':89.1885}
M18={'HCN':{'OFR':(0.08,0.02),'SOB':(0.07,0.02)},'HCO+':{'OFR':(0.08,0.02),'SOB':(0.10,0.02)}}
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
    fm=field_pixels(R,rmin=3.5,rmax=22.0)
    def Wm(lo,hi):
        k=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
    for tag,pos,(lo,hi) in (('OF-Red',OFR,(240,360)),('South-OF-Blue',SOB,(-360,-180))):
        M=Wm(lo,hi); key='OFR' if tag=='OF-Red' else 'SOB'
        for rad in (1.05,1.5,2.1):
            h=(np.hypot(*np.mgrid[-int(rad/ps+1):int(rad/ps+1)+1,-int(rad/ps+1):int(rad/ps+1)+1])*ps<=rad).astype(float)/R['ba']
            A=fftconvolve(M,h,mode='same'); sg=float(np.std(A[fm])); x,y=R['pxy'](*pos); f=float(M[apmask(n,ps,(x,y),rad)].sum()/R['ba'])
            s10,e10=M18[mol][key]; r=(f/s10)*(nu10[mol]/nu32[mol])**2; er=r*np.hypot(sg/f,e10/s10)
            P(f'{mol} {tag} [{lo},{hi}] r={rad}": 3-2 flux {f:.3f}+-{sg:.3f} Jy km/s; M18 1-0 {s10}+-{e10}; Tb ratio r31={r:.2f}+-{er:.2f}')
            res[f'{mol}|{tag}|{rad}']=dict(f32=f,e32=sg,r31=float(r),e_r31=float(er))
# ---- mass bookkeeping
Mw=(2.1e6,7.0e6)                         # lobe dense mass in -360..-180 (Q2r2 report), HCO+ X=1e-8 .. X=3e-9
J=json.load(open(OUT+'/s104b.json'))
for mol,key in (('HCN','HCN|180-320'),('HCO+','HCO+|180-320')):
    q=J[key]; ratio=q['ridge']/q['lobe']; er=ratio*np.hypot(q['ridge_e']/q['ridge'],q['lobe_e']/q['lobe'])
    P(f'{mol}: matched-window (|v|=180..320) ridge/lobe flux ratio {ratio:.2f}+-{er:.2f}')
    res[f'ratio|{mol}']=(ratio,er)
rH=res['ratio|HCN']; rO=res['ratio|HCO+']; w=(1/rH[1]**2,1/rO[1]**2); rm=(rH[0]*w[0]+rO[0]*w[1])/sum(w); em=1/np.sqrt(sum(w))
P(f'weighted mean ridge/lobe (HCN,HCO+; treats the two as independent, flux scale systematic +-8% not included) = {rm:.2f}+-{em:.2f}')
for lab,m in (('X(HCO+)=1e-8',Mw[0]),('X(HCO+)=3e-9',Mw[1])):
    P(f'ridge dense mass in |v|=180..360 window ({lab}): {m*rm:.2e} Msun (lobe window {m:.1e} x {rm:.2f}); bipolar total in the window = {m*(1+rm):.2e}; x1.9 if both counted (S-lobe+ridge)')
# Mdot = M v_los tan i / R_proj (Q2r2 convention), v_los = flux-weighted |<v>| in 150..400 (s104b), R_proj lobe 1.88" peak (s77g) -> 365 pc, ridge 3.3" -> 640 pc
D=194.0; Rl=1.88*D; Rr=3.3*D; vl=253.0; vr=257.0
for i_ in (70,80,85):
    t=np.tan(np.deg2rad(i_)); f=1.0227e-6   # km/s -> pc/yr
    for lab,m in (('X=1e-8',Mw[0]),('X=3e-9',Mw[1])):
        Ml=m; Mr_=m*rm; Mdl=Ml*vl*f*t/Rl; Mdr=Mr_*vr*f*t/Rr
        P(f'i={i_}: {lab}: Mdot lobe {Mdl:.1f} Msun/yr, ridge {Mdr:.1f} (ratio {Mdr/Mdl:.2f}); Pdot ~ M v^2 tan^2... momentum (los) p_lobe={Ml*vl:.2e} p_ridge={Mr_*vr:.2e} Msun km/s (ratio {Mr_*vr/(Ml*vl):.2f}); t_dyn(proj/(v tan i)) lobe {Rl/(vl*f*t)/1e6:.2f} ridge {Rr/(vr*f*t)/1e6:.2f} Myr')
json.dump(res,open(OUT+'/s104c.json','w'),indent=1,default=float); open(OUT+'/s104c.txt','w').write('\n'.join(txt))
