# Q13b: HDO 2(2,0)-3(1,3) 266.1611 GHz (Eu=157 K) is the only catalogue line (HDO, Eu<1000 K, CDMS/JPL dump of Q0) that can blend the HCN 3-2 blue wing; it is
# 310.0 km/s blue of HCN in the HCN frame.  Two tests that do not need a second HDO line:
#  (1) spectral: HCN residual R(v)=HCN-k*HCO+ on the blue side (2100-2575) fitted by (a) zero, (b) an HDO template = HCO+ profile at the same position
#      (gas kinematics, shifted by -310 km/s: v_app = -310.0 + 1.00104*v_gas), (c) a free Gaussian (3 par.), (d) both (b)+(c).
#  (2) spatial/chemical: E_blue (HCN-k*HCO+, 2375-2575) per unit HCO+ emission at the velocities that would host the HDO parent gas (HCO+ window 2685-2885)
#      at Soutflow vs S, N, Bridge, W clump: the HDO/HCO+ ratio a blend would need at the lobe relative to the limits elsewhere.
# Outputs: out/q13/q13b.log, q13b_template_fit.csv, q13b_hdo_ratio.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
from scipy.optimize import least_squares
C=299792.458; DV_HDO=(1-266.1611/265.886434)*C   # -310.0
log=open(O+'q13b.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
pr(f'HDO 266.1611 GHz offset from HCN 3-2 (265.886434): {DV_HDO:.1f} km/s; v_app = {DV_HDO:.1f} + {266.1611/265.886434:.5f}*v_gas')
rob=lambda a: 1.4826*np.median(np.abs(a-np.median(a)))
def binspec(v,s,w):
    n=int(w//abs(np.median(np.diff(v)))); k=len(v)//n; return v[:k*n].reshape(k,n).mean(1),s[:k*n].reshape(k,n).mean(1)
def run(ous,rad,pn='Soutflow',w=22.):
    Mx=wmaps(ous); dx,dy=Mx['dx'],Mx['dy']; x,y=POS[pn]; m=np.hypot(dx-x,dy-y)<=rad
    h=Mx['HCN_spec'][:,m].mean(1); o=Mx['HCO+_spec'][:,m].mean(1); vh=Mx['v_HCN']; vo=Mx['v_HCO+']
    # common velocity grid = HCN channels; interpolate HCO+ onto it
    oi=np.interp(vh,vo[np.argsort(vo)],o[np.argsort(vo)])
    # k from core 2575-2875
    mc=(vh>=2575)&(vh<2875); k=np.sum(h[mc]*oi[mc])/np.sum(oi[mc]**2)
    # empirical noise of w-wide bins on line-free HCN channels (per aperture-average)
    fm=np.zeros(len(vh),bool)
    for a,b in FREE['HCN']: fm|=(vh>=a)&(vh<=b)
    res=h-k*oi
    vb,rb=binspec(vh,res,w); _,hb=binspec(vh,h,w); _,ob=binspec(vh,oi,w)
    mfree=(vb>1380)&(vb<2100); sig=rob(rb[mfree])
    # templates
    tmpl=lambda v: np.interp((v-DV_HDO)/1.00104,vh[np.argsort(vh)],oi[np.argsort(vh)])   # HCO+ gas profile placed at HDO apparent velocity
    sel=(vb>2100)&(vb<2575); V=vb[sel]; Rr=rb[sel]
    T=tmpl(V); T0=T/T.max()
    # models
    chi0=np.sum((Rr/sig)**2)
    a=np.sum(Rr*T0)/np.sum(T0**2); chiH=np.sum(((Rr-a*T0)/sig)**2); ea=sig/np.sqrt(np.sum(T0**2))
    g=lambda p,v: p[0]*np.exp(-0.5*((v-p[1])/p[2])**2)
    f=least_squares(lambda p: (Rr-g(p,V))/sig,[0.8*Rr.max(),2500,90],bounds=([0,2150,30],[10,2650,250])); chiG=np.sum(f.fun**2)
    f2=least_squares(lambda p: (Rr-g(p[:3],V)-p[3]*T0)/sig,[0.8*Rr.max(),2500,90,0.0],bounds=([0,2150,30,-10],[10,2650,250,10])); chiGH=np.sum(f2.fun**2)
    J=f2.jac; cov=np.linalg.inv(J.T@J); eH=np.sqrt(cov[3,3])
    wm=(V>=2375)&(V<2575); sumT=T0[wm].sum(); sumR=Rr[wm].sum()
    fH=f2.x[3]*sumT/sumR; efH=eH*sumT/sumR     # HDO share of the residual flux in 2375-2575 in the joint (Gaussian wing + HDO template) fit
    # HDO template using only the bright narrow component of HCO+ is already part of oi; also a pure Gaussian-gas HDO (sigma 100, v=2696 broad comp) variant
    Tb=np.exp(-0.5*((V-DV_HDO-1.00104*2696)/(1.00104*104))**2); ab=np.sum(Rr*Tb)/np.sum(Tb**2); chiB=np.sum(((Rr-ab*Tb)/sig)**2)
    Tn=np.exp(-0.5*((V-DV_HDO-1.00104*2802)/(1.00104*32))**2); an=np.sum(Rr*Tn)/np.sum(Tn**2); chiN=np.sum(((Rr-an*Tn)/sig)**2)
    r=dict(ous=ous,rad=rad,pos=pn,k=k,bin_kms=w,sig_bin_mJy=sig*1e3,nbins=len(V),chi2_null=chi0,chi2_HDOtemplate=chiH,amp_HDOtemplate_mJy=a*1e3,e_amp=ea*1e3,
           chi2_gauss=chiG,gauss_amp_mJy=f.x[0]*1e3,gauss_v=f.x[1],gauss_sig=f.x[2],chi2_gauss_plus_HDO=chiGH,HDO_amp_in_joint_mJy=f2.x[3]*1e3,
           chi2_HDO_broadGauss=chiB,chi2_HDO_narrowGauss=chiN,amp_narrow_mJy=an*1e3,
           template_peak_v=V[np.argmax(T)],HDO_share_joint=fH,e_HDO_share_joint=efH,HDO_share_alone=a*sumT/sumR)
    return r
rows=[run(o,rad) for o in ('Xb4d','Xb4b') for rad in (0.5,1.0)]
Tt=pd.DataFrame(rows); Tt.round(3).to_csv(O+'q13b_template_fit.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
pr(Tt.round(2).T.to_string())
# ---- (2) spatial/chemical argument: E_blue per HCO+ parent-gas window ----
Mx=wmaps('Xb4d'); dx,dy=Mx['dx'],Mx['dy']; vo=Mx['v_HCO+']; vh=Mx['v_HCN']; dv=Mx['dv']
oi=Mx['HCO+_spec']; H=Mx['HCN_spec']
par=(vo>=2685)&(vo<2885); PH=oi[par].sum(0)*dv       # HCO+ in the parent-gas window (HDO of gas at 2685-2885 appears at 2375-2575 in the HCN frame)
Hh=Mx['HCN_blue']; Oo=Mx['HCO+_blue']
Ex=lambda k: Hh-k*Oo
k=0.574
rows=[]
rng=np.random.default_rng(3); ex=[POS['N'],POS['S'],POS['Bridge'],POS['Soutflow'],(-1.8,-2.4)]
for rad in (1.0,):
    nul=[]
    while len(nul)<300:
        x,y=rng.uniform(-7.5,7.5,2)
        if np.hypot(x,y)>7.5 or min(np.hypot(x-a,y-b) for a,b in ex)<3.0: continue
        if apmean(Mx['HCO+_core'],dx,dy,x,y,rad)>0.25: continue
        nul.append((x,y))
    sE=rob(np.array([apmean(Ex(k),dx,dy,x,y,rad) for x,y in nul]))
    for pn,(x,y) in {**{p:POS[p] for p in ('Soutflow','S','N','Bridge')},'Wclump':(-1.8,-2.4),'NEclump':(6.75,3.9),'Eclump':(5.25,-2.85),'SWarm':(-3.75,-4.5)}.items():
        E=apmean(Ex(k),dx,dy,x,y,rad); P=apmean(PH,dx,dy,x,y,rad)
        k_loc=apmean(Mx['HCN_core'],dx,dy,x,y,rad)/apmean(Mx['HCO+_core'],dx,dy,x,y,rad)
        E2=apmean(Hh-k_loc*Oo,dx,dy,x,y,rad)
        rows.append(dict(pos=pn,E_blue=E,sigE=sE,snr=E/sE,E_lim3=max(E,0)+3*sE,HCO_parent=P,ratio=E/P,ratio_lim3=(max(E,0)+3*sE)/P,k_local=k_loc,E_blue_klocal=E2,snr_klocal=E2/sE))
U=pd.DataFrame(rows); U.round(4).to_csv(O+'q13b_hdo_ratio.csv',index=False)
pr('E_blue (HCN-0.574*HCO+, 2375-2575) per HCO+ in the parent-gas window 2685-2885, r=1.0" at 1.05", Xb4d:'); pr(U.round(3).to_string(index=False))
lo=U[U.pos=='Soutflow'].iloc[0]
for pn in ('S','N','Bridge','Wclump'):
    r=U[U.pos==pn].iloc[0]; pr(f'  lobe/{pn}: E/HCO+parent {lo.ratio:.3f} vs {pn} <{r.ratio_lim3:.3f} (3 sigma) -> a blend with a universal HDO/HCO+ is excluded unless the lobe has >= {lo.ratio/r.ratio_lim3:.1f}x the HDO per HCO+ of {pn}')
