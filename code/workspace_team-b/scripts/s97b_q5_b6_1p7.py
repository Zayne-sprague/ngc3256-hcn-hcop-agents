# Q5 round 2: Band 6 lines at Harada's 1.7" resolution (Xb4d smoothed, s97a): HCO+/HCN 3-2 integrated over +-150 km/s of 2775 km/s as in Harada+2018 (K km/s, linear baseline),
# and matched-filter (HCO+ 3-2 profile template) amplitudes of CH3OH 2(0)-1(-1)E 254.015 and SO 251.826 (with/without the CH3OH A-ladder), with empirical noise from sliding the template over the spw.
# Output: out/q5/harada_compare_1p7.csv (ours), out/q5/harada_vs_ours.csv.   usage: python3 -I s97b_q5_b6_1p7.py
import sys; sys.path.insert(0,'/workspace/scripts')
import s91_q5_lib as L
from s91_q5_lib import *
L.REFLINE['Xb4s']=('spw25',267.55763); L.HCNLINE['Xb4s']=('spw27',265.88618)
REG=['N','S','Bridge','Sblue','W']; O='Xb4s'; TH=1.7
SO=('SO 6_5-5_4',251.82577); CH3=('CH3OH 2_0-1_-1 E',254.01538)
d=pd.read_csv('/workspace/out/linelist.csv'); m=d[(d.sp=='CH3OH')&(d.freq_GHz>251.70)&(d.freq_GHz<251.93)].copy(); m['k']=m.freq_GHz.round(3); m=m.sort_values('sijmu2').groupby('k').tail(1)
ladder=lambda T: list(zip(m.freq_GHz.values,(m.freq_GHz**3*m.sijmu2*np.exp(-m.upper_state_energy_K/T)).values))
def K(S,nu_obs): return 1.222e6*S/(nu_obs**2*TH**2)      # Jy/beam -> K, nu_obs GHz
rows=[]
def mom(spw,reg,f0,vc=2775.,w=150.):
    nu,y=spec(O,spw,reg); v=vel(nu,f0)  # v relative to line rest
    o=np.argsort(v); v,y=v[o],y[o]/1e3; fl=(abs(v-vc)>500)&(abs(v-vc)<1100)&np.isfinite(y); pf=np.polyfit(v[fl],y[fl],1); yb=y-np.polyval(pf,v)
    s=abs(v-vc)<w; dv=abs(np.median(np.diff(v))); I=yb[s].sum()*dv; rms=np.std(yb[fl]); e=rms*dv*np.sqrt(s.sum())
    nuo=np.interp(vc,v,nu[o])/1e9; return I,e,K(I,nuo),K(e,nuo),K(yb[s].max(),nuo)*1e3
for reg in REG:
    r={'pos':reg}
    for nm,spw,f0 in (('HCO+ 3-2','spw25',267.55763),('HCN 3-2','spw27',265.88618)):
        I,e,Ik,ek,pk=mom(spw,reg,f0); r[nm+' Jy']=I; r[nm+' K']=Ik; r[nm+' eK']=ek; r[nm+' Tpk_mK']=pk
    shp,vc,fw,snr=line_shape(O,reg,mode='profile',tname='ref'); r['v0']=vc; r['fwhm']=fw; r['snr_ref']=snr
    nuo=265.1; fo=lambda f0: f0*(1-vc/C)
    # CH3OH 254 (spw29) : joint with SO2/SiS/NS neighbours as in s91c; single-line (+ baseline) too
    SO2=('SO2 6(3,3)-6(2,4)',254.2805); SIS=('SiS 14-13',254.1031); NS=('NS 11/2-9/2 (f)',253.9685)
    rr=joint_fit(O,'spw29',reg,[SO2,CH3,SIS,NS],shp,vc,window=1100,order=1,center=SO2[1]); r['CH3OH254 W']=rr[CH3[0]][0]; r['CH3OH254 sf']=rr[CH3[0]][1]
    rr1=joint_fit(O,'spw29',reg,[CH3],shp,vc,window=500,order=1); r['CH3OH254 W_single']=rr1[CH3[0]][0]; r['CH3OH254 sf_single']=rr1[CH3[0]][1]
    ms=slide_sigma(O,'spw29',reg,CH3[1],shp,vc,window=500); r['CH3OH254 mad']=ms[0]; r['CH3OH254 median']=ms[3]
    # SO
    for T in (None,10):
        if T is None: rr=joint_fit(O,'spw31',reg,[SO],shp,vc,window=700,order=1); k='SO'; w=rr[SO[0]]
        else: rr=joint_fit(O,'spw31',reg,[SO],shp,vc,window=1000,order=1,extra=[('lad',ladder(10))],center=SO[1]); k='SO+lad'; w=rr[SO[0]]
        r[k+' W']=w[0]; r[k+' sf']=w[1]
    ms=slide_sigma(O,'spw31',reg,SO[1],shp,vc,window=700); r['SO mad']=ms[0]; r['SO median']=ms[3]
    rows.append(r)
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/harada_compare_1p7.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',50); print(df.T.to_string())
