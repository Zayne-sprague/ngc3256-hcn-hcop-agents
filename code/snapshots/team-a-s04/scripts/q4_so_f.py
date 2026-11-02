# Q4b-6: SO temperature / excitation limits as a function of the Xb4f (Band 6, 5.7x3.5") flux recovery f of the Band-5-beam (9.4x6.3") SO flux, and of the cross-band calibration.
# Inputs: out/q4/so_matched.csv (SO 4_5-3_4 detection Xb55; SO 6_5-5_4 and 34SO 6_6-5_5 limits Xb4f; SO 8_8-7_7 limit Xb47), out/q4/so_radex_grid.csv (RADEX so@lique grid, eta=1 predictions).
# (1) LTE-thin Trot: W65 (true) = W65_obs/f ; marginalised limit by Monte Carlo over W45 (stat + 14% cross-band calibration), W65 noise, and f (fixed 0.5/0.7/1, or a flat prior in f in [0.5,1], or [0.35,1] stress test); T_max per draw from R(T)=W65/W45.
# (2) 8_8-7_7 with the same f (a Band-7 recovery f7 is taken equal to f, then independent in [0.35,1]); 34SO 6_6-5_5 (isotopologue ratio 32S/34S=22 assumed): reported as the equivalent SO 6_5 limit.
# (3) RADEX (Tkin, n) allowed region (chi2 with eta free in [0,1], f-scaled limits, 14%/10% calibration in quadrature), per f.   Output out/q4b/so_f_trot.csv, so_f_radex.csv
import numpy as np, pandas as pd
O='/workspace/out/q4b/'; S=pd.read_csv('/workspace/out/q4/so_matched.csv').set_index('line')
H,KB,CL=6.62607e-27,1.380649e-16,2.99792458e10
LINES={'SO 4_5-3_4':dict(f=178.605403,Eu=24.42876,logA=-4.15357,J=5),'SO 6_5-5_4':dict(f=251.825770,Eu=50.66106,logA=-3.71555,J=5),'SO 8_8-7_7':dict(f=344.310612,Eu=87.48155,logA=-3.28519,J=8),'34SO 6(6)-5(5)':dict(f=253.208020,Eu=55.69,logA=-3.695,J=6)}
def Wp(key,T): L=LINES[key]; return (2*L['J']+1)*np.exp(-L['Eu']/T)*10**L['logA']/L['f']**2
Tg=np.logspace(np.log10(3),np.log10(300),400)
def Tmax(ratio,key):
    """T at which the thin-LTE ratio W(key)/W(4_5) equals `ratio` (monotonic in T); 0 if ratio<=R(Tmin); inf if above R(Tmax)"""
    R=Wp(key,Tg)/Wp('SO 4_5-3_4',Tg); return np.where(ratio<=R[0],Tg[0],np.where(ratio>=R[-1],Tg[-1],np.interp(ratio,R,Tg)))
rng=np.random.default_rng(5); N=200000
W45=S.loc['SO 4_5-3_4','I_Kkms']; e45=S.loc['SO 4_5-3_4','eI_Kkms']
rows=[]
for key in ('SO 6_5-5_4','SO 8_8-7_7'):
    W,e=S.loc[key,'I_Kkms'],S.loc[key,'eI_Kkms']
    for w45lab,w45 in (('W45=0.63 (HCO+-width Gaussian)',W45),('W45=0.79 (FWHM 186 km/s)',S.loc['SO 4_5-3_4','I_alt_186'])):
        for flab,fd in (('f=1',lambda n:np.ones(n)),('f=0.7',lambda n:np.full(n,0.7)),('f=0.5',lambda n:np.full(n,0.5)),('f~U[0.5,1] flat prior',lambda n:rng.uniform(0.5,1,n)),('f~U[0.35,1] stress',lambda n:rng.uniform(0.35,1,n))):
            for cal in (0.0,0.14):
                f=fd(N); calb=np.exp(rng.normal(0,cal,N)) if cal>0 else np.ones(N)
                w45d=np.clip(rng.normal(w45,e45,N),0.05,None)
                w=rng.normal(W,e,N)/f/calb          # true W(key) in the Band-5-beam flux scale
                T=Tmax(w/w45d,key)
                rows.append(dict(line=key,W45=w45lab,f=flab,cross_band_cal=cal,Trot_max_68=np.percentile(T,68.3),Trot_max_95=np.percentile(T,95),Trot_max_997=np.percentile(T,99.7),frac_T_gt_30=np.mean(T>30)))
R=pd.DataFrame(rows); R.round(2).to_csv(O+'so_f_trot.csv',index=False); pd.set_option('display.width',250,'display.max_rows',200)
print(R[R.W45.str.startswith('W45=0.63')].round(1).to_string())
# 34SO: equivalent SO 6_5-5_4 limit at 32S/34S=22 (same f)
e34=S.loc['34SO 6(6)-5(5)','I_3sig_upper_Kkms']*22; print('34SO 3sigma limit x22 = %.2f K km/s equivalent to a SO 6_5 limit (vs direct SO 6_5 3sigma %.3f): not constraining; 34SO/SO(4_5) ratio limit = %.3f vs LTE-thin expectation 1/22=%.3f x Wratio'%(e34,S.loc['SO 6_5-5_4','I_3sig_upper_Kkms'],S.loc['34SO 6(6)-5(5)','I_3sig_upper_Kkms']/W45,1/22))
# ---- RADEX allowed (Tkin, n) per f
G=pd.read_csv('/workspace/out/q4/so_radex_grid.csv')
e=np.array([np.hypot(e45,0.10*W45),0,0]); out=[]
for f in (1.0,0.7,0.5):
    y=np.array([W45,S.loc['SO 6_5-5_4','I_Kkms']/f,S.loc['SO 8_8-7_7','I_Kkms']/f])
    e=np.array([np.hypot(e45,0.10*W45),np.hypot(S.loc['SO 6_5-5_4','eI_Kkms']/f,0.14*0.047),np.hypot(S.loc['SO 8_8-7_7','eI_Kkms']/f,0.14*0.1)])
    P=G[['pred45','pred65','pred88']].values; eta=np.clip((P*y/e**2).sum(1)/np.maximum((P*P/e**2).sum(1),1e-30),0,1); chi=(((y-eta[:,None]*P)/e)**2).sum(1)
    G['chi_f']=chi; G['eta_f']=eta; m=chi.min()
    for lab,dc in (('2sig(3par)',8.0),):
        s_=G[G.chi_f<=m+dc]
        for T,gg in s_.groupby('Tkin'):
            out.append(dict(f=f,Tkin=T,level=lab,logn_min=gg.log_n.min(),logn_max=gg.log_n.max(),logN_min=gg.log_N.min(),logN_max=gg.log_N.max(),eta_min=gg.eta_f.min(),eta_max=gg.eta_f.max(),chi2min_all=m))
        # highest allowed Tkin at each density
        for ln in (3,4,5,6,7,8):
            q=s_[np.isclose(s_.log_n,ln)]; out.append(dict(f=f,Tkin=-1,level=f'Tkin_max at log n={ln}',logn_min=ln,logn_max=ln,logN_min=np.nan,logN_max=np.nan,eta_min=q.Tkin.min() if len(q) else np.nan,eta_max=q.Tkin.max() if len(q) else np.nan,chi2min_all=m))
RR=pd.DataFrame(out); RR.round(3).to_csv(O+'so_f_radex.csv',index=False); print(RR[RR.Tkin==-1][['f','level','eta_min','eta_max']].rename(columns={'eta_min':'Tkin_min_allowed','eta_max':'Tkin_max_allowed'}).to_string()); print(RR[RR.Tkin>0].groupby('f').apply(lambda g:g.assign(r=g.Tkin.astype(str)+': logn '+g.logn_min.astype(str)+'-'+g.logn_max.astype(str)).r.tolist()).to_string())
