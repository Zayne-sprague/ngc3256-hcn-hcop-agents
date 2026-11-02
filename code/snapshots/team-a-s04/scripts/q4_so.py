# Q4d: SO at ONE matched beam (9.4x6.3" PA -88, the Band-5 beam) at the centre of the field (RA0,DE0; N+S+Bridge unresolved in Band 5):
#   SO 5(4)-4(3) [N_J=4_5-3_4, 178.605 GHz, Eu 24 K, Xb55 spw18]  detected, fitted jointly with the HCO+ 2-1 wing and HOC+ 2-1 (exponential-wing + Gaussians, q4_lib.fit_target)
#   SO 5(6)-4(5) [6_5-5_4, 251.826 GHz, Eu 51 K, Xb4f spw22]       limit (CH3OH A-cluster blend noted)      ; SO 8(8)-7(7) [344.311 GHz, Eu 87 K, Xb47 spw20] limit, HC15N 4-3 blend modelled
#   34SO 6(6)-5(5) [253.208, Eu 56 K, Xb4f spw22] limit.
# Then: LTE (thin) rotational temperature limit from the line ratios, and RADEX (LAMDA so@lique.dat) (Tkin, n, N) allowed region.  Outputs out/q4/so_matched.csv, so_lte.csv, so_radex.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_lib import *
from scipy.optimize import least_squares
RA0,DE0=156.96346,-43.90461
TB=(9.4,6.3,-88.)
H,KB,CL=6.62607e-27,1.380649e-16,2.99792458e10
LINES={'SO 4_5-3_4':dict(f=178.605403,Eu=24.42876,logA=-4.15357,J=5),'SO 6_5-5_4':dict(f=251.825770,Eu=50.66106,logA=-3.71555,J=5),'SO 8_8-7_7':dict(f=344.310612,Eu=87.48155,logA=-3.28519,J=8),'34SO 6(6)-5(5)':dict(f=253.208020,Eu=55.69,logA=-3.695,J=6)}
def gfit(v,s):   # Gaussian fit (core of HCO+ 2-1) for centroid / FWHM
    from scipy.optimize import curve_fit
    m=np.abs(v-2780)<500; p,_=curve_fit(lambda x,a,v0,w:a*np.exp(-0.5*((x-v0)/(w/2.355))**2),v[m],s[m],p0=[s[m].max(),2780,150]); return p
rows=[]
# --- SO 4_5-3_4 in Xb55 spw18
fr,s,_,nb=spectrum('Xb55','spw18',None,TB,pos=(RA0,DE0)) if False else spectrum('Xb55','spw18','N',TB,pos=(RA0,DE0))
fh=178.375056; v=vel(fr,fh); a,v0,w=gfit(v,s*1e3)
off=(1-LINES['SO 4_5-3_4']['f']/fh)*C; offH=(1-178.972051/fh)*C
res={}
for lab,fw in (('FWHM=HCO+2-1 core',w),('FWHM=186 (blind Gaussian)',186.)):
    Iv,ev,Iw,_,nch,tau=fit_target(v,s*1e3,v0,fw,off,extras=(offH,),vmax=1300)
    res[lab]=(Iv,ev); 
Iv,ev=res['FWHM=HCO+2-1 core']; k=kb(LINES['SO 4_5-3_4']['f'],TB[0],TB[1])*1e-3
rows.append(dict(line='SO 4_5-3_4',ous='Xb55',beam=f'{TB[0]}x{TB[1]}',v0=v0,fwhm_HCO=w,I_Kkms=Iv*k,eI_Kkms=ev*k,I_alt_186=res['FWHM=186 (blind Gaussian)'][0]*k,snr=Iv/ev,rms_note='stat+tau spread'))
# HCO+ 2-1 integral at the same position/beam for ratios (Gaussian-core window +-1.5 FWHM, local baseline)
vm=(v>v0-1.5*w)&(v<v0+1.5*w); dv=abs(np.median(np.diff(v))); IH=np.sum(s[vm]*1e3)*dv*kb(fh,TB[0],TB[1])*1e-3
rows.append(dict(line='HCO+ 2-1',ous='Xb55',beam=f'{TB[0]}x{TB[1]}',v0=v0,fwhm_HCO=w,I_Kkms=IH,eI_Kkms=0.03*IH,snr=np.nan,rms_note='window +-1.5 FWHM, no baseline (bright)'))
# --- limits: SO 6_5-5_4 (Xb4f spw22) and 34SO; window +-1 FWHM (=w) about v0, linear baseline from +-(2.2..6) FWHM  [CH3OH A-cluster/other lines in the baseline excluded by the robust MAD]
def lim_line(ous,spw,key,extra_gauss=None):
    f0=LINES[key]['f']; fr,s,_,nb=spectrum(ous,spw,'N',TB,pos=(RA0,DE0)); v=vel(fr,f0); s=s*1e3
    x=v-v0; m=(np.abs(x)>2.2*w)&(np.abs(x)<1000)
    if extra_gauss:   # model blends: Gaussians at fixed offsets, same width (linear fit incl. baseline), only window channels used for the line
        sig=w/2.355; mm=np.abs(x)<1000; A=np.vstack([np.ones(mm.sum()),x[mm],np.exp(-0.5*(x[mm]/sig)**2)]+[np.exp(-0.5*((x[mm]-e)/sig)**2) for e in extra_gauss]).T
        c,*_=np.linalg.lstsq(A,s[mm],rcond=None); r=s[mm]-A@c; sd=robust_sigma(r); cov=np.linalg.inv(A.T@A)*sd**2
        dvv=abs(np.median(np.diff(v))); n=sig*np.sqrt(2*np.pi)/dvv; Iv=c[2]*n*dvv; ev=np.sqrt(cov[2,2])*n*dvv
    else:
        p=np.polyfit(x[m],s[m],1); r=s-np.polyval(p,x); sd=robust_sigma(r[m]); k_=np.abs(x)<=w; dvv=abs(np.median(np.diff(v))); Iv=r[k_].sum()*dvv; ev=sd*dvv*np.sqrt(k_.sum())
    kk=kb(f0,TB[0],TB[1])*1e-3; return Iv*kk,ev*kk,sd
for ous,spw,key,ex in (('Xb4f','spw22','SO 6_5-5_4',None),('Xb4f','spw22','34SO 6(6)-5(5)',None),('Xb47','spw20','SO 8_8-7_7',[(1-344.200320/344.310612)*C])):
    I,e,sd=lim_line(ous,spw,key,ex); rows.append(dict(line=key,ous=ous,beam=f'{TB[0]}x{TB[1]}',v0=v0,fwhm_HCO=w,I_Kkms=I,eI_Kkms=e,I_alt_186=np.nan,snr=I/e,rms_note=f'sigma_ch={sd:.2f} mJy/beam; 3sig limit {max(I,0)+3*e:.3f} K km/s'+('; HC15N 4-3 at -96 km/s fitted jointly' if ex else '')))
R=pd.DataFrame(rows); R['I_3sig_upper_Kkms']=np.where(R.line.isin(['SO 4_5-3_4','HCO+ 2-1']),np.nan,np.maximum(R.I_Kkms,0)+3*R.eI_Kkms); R.round(4).to_csv(OUT+'so_matched.csv',index=False)
pd.set_option('display.width',250); print(R.round(3).to_string())
# --- LTE thin: ratios of integrated Tb vs Trot
def Wpred(key,T): L=LINES[key]; g=2*L['J']+1; return g*np.exp(-L['Eu']/T)*10**L['logA']/L['f']**2
Tg=np.array([4,5,6,8,10,12,15,20,25,30,40,50,70,100.]); o=R.set_index('line')
obs=o.loc['SO 4_5-3_4','I_Kkms']; eobs=o.loc['SO 4_5-3_4','eI_Kkms']
ltes=[]
for key in ('SO 6_5-5_4','SO 8_8-7_7','34SO 6(6)-5(5)'):
    lim=o.loc[key,'I_3sig_upper_Kkms']
    for T in Tg:
        pred=obs*Wpred(key,T)/Wpred('SO 4_5-3_4',T) if key!='34SO 6(6)-5(5)' else np.nan
        ltes.append(dict(line=key,Trot=T,W_pred_Kkms=pred,W_3sig_limit=lim,excluded=bool(pred>lim) if np.isfinite(pred) else np.nan))
L=pd.DataFrame(ltes); L.to_csv(OUT+'so_lte.csv',index=False)
for key in ('SO 6_5-5_4','SO 8_8-7_7'):
    s_=L[(L.line==key)]; ok=s_[~s_.excluded.astype(bool)]; print(key,'LTE Trot allowed (thin, SO 4_5 detected value):',ok.Trot.min() if len(ok) else None,'-',ok.Trot.max() if len(ok) else None,'K ; excluded for Trot >=',s_[s_.excluded.astype(bool)].Trot.min() if s_.excluded.any() else None)
# SO column (LTE, Trot=T) from the 4_5-3_4 line:  N_u=8 pi k nu^2 W/(h c^3 A) ; N=N_u/g*Q*exp(Eu/T);  Q(SO, 3Sigma-) ~ 3*(kT/hB)*(1+hB/(3kT)),  B=21.5 GHz -> hB/k=1.033 K
def Nsol(T,W):
    Lq=LINES['SO 4_5-3_4']; Nu=8*np.pi*KB*(Lq['f']*1e9)**2*W*1e5/(H*CL**3*10**Lq['logA']); Q=3*(T/1.033)*(1+1.033/(3*T)); return Nu/(2*Lq['J']+1)*Q*np.exp(Lq['Eu']/T)
for T in (6,8,10,15,20,30): print('N(SO) LTE Trot',T,'= %.2e cm-2 (thin, beam-filling 1, from W=%.3f)'%(Nsol(T,obs),obs))
# --- RADEX so@lique: (Tkin, n, N) grid; predicted Tb of the 3 lines, filling 1; chi2 vs detection + limits (limits as 0 +- sigma)
import subprocess,os,re,itertools
from concurrent.futures import ProcessPoolExecutor
def rdx(args):
    T,n,N,dv=args; wd=f'/tmp/rxso_{os.getpid()}'; os.makedirs(wd,exist_ok=True)
    inp=f'/opt/lamda/so@lique.dat\n{wd}/o.txt\n150 400\n{T}\n1\nH2\n{n:.4e}\n2.73\n{N:.4e}\n{dv}\n0\n'; subprocess.run(['radex'],input=inp,text=True,capture_output=True); out={}
    for l in open(f'{wd}/o.txt').read().split('\n'):
        p=l.replace('--',' ').split()
        if len(p)==12 and p[0][0].isdigit():
            try: out[round(float(p[3]),1)]=(float(p[7]),float(p[6]))   # T_R, tau
            except Exception: pass
    return out
if __name__=='__main__':
    TKs=[8,10,15,20,30,50,80,150]; LNs=np.arange(3,8.01,0.5); LNc=np.arange(12,16.01,0.25); DVs=160.
    jobs=list(itertools.product(TKs,10**LNs,10**LNc,[DVs]))
    with ProcessPoolExecutor(4) as ex: res=list(ex.map(rdx,jobs,chunksize=40))
    keys={'SO 4_5-3_4':178.6,'SO 6_5-5_4':251.8,'SO 8_8-7_7':344.3}
    o=R.set_index('line'); y=np.array([o.loc[k,'I_Kkms'] for k in keys]); e=np.array([np.hypot(o.loc[k,'eI_Kkms'],0.1*abs(o.loc['SO 4_5-3_4','I_Kkms'])*(k=='SO 4_5-3_4')) for k in keys])
    dvk=DVs*np.sqrt(2*np.pi)/2.355   # Gaussian line of FWHM dv: W = TR * 1.064*FWHM
    rows=[]
    for (T,n,N,_),r in zip(jobs,res):
        try: pr=np.array([r[round(f0,1)][0]*1.0645*DVs for f0 in (178.6,251.8,344.3)]); ta=np.array([r[round(f0,1)][1] for f0 in (178.6,251.8,344.3)])
        except Exception: continue
        # RADEX T_R already includes the background-subtracted Tb for a uniform sphere; W = T_R * 1.0645 * dv (Gaussian profile of FWHM dv)
        eta=np.clip(np.sum(pr*y/e**2)/max(np.sum(pr*pr/e**2),1e-30),0,1.); chi=np.sum(((y-eta*pr)/e)**2)
        rows.append(dict(Tkin=T,log_n=np.log10(n),log_N=np.log10(N),eta=eta,chi2=chi,tau45=ta[0],tau65=ta[1],pred45=pr[0],pred65=pr[1],pred88=pr[2]))
    G=pd.DataFrame(rows); G.to_csv(OUT+'so_radex_grid.csv',index=False); m=G.chi2.min(); print('RADEX SO: min chi2 %.2f of 3 data, 3 params'%m)
    for lab,dc in (('1sig(3par)',3.53),('2sig',8.0)):
        s_=G[G.chi2<=m+dc]; print(lab,'Tkin',sorted(s_.Tkin.unique()),'log n %.1f-%.1f'%(s_.log_n.min(),s_.log_n.max()),'log N %.2f-%.2f'%(s_.log_N.min(),s_.log_N.max()),'eta %.2f-%.2f'%(s_.eta.min(),s_.eta.max()),'tau(4_5) max %.2f'%s_.tau45.max())
    # allowed Tkin vs density (eta free): min chi2 per (Tkin, log n)
    P=G.groupby(['Tkin','log_n']).chi2.min().unstack(); print((P-m).round(1).to_string())
    G.groupby(['Tkin','log_n']).chi2.min().unstack().to_csv(OUT+'so_radex_chi2_T_n.csv')
