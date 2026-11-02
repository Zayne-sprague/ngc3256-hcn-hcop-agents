# Q1 fix 3+6: HCN/HCO+ 3-2 (Tb) as a function of velocity at N, S, Bridge, Wclump, Soutflow (Xb4d convolved to a common 1.05" circular beam),
# window-width dependence (+-200/+-350/+-500, HCN v2=1f guard), line profile comparison (centroid, FWHM, dips), and the Soutflow wing ratio:
# window variants, Xb4b convolved to 1.05", two-Gaussian (core+blue wing) spectral fits and a residual bootstrap.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
from matched import smoothed_spectrum
from scipy.optimize import curve_fit
import pandas as pd
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15),'Wclump':(-1.8,-2.4),'Soutflow':(-0.15,-4.5)}
TB=(1.05,1.05,0.)
FH,FO=265.886434,267.557626
fH2=267.1993                                   # HCN v2=1f J=3-2 (l-type doublet, 267.199 GHz); position in the HCO+ velocity frame:
V2=(1-fH2/FO)*C; print('HCN v2=1f in the HCO+ radio-velocity frame: %+.0f km/s relative to HCO+ (positive = redshifted)'%V2)
cT=Tb_factor(FH/(1+VSYS/C),1.05)/Tb_factor(FO/(1+VSYS/C),1.05)   # HCN Jy/beam -> Tb relative to HCO+ Jy/beam (same beam)
def spectra(ous,spwH,spwO,pos):
    out={}
    for nm,(ra_off,de_off) in pos.items():
        ra=RA0+ra_off/3600/np.cos(np.deg2rad(DE0)); de=DE0+de_off/3600; S={}
        for sp,spw,f0 in (('HCN',spwH,FH),('HCO+',spwO,FO)):
            p=f(ous,spw,kind='cube.I.pbcor')
            sm,nat,nb=smoothed_spectrum(p,TB,ra,de); hd=fits.open(p,memmap=True)[0].header
            v=(1-freqaxis(hd)/1e9/f0)*C; S[sp]=(v,sm)
        out[nm]=S
    return out
def noise(v,s,vc=VSYS,excl=700):
    m=np.abs(v-vc)>excl; x=s[m]; return 1.4826*np.median(np.abs(x-np.median(x)))
SP=spectra('Xb4d','spw27','spw25',POS)
rows=[];bins=[];prof=[]
dv=4.4
for nm,S in SP.items():
    (vH,sH),(vO,sO)=S['HCN'],S['HCO+']; eH,eO=noise(vH,sH),noise(vO,sO)
    # ---- window dependence ----
    for lab,(lo,hi,cut) in {'+-200':(-200,200,None),'+-350':(-350,350,None),'+-500 (raw)':(-500,500,None),'+-500 minus v2=1f channels':(-500,500,(V2-45,V2+45)),
                            '-350..+350 minus v2=1f':(-350,350,(V2-45,V2+45))}.items():
        def integ(v,s,e):
            m=(v-VSYS>=lo)&(v-VSYS<=hi)
            if cut: m&=~((v-VSYS>=cut[0])&(v-VSYS<=cut[1]))
            return s[m].sum()*dv,e*dv*np.sqrt(m.sum()),m.sum()
        a,ea,n=integ(vH,sH,eH); b,eb,_=integ(vO,sO,eO); R=a/b*cT
        rows.append(dict(region=nm,window=lab,I_HCN=a,e_HCN=ea,I_HCO=b,e_HCO=eb,R=R,eR=abs(R)*np.hypot(ea/a,eb/b),snr_HCN=a/ea,snr_HCO=b/eb,nch=n))
    # ---- 50 km/s bins ----
    for lo in np.arange(-500,500,50):
        m=lambda v:(v-VSYS>=lo)&(v-VSYS<lo+50)
        a=sH[m(vH)].sum()*dv; b=sO[m(vO)].sum()*dv; nn=m(vH).sum(); ea=eH*dv*np.sqrt(nn); eb=eO*dv*np.sqrt(m(vO).sum())
        row=dict(region=nm,v_lo=lo,v_hi=lo+50,I_HCN=a,e_HCN=ea,I_HCO=b,e_HCO=eb,snr_HCN=a/ea,snr_HCO=b/eb)
        if a/ea>3 and b/eb>3: R=a/b*cT; row.update(R=R,eR=R*np.hypot(ea/a,eb/b),kind='value')
        elif b/eb>3: row.update(R=3*ea/b*cT,kind='upper limit')
        elif a/ea>3: row.update(R=a/(3*eb)*cT,kind='lower limit')
        else: row.update(R=np.nan,kind='undetermined')
        bins.append(row)
    # ---- profile comparison ----
    for sp,(v,s,e) in (('HCN',(vH,sH,eH)),('HCO+',(vO,sO,eO))):
        m=np.abs(v-VSYS)<500; vv=v[m]; ss=s[m]; pk=ss.max(); vp=vv[ss.argmax()]
        mm=(vv-VSYS>-300)&(vv-VSYS<300); cen=(vv[mm]*ss[mm]).sum()/ss[mm].sum()
        above=np.where(ss>pk/2)[0]; fw=vv[above[-1]]-vv[above[0]]
        # deepest negative dip inside the line (+-300 km/s) in units of the channel noise (absorption / self-absorption would show as dips relative to the neighbouring HCN/HCO+ shape)
        prof.append(dict(region=nm,line=sp,peak_mJy=pk*1e3,v_peak_minus_vsys=vp-VSYS,centroid_minus_vsys=cen-VSYS,fwhm_halfmax_span=fw,chan_noise_mJy=e*1e3,min_in_line_over_noise=ss[mm].min()/e))
    # shape difference: HCN scaled (to HCO+ at the peak channel, Tb) minus HCO+, per 50 km/s, in sigma (dips of one line relative to the other)
    pkO=sO[np.abs(vO-VSYS)<300].max(); pkH=sH[np.abs(vH-VSYS)<300].max(); sc=pkO/pkH/cT
    # (kept in bin table via ratio)
B=pd.DataFrame(bins); B.to_csv('/workspace/out/q1/q1_ratio_vs_velocity.csv',index=False)
Wd=pd.DataFrame(rows); Wd.to_csv('/workspace/out/q1/q1_ratio_windows.csv',index=False)
Pf=pd.DataFrame(prof); Pf.to_csv('/workspace/out/q1/q1_line_profiles.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',400)
print(Wd.round(3).to_string()); print(Pf.round(2).to_string())
print(B[B.region.isin(['N','S','Bridge'])].pivot(index='v_lo',columns='region',values='R').round(2).to_string())
# --------------------- Soutflow wing ratio ---------------------
def gfit(v,s,e,vcore0,excl=None):
    m=np.abs(v-VSYS)<700
    if excl: m&=~((v-VSYS>=excl[0])&(v-VSYS<=excl[1]))
    x,y=v[m],s[m]
    g=lambda x,a1,c1,w1,a2,c2,w2: a1*np.exp(-0.5*((x-c1)/w1)**2)+a2*np.exp(-0.5*((x-c2)/w2)**2)
    p0=[y.max(),vcore0,80,0.3*y.max(),vcore0-250,120]
    try:
        p,cv=curve_fit(g,x,y,p0=p0,sigma=np.full(len(x),e),bounds=([0,vcore0-150,20,0,vcore0-600,40],[1e3,vcore0+150,300,1e3,vcore0-100,300]),maxfev=20000)
    except Exception as ex: return None
    return x,y,p,g,m
SO=SP['Soutflow']; W=[]
for sp in ('HCN','HCO+'):
    pass
# wing windows (absolute radio velocity, as in p0001)
def wing_window(S,lo,hi,eH,eO):
    (vH,sH),(vO,sO)=S['HCN'],S['HCO+']; mH=(vH>=lo)&(vH<=hi); mO=(vO>=lo)&(vO<=hi)
    a=sH[mH].sum()*dv; b=sO[mO].sum()*dv; ea=eH*dv*np.sqrt(mH.sum()); eb=eO*dv*np.sqrt(mO.sum()); return a,ea,b,eb
(vH,sH),(vO,sO)=SO['HCN'],SO['HCO+']; eH,eO=noise(vH,sH),noise(vO,sO)
wrows=[]
for lo,hi in ((2272,2572),(2300,2500),(2200,2600),(2350,2550),(2272,2472),(2372,2572)):
    a,ea,b,eb=wing_window(SO,lo,hi,eH,eO); R=a/b*cT
    # note: the blue wing sits at |dv|=200-500 km/s; HCN v2=1f falls at %+.0f km/s in the HCO+ frame, i.e. outside blue windows
    wrows.append(dict(set='Xb4d 1.05"',window=f'{lo}-{hi}',I_HCN=a,e_HCN=ea,I_HCO=b,e_HCO=eb,snr_HCN=a/ea,snr_HCO=b/eb,R=R,eR=abs(R)*np.hypot(ea/a,eb/b)))
# Xb4b convolved to 1.05"
SB=spectra('Xb4b','spw27','spw25',{'Soutflow':POS['Soutflow'],'N':POS['N'],'S':POS['S']})
for nm in ('Soutflow','N','S'):
    S=SB[nm]; (vH2,sH2),(vO2,sO2)=S['HCN'],S['HCO+']; e1,e2=noise(vH2,sH2),noise(vO2,sO2)
    for lo,hi in ((2272,2572),(2300,2500),(2200,2600)):
        a,ea,b,eb=wing_window(S,lo,hi,e1,e2); R=a/b*cT
        wrows.append(dict(set=f'Xb4b->1.05" {nm}',window=f'{lo}-{hi}',I_HCN=a,e_HCN=ea,I_HCO=b,e_HCO=eb,snr_HCN=a/ea,snr_HCO=b/eb,R=R,eR=abs(R)*np.hypot(ea/a,eb/b)))
# same windows in the Xb4d spectrum at N and S for reference (blue wing at the nuclei)
for nm in ('N','S'):
    S=SP[nm]; (vH2,sH2),(vO2,sO2)=S['HCN'],S['HCO+']; e1,e2=noise(vH2,sH2),noise(vO2,sO2)
    a,ea,b,eb=wing_window(S,2272,2572,e1,e2); R=a/b*cT
    wrows.append(dict(set=f'Xb4d 1.05" {nm}',window='2272-2572',I_HCN=a,e_HCN=ea,I_HCO=b,e_HCO=eb,snr_HCN=a/ea,snr_HCO=b/eb,R=R,eR=abs(R)*np.hypot(ea/a,eb/b)))
WR=pd.DataFrame(wrows)
# Gaussian (core + blue wing) fits and bootstrap of the wing ratio
gfits=[];
vc0=2722.
fits_={}
for sp,v,s,e in (('HCN',vH,sH,eH),('HCO+',vO,sO,eO)):
    r=gfit(v,s,e,vc0)
    if r is None: print('gfit failed',sp); continue
    x,y,p,g,m=r; fits_[sp]=(x,y,p,g,e); res=y-g(x,*p)
    wing=p[3]*p[5]*np.sqrt(2*np.pi)/dv*dv   # Jy/beam km/s  (a2*w2*sqrt(2pi))
    core=p[0]*p[2]*np.sqrt(2*np.pi)
    gfits.append(dict(line=sp,core_amp_mJy=p[0]*1e3,core_v=p[1],core_sigma=p[2],wing_amp_mJy=p[3]*1e3,wing_v=p[4],wing_sigma=p[5],wing_flux=wing,core_flux=core,chi2r=float(np.sum((res/e)**2)/(len(x)-6))))
rng=np.random.default_rng(3); boot=[]
if len(fits_)==2:
    for _ in range(300):
        wr={}
        for sp,(x,y,p,g,e) in fits_.items():
            res=y-g(x,*p); yb=g(x,*p)+rng.choice(res,len(res),replace=True)
            try: pb,_=curve_fit(g,x,yb,p0=p,sigma=np.full(len(x),e),bounds=([0,vc0-150,20,0,vc0-600,40],[1e3,vc0+150,300,1e3,vc0-100,300]),maxfev=20000)
            except Exception: pb=None
            wr[sp]=pb
        if all(wr[k] is not None for k in wr):
            fl={k:wr[k][3]*wr[k][5]*np.sqrt(2*np.pi) for k in wr}; cr={k:wr[k][0]*wr[k][2]*np.sqrt(2*np.pi) for k in wr}
            boot.append((fl['HCN']/fl['HCO+']*cT,cr['HCN']/cr['HCO+']*cT))
boot=np.array(boot); G=pd.DataFrame(gfits)
if len(boot):
    G.loc[len(G)]=dict(line='RATIO wing HCN/HCO+ (Tb) [boot med,p16,p84]',core_amp_mJy=np.median(boot[:,0]),core_v=np.percentile(boot[:,0],16),core_sigma=np.percentile(boot[:,0],84),
                       wing_amp_mJy=np.median(boot[:,1]),wing_v=np.percentile(boot[:,1],16),wing_sigma=np.percentile(boot[:,1],84),wing_flux=len(boot),core_flux=np.nan,chi2r=np.nan)
    G.to_csv('/workspace/out/q1/q1_soutflow_gauss.csv',index=False)
WR.to_csv('/workspace/out/q1/q1_soutflow_wing_variants.csv',index=False)
print(WR.round(3).to_string()); print(G.round(3).to_string())
# compare the position of the wing at the nuclei with a one-component null test: wing ratio at N,S in same window printed above
