# Unified per-line / per-region / per-OUS measurement table (one convention).
# Spectrum: beam-sized aperture ('ap') and peak pixel ('pix') from out/regspec (s7a). Window/template = the same region's HCO+ profile
# (same OUS, same spectral type) in absolute radio velocity; weak lines: matched-filter amplitude I_T [Jy/beam km/s] (= integrated intensity if
# the line has the template profile) and template-window box sum W. Errors: stat (line-free rms, lag-1 inflation), empirical
# (template slid over line-free channels of the same spectrum = includes baseline structure), flux-scale and pb terms added separately.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
from regions import REGIONS
import pandas as pd
Z={}
def spec(ous,spw,reg,st):
    k=(ous,spw)
    if k not in Z: Z[k]=np.load(f'/workspace/out/regspec/{ous}.{spw}.cube.pbcor.npz')
    z=Z[k]; key=f'{reg}_{st}'
    return (z['nu'],z[key]) if key in z else (None,None)
MAP={'spw16':'spw25','spw18':'spw27','spw20':'spw29','spw22':'spw31'}
TEMPL={'Xb55':('spw18',178.37506),'Xb47':('spw18',356.73424),'Xb4f':('spw16',267.55763),'Xb4d':('spw25',267.55763),'Xb4b':('spw25',267.55763)}
# (label, rest GHz, [(ous,spw)], note)
def b6(sp): return [('Xb4f',sp),('Xb4d',MAP[sp]),('Xb4b',MAP[sp])]
LINES=[
 ('HCN 2-1',177.26111,[('Xb55','spw16')],'strong'),('HCO+ 2-1',178.37506,[('Xb55','spw18')],'strong'),
 ('HOC+ 2-1',178.97205,[('Xb55','spw18')],'tentative; no competing species with Eu<300 K in curated list'),
 ('H41beta RRL',177.72284,[('Xb55','spw16')],'near HCN v2=1e J=2-1 (177.699)'),
 ('CH3OH J-1-J0 E blend',165.06113,[('Xb55','spw22')],'1_-1-1_0/2_-1-2_0/3_-1-3_0 E at 165.050/.061/.099'),
 ('CH3CN 9-8 K=0,1',165.56825,[('Xb55','spw22')],'marginal'),
 ('HCN 4-3',354.50547,[('Xb47','spw16')],'strong'),('HCO+ 4-3',356.73424,[('Xb47','spw18')],'strong'),
 ('CS 7-6',342.88285,[('Xb47','spw22')],'blended with H2CS 10_0,10-9_0,9 (342.946) and 29SiO 8-7 (342.981)'),
 ('H2CS 10(0,10)-9(0,9)',342.94644,[('Xb47','spw22')],'blend with CS 7-6'),
 ('HC15N 4-3',344.20016,[('Xb47','spw20')],''),('SO 8_8-7_7',344.31061,[('Xb47','spw20')],''),
 ('HCN v2=1f 4-3',356.2559,[('Xb47','spw18')],'-400 km/s from HCO+ 4-3'),
 ('HCN 3-2',265.88618,b6('spw18'),'strong'),('HCO+ 3-2',267.55763,b6('spw16'),'strong'),
 ('SO 5_6-4_5',251.82577,b6('spw22'),'tentative at S (Xb4f)'),
 ('NS 11/2-9/2 f (253.968; "253.9 bump")',253.96793,b6('spw20'),'unconfirmed; SO2 253.957, c-HCCCH 254.013 competing'),
 ('NS 11/2-9/2 e (253.571)',253.57052,b6('spw20'),'partner of 253.968'),
 ('SiS 14-13',254.10321,b6('spw20'),''),('HC3N 28-27',254.69970,b6('spw20'),'tentative at N'),
 ('HCN v2=1f 3-2',267.1993,b6('spw16'),'-400 km/s from HCO+ 3-2'),
]
EXTRA={('Xb47','spw20'):[345.79599,345.34399],('Xb55','spw22'):[165.19049,165.6798],('Xb47','spw16'):[353.6]}   # CO 3-2 (band edge) + H13CN 4-3 (blend on CO wing); CH3OH 4_-1-4_0 E; H26alpha-ish
KNOWN={}   # known line velocities (absolute radio) per (ous,spw) for line-free exclusion
for lab,f0,lst,_ in LINES:
    for key in lst: KNOWN.setdefault(key,[]).append(f0)
def vel(nu,f0): return C*(1-nu/(f0*1e9))
def clean(y,v,fl):
    yb=y[fl]
    for _ in range(4):
        med=np.median(yb); s=1.4826*np.median(np.abs(yb-med)); yb=yb[abs(yb-med)<3.5*s]
    return np.median(yb)
def make_template(ous,reg,st):
    tspw,tf0=TEMPL[ous]; nu,y=spec(ous,tspw,reg,st)
    if nu is None: return None
    v=vel(nu,tf0); o=np.argsort(v); v,y=v[o],y[o]; ok=np.isfinite(y); v,y=v[ok],y[ok]
    fl=abs(v-VSYS)>500; 
    if fl.sum()<30: return None
    base=clean(y,v,fl); rms=np.std(y[fl&(abs(v-VSYS)<1000)]); y=y-base
    sm=np.convolve(y,np.ones(3)/3,'same'); sm[abs(v-VSYS)>400]=0
    ip=np.argmax(sm)
    if sm[ip]<6*rms/np.sqrt(3): return None
    # contiguous support above 10% of peak around the peak
    i0=ip
    while i0>0 and sm[i0-1]>0.1*sm[ip]: i0-=1
    i1=ip
    while i1<len(sm)-1 and sm[i1+1]>0.1*sm[ip]: i1+=1
    return dict(v=v,t=np.clip(sm,0,None),lo=v[i0],hi=v[i1],vpk=v[ip])
def rho_infl(yb):
    d=yb-yb.mean(); rho=np.sum(d[1:]*d[:-1])/np.sum(d*d); r=max(rho,0); return np.sqrt((1+r)/(1-r))
CAL={'Xb55':0.10,'Xb47':0.10,'Xb4f':0.05,'Xb4d':0.05,'Xb4b':0.05}   # flux-scale (fractional, ALMA Cycle-5 nominal B5/B7 10%, B6 5%)
FWHM12={'Xb55':33.0,'Xb47':16.6,'Xb4f':22.0,'Xb4d':22.0,'Xb4b':22.0}   # 12 m primary beam FWHM [arcsec] (1.13 lambda/D at the mfs freq)
PH=(156.96346,-43.90461)
def pb_of(ous,reg):
    ra,dec=REGIONS[reg][:2]; r=np.hypot((ra-PH[0])*np.cos(np.deg2rad(PH[1])),dec-PH[1])*3600
    pb=np.exp(-4*np.log(2)*(r/FWHM12[ous])**2); return r,pb
def blend_flag(ous,reg,bmaj,bmin):
    # response of the beam (taken circular, geometric mean FWHM) at the position of the nearest other bright region (N,S,Bridge,E1)
    th=np.sqrt(bmaj*bmin); worst=0; who=''
    for k in ('N','S','Bridge','E1','E2'):
        if k==reg: continue
        d=np.hypot((REGIONS[k][0]-REGIONS[reg][0])*np.cos(np.deg2rad(PH[1])),REGIONS[k][1]-REGIONS[reg][1])*3600
        r=np.exp(-4*np.log(2)*(d/th)**2)
        if r>worst: worst,who=r,k
    return worst,who
STRONGK={}   # lines excluded from every line-free mask of the spw: the strong/detected ones (weak candidates are not excluded: conservative)
for lab,f0,lst,note in LINES:
    if note=='strong' or lab in ('CS 7-6','H2CS 10(0,10)-9(0,9)','HOC+ 2-1','CH3OH J-1-J0 E blend'):
        for key in lst: STRONGK.setdefault(key,[]).append(f0)
rows=[]
for lab,f0,lst,note in LINES:
    strong=note=='strong'
    for ous,spw in lst:
        for reg in REGIONS:
            if (ous,spw) not in Z: spec(ous,spw,reg,'ap')
            z=Z[(ous,spw)]; bmaj,bmin=float(z['bmaj']),float(z['bmin'])
            for st in ('ap','pix'):
                nu,y=spec(ous,spw,reg,st)
                if nu is None: continue
                T=make_template(ous,reg,st)
                if T is None: continue
                v=vel(nu,f0); o=np.argsort(v); v,y=v[o],y[o]; ok=np.isfinite(y); v,y=v[ok],y[ok]; dv=abs(np.median(np.diff(v)))
                tfun=lambda vv,s=0.0: np.interp(vv-s,T['v'],T['t'],left=0,right=0)
                half=0.5*(T['hi']-T['lo']); vc=0.5*(T['hi']+T['lo'])
                # line-free mask: exclude +-(half+150 km/s) around EVERY known line of this spw (own and other lines), each located at its
                # velocity in this line's frame, assuming the region's systemic velocity vc (template centre of the region's HCO+ profile)
                free=np.ones(len(v),bool)
                for fk in [f0]+STRONGK.get((ous,spw),[])+EXTRA.get((ous,spw),[]):
                    nuk=fk*1e9*(1-vc/C); vk=C*(1-nuk/(f0*1e9)); free&=abs(v-vk)>half+150
                if free.sum()<30: continue
                base=clean(y,v,free); yy=y-base; rms=np.std(yy[free]); infl=rho_infl(yy[free])
                t=tfun(v); n_t=np.sum(t**2)
                if n_t<=0: continue
                win=(v>=T['lo'])&(v<=T['hi']); nw=win.sum()
                if win.sum()<4 or (v.min()>T['lo']+50) or (v.max()<T['hi']-50): continue
                IT=np.sum(t*yy)/n_t*1.0       # amplitude in units of template (template normalised to unit peak? -> rescale to unit integral below)
                tint=np.sum(t)*dv
                IT=np.sum(t*yy)/n_t/ (1.0) ; # template has peak-unit scaling of flux density; amplitude a means y=a*t(v); integral = a*tint
                # make template unit-peak so that integral I = a*tint
                tn=t/ t.max(); IT=np.sum(tn*yy)/np.sum(tn**2); I_T=IT*np.sum(tn)*dv; sstat=rms/np.sqrt(np.sum(tn**2))*np.sum(tn)*dv*infl
                W=np.sum(yy[win])*dv; sW=rms*np.sqrt(nw)*dv*infl
                # empirical: slide template over line-free positions
                emp=[]
                for s in np.arange(-3500,3500,max(dv,half/2)):
                    tt=np.interp(v-vc-s,T['v']-vc,T['t'],left=0,right=0); tt=tt/T['t'].max()
                    sup=tt>0.02
                    if sup.sum()<4 or not free[sup].all(): continue
                    # support must lie fully inside the sampled channels (no truncated template)
                    cvr=(v.min()<=vc+s-half) and (v.max()>=vc+s+half)
                    if not cvr or abs(s)<2*half*0.8: continue
                    emp.append(np.sum(tt*yy)/np.sum(tt**2)*np.sum(tt)*dv)
                emp=np.array(emp); n_emp=len(emp)
                semp=np.std(emp) if n_emp>=4 else np.nan
                sb=np.sqrt(max(semp**2-sstat**2,0)) if n_emp>=4 else np.nan
                stot=np.sqrt(sstat**2+(sb if n_emp>=4 else 0)**2)
                pk_sm=np.convolve(yy,np.ones(3)/3,'same'); wmask=win; pk=pk_sm[wmask].max()
                nu_obs=f0*(1-VSYS/C); bl,who=blend_flag(ous,reg,bmaj,bmin); roff,pb=pb_of(ous,reg)
                sig_cal=np.hypot(CAL[ous]*abs(I_T),(0.10 if pb<0.75 else 0.05 if pb<0.9 else 0.0)*abs(I_T))
                rows.append(dict(ous=ous,beam=f'{bmaj:.2f}x{bmin:.2f}',line=lab,rest_GHz=f0,region=reg,spec=st,dv_kms=dv,win_lo=T['lo'],win_hi=T['hi'],v_template_peak=T['vpk'],
                    rms_chan_mJy=rms*1e3,peak_sm3_mJy=pk*1e3,I_T=I_T,sig_stat=sstat,n_emp=n_emp,sig_emp=semp,sig_base=sb,sig_tot=stot,snr=I_T/stot,
                    UL3=3*stot if I_T/stot<3 else np.nan,W_sum=W,sig_W=sW,sig_cal=sig_cal,
                    I_T_Kkms=1.222e3*(I_T*1e3)/(nu_obs**2*bmaj*bmin),pb=pb,r_off_arcsec=roff,blend_resp=bl,blend_with=who,note=note))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s12_line_table.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',900,'display.max_colwidth',30)
w=df[(df.spec=='ap')&(~df.note.eq('strong'))]
print(w[['ous','line','region','win_lo','win_hi','I_T','sig_stat','sig_emp','n_emp','sig_tot','snr','UL3','blend_resp']].round(3).to_string())
print(len(df),'rows')
