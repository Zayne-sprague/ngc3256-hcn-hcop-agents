# Q4b: HCN v2=1f (IR-pumping tracer, Eu~1050 K) vs ground-state HCN, same J, same OUS/beam, at N, S, Bridge, Soutflow. J=3-2 (Xb4b 0.24", Xb4d 0.95", Xb4f 5.7x3.5") and 4-3 (Xb47 4.1x2.7").
# v2=1f lies 403 km/s from HCO+ (blue wing of HCO+): the HCO+ wing integral in the identical offset window is measured and subtracted; line-free linear baseline.
# Window = +-1 FWHM of the ground-state HCN Gaussian fit (out/line_table_v2.csv) about its centroid. Output out/q4/hcn_vib.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_lib import *
H,K,CL=6.62607e-27,1.380649e-16,2.99792458e10
LT=pd.read_csv('/workspace/out/line_table_v2.csv'); LT=LT[LT.kind=='peakpix']
SP={('Xb4b',3):dict(gnd=('spw27',265.886434),hco=('spw25',267.557626),vib=('spw25',267.199390)),
    ('Xb4d',3):dict(gnd=('spw27',265.886434),hco=('spw25',267.557626),vib=('spw25',267.199390)),
    ('Xb4f',3):dict(gnd=('spw18',265.886434),hco=('spw16',267.557626),vib=('spw16',267.199390)),
    ('Xb47',4):dict(gnd=('spw16',354.505478),hco=('spw18',356.734223),vib=('spw18',356.255550))}
EVIB=1024.4          # K, v2=1 (712 cm-1) energy above ground; rotational energy of J_u the same in both states (B~2.13 K) -> cancels in the population ratio
def Aul(J,fv,fg,Ag):  # A(vib f, J->J-1) from the v=0 value: nu^3 scaling and l=1 line strength (J^2-l^2)/J vs J
    return Ag*(fv/fg)**3*((J*J-1.)/J)/J
AG={3:10**-3.07662*1.0,4:10**-2.68602*1.0}   # curated Splatalogue/JPL log A (v=0 3-2: -3.077 ; 4-3: -2.686 s-1)
def win(v,s,v0,w,free,order=1):
    m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    if m.sum()<8: return np.nan,np.nan,np.nan,0
    p=np.polyfit(v[m]-v0,s[m],order); r=s-np.polyval(p,v-v0); sg=robust_sigma(r[m]); k=(v>=v0-w)&(v<=v0+w)
    dv=abs(np.median(np.diff(v))); return r[k].sum()*dv,sg*dv*np.sqrt(k.sum()),sg,k.sum()
rows=[]
for (ous,J),d in SP.items():
    for reg in ['N','S','Bridge']+(['Soutflow'] if ous in ('Xb4f','Xb47') else []):   # Soutflow ground FWHM 250-450 km/s at Xb4d/Xb4b: wing model not separable
        g=LT[(LT.ous==ous)&(LT.region==reg)&(LT.line==f'HCN {J}-{J-1}')].iloc[0]; v0=g.g_v0; fw=g.g_fwhm; w=fw
        bmaj,bmin=g.bmaj,g.bmin; out={}
        for key,(spw,f0) in d.items():
            fr,s,_,nb=spectrum(ous,spw,reg); out[key]=(vel(fr,f0),s*1e3,f0)   # mJy/beam
        fv=d['vib'][1]; fg=d['gnd'][1]; fh=d['hco'][1]
        vv,sv,_=out['vib']; vg,sg_,_=out['gnd']; vh,sh,_=out['hco']
        # ground-state HCN integral (K km/s): window +-1 FWHM, baseline from channels outside +-3 FWHM within +-900 km/s
        Ig,eg,_,_=win(vg,sg_,v0,w,[(v0-900,v0-3*fw),(v0+3*fw,v0+900)])
        # vib line sits on the red wing of HCO+ (offset off km/s): joint linear fit in the HCO+ frame of  baseline(1) + A_w*exp(-|v-v0|/tau) (|v-v0|>1.5 FWHM, HCO+ core masked)
        # + A_v*Gauss(v0+off, sigma=FWHM_gnd/2.355); tau scanned 40-300 km/s (profile chi2, error from the chi2 minimum + 1 within-tau spread)
        off=(1-fv/fh)*C
        Iv,ev,Iw,ew,nch,tau_b=fit_target(vh,sh,v0,fw,off)
        # convert mJy/beam km/s -> K km/s with the native beam
        kvib=kb(fv,bmaj,bmin)*1e-3; kg=kb(fg,bmaj,bmin)*1e-3; khc=kb(fh,bmaj,bmin)*1e-3
        Ivk,evk=Iv*kvib,ev*kvib; Igk,egk=Ig*kg,eg*kg; Iwk,ewk=Iw*khc,ew*khc
        exc=Ivk; eexc=evk   # the joint fit already separates the HCO+ wing; Iwk (wing level in the window) is informational
        Ag=AG[J]; Av=Aul(J,fv,fg,Ag)
        Nu=lambda W,f,A: 8*np.pi*K*(f*1e9)**2*W*1e5/(H*CL**3*A)
        Nv_lim=Nu(max(exc,0)+3*eexc,fv,Av); N0=Nu(Igk,fg,Ag); ratio_lim=Nv_lim/N0
        Tvib_lim=EVIB/np.log(1./ratio_lim) if ratio_lim<1 else np.nan
        rows.append(dict(ous=ous,J=f'{J}-{J-1}',region=reg,beam=f'{bmaj:.2f}x{bmin:.2f}',v0=v0,fwhm=fw,I_gnd_Kkms=Igk,eI_gnd=egk,I_vib_Kkms=Ivk,eI_vib=evk,I_HCOpwing_Kkms=Iwk,eI_wing=ewk,I_vib_minus_wing=exc,e_excess=eexc,
                         vib_over_gnd_Tb_3sig_upper=(max(exc,0)+3*eexc)/Igk*(fg/fv)**2,Tb_ratio_best=exc/Igk*(fg/fv)**2,N_vib_over_N_gnd_3sig=ratio_lim,Tvib_3sig_upper_K=Tvib_lim,nch_baseline=nch))
R=pd.DataFrame(rows); R.round(4).to_csv(OUT+'hcn_vib.csv',index=False); pd.set_option('display.width',280); print(R.round(3).to_string())
