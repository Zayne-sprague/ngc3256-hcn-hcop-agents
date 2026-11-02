# Q14l (fix 1,5): final Harada+2018 matched test with (a) whole-line template-deblended 254.015 intensities (q14h) at Xb4d 1.7" (+ feathering with the Xb4f-Xb4d deficit at 6.3x3.8),
# (b) the +-150 km/s window fraction of the line profile (Harada integrate +-150 km/s about 2775 km/s), (c) LTE conversion at Harada's own Trot (241.79/96.74) and (d) a non-LTE (RADEX) conversion constrained by
# R(254/Q-blend) (q14b) and Harada's R(241.79/96.74).  Output out/q14e/harada_final.csv, radex_posterior.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q14e_lib import *
from scipy.interpolate import RegularGridInterpolator as RGI
rng=np.random.default_rng(5)
D=pd.read_csv('/workspace/out/q14e/deblend_all.csv'); G=pd.read_csv('/workspace/out/q14e/radex_ch3oh_grid.csv'); L=pd.read_csv('/workspace/out/q14f/lte_ratios.csv'); L=L[L.cmb_contrast].set_index('T')
def dg(ous,beam,reg,comp='CH3OH',tmpl='HCO+',model='all'):
    r=D[(D.cluster=='A')&(D.ous==ous)&(D.beam==beam)&(D.region==reg)&(D.comp==comp)&(D.template==tmpl)&(D.model==model)]; return float(r.I.iloc[0]),float(r.eI.iloc[0])
# (b) +-150 km/s fraction of the HCO+ 3-2 template at 1.7" (Xb4d)
fr,a,cw,pix,beam,nb=cut('Xb4d','spw25',(1.7,1.7,0.),half=8.); v=vel(fr,267.557626); F150={}
for reg in ('N','S'):
    s=sample(a,cw,REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']); p=hco_fit(v,s); v0,fw=p[1],p[2]
    side=(np.abs(v-v0)>2.0*fw)&(np.abs(v-v0)<900); b=np.polyfit(v[side]-v0,s[side],1); t=s-np.polyval(b,v-v0); t=np.where(np.abs(v-v0)<2.2*fw,t,0.)
    F150[reg]=t[np.abs(v-2775)<150].sum()/t.sum()
print('fraction of the line within +-150 km/s of 2775:',F150)
# (d) RADEX posterior on a (Tkin,log n) grid; N_tot such that I254 = I_obs (Xb4f 6.3x3.8 deblended)
Tk=np.array(sorted(G.Tkin.unique())); ln=np.array(sorted(G.log_n.unique())); lN=np.array(sorted(G.log_Ntot.unique()))
def cube(col): return np.array([[[G[(G.Tkin==t)&(G.log_n==n)&(G.log_Ntot==N)][col].iloc[0] for N in lN] for n in ln] for t in Tk])
CUBES={c:cube(c) for c in ('I254','Iblend','I96','I241')}
def post(Iobs,Rb_obs,Rb_err,R241_obs,R241_err,tk_lo=20.):
    tf=np.linspace(Tk.min(),Tk.max(),41); nf=np.linspace(ln.min(),ln.max(),71); out=[]
    itp={c:RGI((Tk,ln,lN),np.log(np.clip(CUBES[c],1e-12,None)),bounds_error=False,fill_value=None) for c in CUBES}
    for t in tf:
        for n in nf:
            Ns=np.linspace(lN.min(),lN.max(),200); pts=np.array([[t,n,x] for x in Ns])
            I254=np.exp(itp['I254'](pts)); k=np.argmin(abs(np.log(I254/Iobs))); 
            if abs(np.log(I254[k]/Iobs))>0.1: continue
            pt=pts[k:k+1]; e=lambda c: float(np.exp(itp[c](pt))[0])
            Rb=e('I254')/e('Iblend'); R96=e('I96')/e('I254'); R241=e('I241')/e('I96')
            chi=((Rb-Rb_obs)/Rb_err)**2+((R241-R241_obs)/R241_err)**2
            out.append(dict(Tkin=t,log_n=n,log_Ntot=Ns[k],Rblend=Rb,R96_254=R96,R241_96=R241,chi2=chi))
    P=pd.DataFrame(out); P['w']=np.exp(-0.5*(P.chi2-P.chi2.min())); return P
def wq(x,w,q): i=np.argsort(x); c=np.cumsum(w[i])/w.sum(); return np.interp(q,c,np.asarray(x)[i])
rows=[]
CONS={'S':dict(Iobs=0.569,Rb=0.31,Rbe=0.075,R241=0.70,R241e=0.10,I1p7=dg('Xb4d','1.7','S'),def63=(dg('Xb4f','6.3x3.8','S')[0]-dg('Xb4d','6.3x3.8','S')[0],np.hypot(dg('Xb4f','6.3x3.8','S')[1],dg('Xb4d','6.3x3.8','S')[1])),harada=(6.22,0.66),hco=(43.69,2.0),trotH=11.05),
      'N':dict(Iobs=0.147,Rb=0.12,Rbe=0.05,R241=1.89,R241e=0.45,I1p7=dg('Xb4d','1.7','N'),def63=(dg('Xb4f','6.3x3.8','N')[0]-dg('Xb4d','6.3x3.8','N')[0],np.hypot(dg('Xb4f','6.3x3.8','N')[1],dg('Xb4d','6.3x3.8','N')[1])),harada=(3.09,0.66),hco=(91.56,2.0),trotH=17.3)}
POST={}
for reg,c in CONS.items():
    P=post(c['Iobs'],c['Rb'],c['Rbe'],c['R241'],c['R241e']); P.to_csv(f'/workspace/out/q14e/radex_posterior_{reg}.csv',index=False); POST[reg]=P
    w=P.w.values; print(reg,'chi2min',P.chi2.min().round(2),'n_cells',len(P))
    for col in ('log_n','Tkin','log_Ntot','R96_254','Rblend','R241_96'): print('  ',col,[round(wq(P[col],w,q),3) for q in (0.025,0.16,0.5,0.84,0.975)])
    I1,e1=c['I1p7']; d,ed=c['def63']; ff=F150[reg]
    for variant,(I,e) in {'Xb4d 1.7" template':(I1,e1),'feathered (Xb4d 1.7"+[Xb4f-Xb4d]6.3x3.8)':(I1+d,np.hypot(e1,ed))}.items():
        Iw=I*ff; ew=e*ff
        # LTE at Harada's Trot and at 8-30 K
        for tag,Rl in (('LTE Trot(Harada 241/96)=%.0f K'%c['trotH'],np.interp(c['trotH'],L.index.values,L.R96_254.values)),('LTE 8-30 K (6.6-4.9)',5.7)):
            rows.append(dict(region=reg,variant=variant,conv=tag,I254_win=Iw,e=ew,R96_254=Rl,pred96=Iw*Rl,e_pred=ew*Rl,harada=c['harada'][0],e_harada=c['harada'][1],pred_over_harada=Iw*Rl/c['harada'][0],nsig=(Iw*Rl-c['harada'][0])/np.hypot(ew*Rl,c['harada'][1]),ratio_to_HCO_Harada=Iw*Rl/c['hco'][0]))
        Rs=wq(P.R96_254,w,[0.16,0.5,0.84]); 
        for lab,Rl in zip(('16%','50%','84%'),Rs):
            rows.append(dict(region=reg,variant=variant,conv='RADEX posterior R96/254 '+lab,I254_win=Iw,e=ew,R96_254=Rl,pred96=Iw*Rl,e_pred=ew*Rl,harada=c['harada'][0],e_harada=c['harada'][1],pred_over_harada=Iw*Rl/c['harada'][0],nsig=(Iw*Rl-c['harada'][0])/np.hypot(ew*Rl,c['harada'][1]),ratio_to_HCO_Harada=Iw*Rl/c['hco'][0]))
HF=pd.DataFrame(rows); HF.to_csv('/workspace/out/q14e/harada_final.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_columns',30); print(HF.round(3).to_string())
