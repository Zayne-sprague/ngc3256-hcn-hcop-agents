# Q3 step 7: dense-gas luminosity / mass of the Soutflow blue wing (2300-2575 km/s, i.e. Vsys-475..-200) and non-detection limits at N (wing/core). Assumptions printed.
import sys; sys.path.insert(0,'/workspace/scripts')
from q3_pos import *
import pandas as pd
DL=40.; z=0.0093; ALPHA=10.   # Mpc (Pereira-Santaella+2024); Gao & Solomon 2004 alpha_HCN(1-0) = 10 Msun/(K km/s pc^2), assumed also for HCO+
Lp=lambda S,nu_obs: 3.25e7*S/nu_obs**2*DL**2*(1+z)**-3
O=pd.read_csv('/workspace/out/q3/q3_r31_vs_M18.csv'); b=O[(O.pos=='M18 South-OF-Blue')]
res=[]
for sp,nuo in (('HCN',265.886434/(1+VS/C)),('HCO+',267.557626/(1+VS/C))):
    cc=cube(sp,'Xb4d',crop=8.); m0,n=wsum(cc,2300,2575); sg=wnoise(cc,2300,2575)[0]
    box=(cc['dx']>-1.7)&(cc['dx']<1.3)&(cc['dy']>-6.0)&(cc['dy']<-3.3); ppb=OM(cc['beam'])/cc['pix']**2
    S=m0[box].sum()/ppb; eS=np.sqrt(box.sum()/ppb)*np.median(sg[box]); L=Lp(S,nuo)
    r=b[(b.line==sp)&(b.vsys==2775)].iloc[0]; r31,er=r.r31,np.hypot(r.e_r31,(b[b.line==sp].r31.max()-b[b.line==sp].r31.min())/2)
    L10=L/r31; M=ALPHA*L10
    res.append(dict(line=sp,S_Jykms=S,eS=eS,Lprime_3to2=L,eL=Lp(eS,nuo),r31=r31,e_r31=er,Lprime_1to0=L10,Mdense_Msun=M,e_M_stat=M*np.hypot(eS/S,er/r31)))
    # N nucleus wing limits
    for wn,(lo,hi) in {'blue 200-350':(VS-350,VS-200),'red 200-350':(VS+200,VS+350)}.items():
        w,_=wsum(cc,lo,hi); sgw=wnoise(cc,lo,hi)[0]; core,_=wsum(cc,VS-150,VS+150)
        x,y=POS['N']; res.append(dict(line=sp,note=f'N wing {wn} native Xb4d',I=samp(cc,w,x,y),e=samp(cc,sgw,x,y),core=samp(cc,core,x,y),lim3_over_core=3*samp(cc,sgw,x,y)/samp(cc,core,x,y)))
R=pd.DataFrame(res);pd.set_option('display.width',250);print(R.to_string());R.to_csv('/workspace/out/q3/q3_mass_limits.csv',index=False)
for k in ('M18 South-OF-Blue','M18 OF-Red','M18 North-OF-Blue','Soutflow'): print(k,[round(a,2) for a in POS[k]])
