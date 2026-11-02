# LTE (optically thin, common Trot) fit to the matched-beam CH3OH line intensities (out/ch3oh_matched.csv); profile chi2 in Trot; N_tot from A.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from ch3oh import c
df=pd.read_csv('/workspace/out/ch3oh_matched.csv')
use=df[((df.ous=='Xb55'))|((df.ous=='Xb4f')&(df.line.isin(['2(0)-1(-1) E','5(-2,4)-4(-1,3) E'])))].copy()   # drop 6(-1,5)-5(-2,4): 6 MHz from spw edge
CAL=0.10
use['e']=np.sqrt(use.eI_Kkms**2+(CAL*use.I_Kkms.clip(lower=0))**2)
def pred(T):
    out=[]
    for r in use.itertuples():
        fr0s={'Q J=1-3':(165.050175,165.061130,165.099240)}.get(r.line,(r.nu_rest,))
        w=0
        for f0 in fr0s:
            l=c[(c.f0-f0).abs()<2e-4].iloc[0]; w+=l.f0*l.sijmu2*np.exp(-l.Eu/T)
        out.append(w)
    return np.array(out)
Ts=np.array([3,4,5,6,7,8,9,10,12,14,17,20,25,30,40,60,100.]); res=[]
for T in Ts:
    w=pred(T); y=use.I_Kkms.values; e=use.e.values
    A=np.sum(w*y/e**2)/np.sum(w*w/e**2); chi2=np.sum(((y-A*w)/e)**2); res.append((T,A,chi2))
res=np.array(res); i=res[:,2].argmin(); print('chi2 profile (T, A, chi2):'); print(np.round(res,3))
dm=res[:,2]-res[i,2]; ok=Ts[dm<1.0]; ok3=Ts[dm<9.0]
print(f'Trot best {Ts[i]} K; 1-sigma range {ok.min()}-{ok.max()}; 3-sigma range {ok3.min()}-{ok3.max()}; chi2min={res[i,2]:.1f} for {len(use)-2} dof')
T=Ts[i]; w=pred(T); use['pred_Kkms']=res[i,1]*w
print(use[['line','ous','Eu','I_Kkms','e','pred_Kkms']].round(3).to_string())
# column density: N_tot = 8 pi k nu^2 I/(h c^3 A_ul) * Q/g_u * exp(Eu/T) for the 2(0)-1(-1) line (CGS); partition fn of E-CH3OH ~ half of total; use Q(T) from approx. Q_total(T)=1.2*T^1.5 ... give N_u only
k=1.380649e-16;h=6.62607e-27;cc=2.99792458e10
l=c[(c.f0-254.015377).abs()<2e-4].iloc[0]; Aij=10**l.logA; nu=l.f0*1e9
I254=use[(use.line=='2(0)-1(-1) E')].I_Kkms.iloc[0]*1e5   # K cm/s
Nu=8*np.pi*k*nu**2*I254/(h*cc**3*Aij)
print(f'N_u(2(0)_E, Eu=20K) = {Nu:.2e} cm^-2 (beam-averaged over 9.4x6.3"), g_u=5')
Q123=use[(use.line=='Q J=1-3')].I_Kkms.iloc[0]
use.to_csv('/workspace/out/ch3oh_ltefit_points.csv',index=False)
pd.DataFrame(res,columns=['Trot','A','chi2']).to_csv('/workspace/out/ch3oh_ltefit_profile.csv',index=False)
