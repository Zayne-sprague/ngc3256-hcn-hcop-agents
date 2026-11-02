# Q14p: SO 8_8-7_7 (344.3106 GHz, Eu 87.5 K) 3-sigma limits (6.3x3.8", Xb47, SO/HCO+4-3 limits of q14a) -> N(SO) limits (LAMDA so@lique.dat, thin LTE) and N(SO)/N(CS) at common Tex. python3 -I q14p_so.py
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q14f_lte as q
lv,tr=q._lam('/opt/lamda/so@lique.dat'); r=tr[np.abs(tr[:,4]-344.310612)<0.01][0]; up=int(r[1])-1
def WSO(T):
    Q=np.sum(lv[:,1]*np.exp(-lv[:,0]/T)); Nu=lv[up,1]*np.exp(-lv[up,0]/T)/Q; nu=r[4]
    w=q.hP*q.cc**3*r[3]*Nu/(8*np.pi*q.kB*(nu*1e9)**2)/1e5; x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; return w*(1-(np.exp(x)-1)/(np.exp(xb)-1))
B=pd.read_csv('/workspace/out/q14/beamseq.csv'); L=pd.read_csv('/workspace/out/q14/beamseq_lines.csv')
exec(open('/workspace/scripts/q14i_abund.py').read().split("import q14f_lte as q")[0].split("def lamda_W")[0].replace("from q14e_lib import *","import numpy as np, pandas as pd"))
I_=pd.read_csv('/workspace/out/q14e/abund_intensities.csv'); rows=[]
for reg in ('N','S','Bridge'):
    b=B[(B.beam=='6.3x3.8')&(B.ous=='Xb47')&(B.region==reg)&(B.num=='SO 8_8-7_7')&(B.den=='HCO+ 4-3')].iloc[0]
    h=L[(L.beam=='6.3x3.8')&(L.ous=='Xb47')&(L.region==reg)&(L.line=='HCO+ 4-3')].iloc[0].I
    lim3=max(b.ratio,0)*h+3*b.eratio*h; cs=I_[(I_.kind=='CS 7-6')&(I_.beam=='6.3x3.8')&(I_.region==reg)].iloc[0].I
    for T in (20.,30.,50.):
        NSO=lim3/WSO(T); NCS=cs/W('CS76',T); rows.append(dict(region=reg,Tex=T,I_SO_lim3=lim3,N_SO_lim3=NSO,N_CS=NCS,SO_over_CS_lim3=NSO/NCS))
S=pd.DataFrame(rows); S.to_csv('/workspace/out/q14e/summary_SO.csv',index=False); print(S.round(3).to_string())
