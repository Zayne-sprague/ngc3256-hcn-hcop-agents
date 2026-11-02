# Q14n (fix 3,4): abundance proxies and SiO informativeness at the common 6.3x3.8" beam (all inputs from q14e/q14h/q14i/q14j/q14k outputs).
# N(CH3OH) (A+E): LTE at Tex 10-30 K and non-LTE (RADEX, n 10^5.5-10^6.5, Tk 20-50 K, scaled in I254; thin) ; N(H2) dust (Td 20-50 K); N(H13CN) limit from alpha (Q4b) ; SiO: Harada-predicted line intensities vs our limits.
# Output out/q14e/summary_*.csv. python3 -I q14n_summary.py
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q14f_lte as q
exec(open('/workspace/scripts/q14i_abund.py').read().split("import q14f_lte as q")[0].split("# validation vs LAMDA")[0].split("def lamda_W")[0].replace("from q14e_lib import *","import numpy as np, pandas as pd"))
rng=np.random.default_rng(9)
I_=pd.read_csv('/workspace/out/q14e/abund_intensities.csv'); H2=pd.read_csv('/workspace/out/q14e/h2_dust.csv'); H2=H2[H2.beam=='6.3x3.8']
RS=pd.read_csv('/workspace/out/q14e/radex_ch3oh_solutions.csv'); G=pd.read_csv('/workspace/out/q14e/radex_ch3oh_grid.csv')
def Ig(kind,ous,beam,reg): r=I_[(I_.kind==kind)&(I_.ous==ous)&(I_.beam==beam)&(I_.region==reg)].iloc[0]; return r.I,r.eI
RG=['N','S','Bridge','Wclump','Soutflow']
# non-LTE column per K km/s of 254: thin, depends on (Tk,n); take solutions at I254=0.57 -> N/I; band n 10^5.5-10^6.5, Tk 20-50
S=RS[(RS.I254_obs==0.57)&(RS.log_n>=5.5)&(RS.log_n<=6.5)&(RS.Tkin>=20)&(RS.Tkin<=50)]; NperI=10**S.log_Ntot/0.57
NperI_lo,NperI_hi=NperI.min(),NperI.max(); print('RADEX N_tot per K km/s of 254 (n 1e5.5-1e6.5, Tk 20-50): %.3g - %.3g ; LTE 10/20/30 K: %s'%(NperI_lo,NperI_hi,[1/q.W254(T) for T in (10,20,30)]))
rows=[]
for reg in RG:
    I,e=Ig('CH3OH','Xb4f','6.3x3.8',reg); Id,ed=Ig('CH3OH','Xb4d','6.3x3.8',reg)
    for Td in (20.,30.,50.):
        h=H2[(H2.region==reg)&(H2.Td==Td)].iloc[0]
        for lab,npi in (('LTE 10K',1/q.W254(10.)),('LTE 20K',1/q.W254(20.)),('LTE 30K',1/q.W254(30.)),('nonLTE lo',NperI_lo),('nonLTE hi',NperI_hi)):
            N=I*npi; Nd=Id*npi
            rows.append(dict(region=reg,Td=Td,conv=lab,I254_Xb4f=I,eI=e,N_CH3OH=N,eN=e*npi,NH2=h.NH2_med,NH2_16=h.NH2_16,NH2_84=h.NH2_84,X_CH3OH=N/h.NH2_med,eX_stat=e*npi/h.NH2_med,X_lo=N/h.NH2_84,X_hi=N/h.NH2_16,X_Xb4d=Nd/h.NH2_med))
X=pd.DataFrame(rows); X.to_csv('/workspace/out/q14e/summary_X_CH3OH.csv',index=False); pd.set_option('display.width',250); pd.set_option('display.max_rows',500)
print(X[(X.Td==30)].round(12).to_string())
# S/N contrast of X(CH3OH)=N/N(H2) with MC: Tex_N,Tex_S choices (Td same for both, or Td_N>Td_S)
IN,eN=Ig('CH3OH','Xb4f','6.3x3.8','N'); IS,eS=Ig('CH3OH','Xb4f','6.3x3.8','S'); out=[]
for (Tn,Ts,TdN,TdS) in [(10,10,30,30),(20,10,30,30),(30,10,30,30),(30,10,40,25),(20,20,30,30),(30,30,30,30),(10,10,50,20)]:
    r=[]
    for k in range(20000):
        iN=max(rng.normal(IN,eN),1e-4); iS=rng.normal(IS,eS)
        hN=H2[(H2.region=='N')&(H2.Td==float(TdN if TdN in(20,30,50) else 30))].iloc[0]; hS=H2[(H2.region=='S')&(H2.Td==float(TdS if TdS in(20,30,50) else 30))].iloc[0]
        # N(H2) scaling for Td not on grid: B_nu(Td) ratio
        def nh2(h,Td,Tg): 
            B=lambda T: 1./np.expm1(4.799e-11*346.53e9/T); return rng.normal(h.NH2_med,0.5*(h.NH2_84-h.NH2_16))*B(Td)**-1*B(Tg)
        r.append((iS*(1/q.W254(Ts))/nh2(hS,TdS,float(TdS if TdS in(20,30,50) else 30)))/(iN*(1/q.W254(Tn))/nh2(hN,TdN,float(TdN if TdN in(20,30,50) else 30))))
    r=np.array(r); out.append(dict(Tex_N=Tn,Tex_S=Ts,Td_N=TdN,Td_S=TdS,XS_over_XN_med=np.median(r),lo16=np.percentile(r,16),hi84=np.percentile(r,84),lo5=np.percentile(r,5),P_gt1=np.mean(r>1),P_gt3=np.mean(r>3)))
C=pd.DataFrame(out); C.to_csv('/workspace/out/q14e/summary_XS_over_XN.csv',index=False); print(C.round(3).to_string())
# CS columns and N(CS)/N(H2)
cs=[]
for reg in RG:
    I,e=Ig('CS 7-6','Xb47','6.3x3.8',reg)
    for Td in (30.,):
        h=H2[(H2.region==reg)&(H2.Td==Td)].iloc[0]
        for T in (20.,30.,50.):
            N=I/W('CS76',T); cs.append(dict(region=reg,Tex=T,Td=Td,I_CS76=I,eI=e,N_CS=N,eN=e/W('CS76',T),NH2=h.NH2_med,X_CS=N/h.NH2_med))
CSd=pd.DataFrame(cs); CSd.to_csv('/workspace/out/q14e/summary_X_CS.csv',index=False); print(CSd.round(12).to_string())
# H13CN-based lower limits on N(CH3OH)/N(H13CN)
alpha={'N':0.105,'S':0.25,'Bridge':0.15}; hcn43={'N':(7.30,0.22),'S':(3.03,0.08),'Bridge':(4.39,0.08)}; h13=[]
for reg in ('N','S','Bridge'):
    I,e=Ig('CH3OH','Xb4f','6.3x3.8',reg); IH=alpha[reg]*hcn43[reg][0]
    for TexC in (10.,20.,30.):
        for TexH in (20.,30.,50.,100.):
            Nc=I/q.W254(TexC); Nh=IH/W('H13CN43',TexH)
            h13.append(dict(region=reg,alpha3=alpha[reg],I_H13CN_lim3=IH,Tex_CH3OH=TexC,Tex_H13CN=TexH,N_CH3OH=Nc,N_H13CN_lim3=Nh,CH3OH_over_H13CN_gt=Nc/Nh,CH3OH_over_HCN_gt_12C13C70=Nc/(Nh*70)))
Hh=pd.DataFrame(h13); Hh.to_csv('/workspace/out/q14e/summary_H13CN.csv',index=False); print(Hh.round(4).to_string())
# SiO informativeness: predicted intensity if N(SiO)/N(CS)=Harada ratio (N 0.062, S 0.057)
sio=[]; HR={'N':0.062,'S':0.057}
for reg in ('N','S'):
    ics,_=Ig('CS 7-6','Xb47','6.3x3.8',reg); l29=Ig('29SiO 8-7 (H2CS=0)','Xb47','6.3x3.8',reg); l29j=Ig('29SiO 8-7 (joint CS+H2CS+29SiO)','Xb47','6.3x3.8',reg)
    l30=Ig('30SiO 6-5','Xb4f','6.3x3.8',reg); l30d=Ig('30SiO 6-5','Xb4d','6.3x3.8',reg)
    for T in (10.,20.,30.,50.):
        NCS=ics/W('CS76',T); NSiO=HR[reg]*NCS
        for nm,key,iso,(I,e) in (('29SiO 8-7 (H2CS=0 fit)','29SiO87',20,l29),('29SiO 8-7 (joint fit)','29SiO87',20,l29j),('30SiO 6-5 Xb4f','30SiO65',30,l30),('30SiO 6-5 Xb4d','30SiO65',30,l30d)):
            Ip=NSiO/iso*W(key,T); sio.append(dict(region=reg,Tex=T,line=nm,I_pred_K_kms=Ip,sigma=e,pred_over_sigma=Ip/e,lim3=max(I,0)+3*e,factor_lim3_over_pred=(max(I,0)+3*e)/Ip))
Sd=pd.DataFrame(sio); Sd.to_csv('/workspace/out/q14e/summary_SiO.csv',index=False); print(Sd.round(4).to_string())
