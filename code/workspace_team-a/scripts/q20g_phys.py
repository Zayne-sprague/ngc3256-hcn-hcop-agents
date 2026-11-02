# Q20g: (1) aperture-referenced filling-factor chain (eta_core, eta_1.7, eta_10) incl. Q5 size interval and unresolved floor;
#       (2) physical quantities of the S/N core solutions from the Q20b SMC posteriors (n, N, X-assumed N_H2, path length, mass, nT, tau);
#       (3) CO(2-1) Tpk (Brunetti+2021 Table 2) consistency; (4) dust-continuum column cross-check.
import numpy as np, pandas as pd, glob
from scipy.ndimage import map_coordinates
hv=12.8
def J(T,Tbg=2.73): return hv/np.expm1(hv/T)-hv/np.expm1(hv/Tbg)
def Jco(T,Tbg=2.73): h=11.07; return h/np.expm1(h/T)-h/np.expm1(h/Tbg)   # CO(2-1) 230.5 GHz h nu/k=11.07 K
BEAM=0.24; PC=194.
print('== (1) filling chain: eta_beam(b)=eta_core*s^2/(s^2+b^2)  (Gaussian core FWHM s in an elliptical-gaussian-equivalent circular beam b)')
rows=[]
for reg,sizes in (('S',(0.25,0.30,0.34)),('N',(0.26,0.27,0.29))):
  for s in sizes+(0.10,):
    for ec in (0.2,0.6,1.0):
        rows.append(dict(region=reg,s_arcsec=s,s_pc=s*PC,eta_core=ec,eta_1p7=ec*s**2/(s**2+1.7**2),eta_10=ec*s**2/(s**2+100.),dil_025=s**2/(s**2+BEAM**2)))
T=pd.DataFrame(rows); T.to_csv('/workspace/out/q20/q20g_eta_chain.csv',index=False); print(T.round(4).to_string())
print('geometric dilution only (eta_core=1): s=0.25-0.34 -> eta_1.7=%.3f-%.3f, eta_10=%.1e-%.1e; floor s=0.10 -> %.4f, %.1e'%(0.25**2/(0.25**2+1.7**2),0.34**2/(0.34**2+1.7**2),0.25**2/(0.25**2+100),0.34**2/(0.34**2+100),0.1**2/(0.1**2+1.7**2),0.1**2/(0.1**2+100)))
print('eta_core needed for eta_1.7=0.01-0.02: s=0.25 ->',[round(e*(0.25**2+1.7**2)/0.25**2,2) for e in (0.01,0.02)],' s=0.34 ->',[round(e*(0.34**2+1.7**2)/0.34**2,2) for e in (0.01,0.02)])
print('\n== (3) CO(2-1) Tpk (Brunetti+2021 Table 2, 55 pc) vs eta_CO>=Tpk/[J_CO(Tk)]')
for reg,tp in (('N',(5.1,11,20,37)),('S',(6.5,13,18,37))):
    for Tk in (17,20,25,30,50,80):
        print(reg,'Tk',Tk,'J_CO',round(Jco(Tk),1),' eta_CO,min for Tpk=',tp,':',[round(t/Jco(Tk),2) for t in tp])
print('\n== (1b) unresolved/floor case: T_core=Tb_025/dil(s); Tk,min for eta_core=1 (J(Tk)>=Tcore); eta_core,min at Tk=20 and 30 K')
from scipy.optimize import brentq
for reg,tb in (('N',(2.635,2.822)),('S',(1.343,2.942))):
    for s in (0.10,0.15,0.20,0.25,0.34):
        dil=s**2/(s**2+BEAM**2); tc=[t/dil for t in tb]
        tkm=[brentq(lambda T:J(T)-t,2.8,5000) for t in tc]
        print(reg,'s=%.2f"'%s,'dil %.2f'%dil,'Tcore HCN/HCO+ %.1f/%.1f K'%tuple(tc),' Tk,min(eta=1) %.0f/%.0f K'%tuple(tkm),' eta_min@20K %.2f/%.2f  @30K %.2f/%.2f'%(tc[0]/J(20),tc[1]/J(20),tc[0]/J(30),tc[1]/J(30)))
print('\n== (2) physical quantities from Q20b posteriors (1z, Tbg 2.73, lq0)')
G={m:np.load(f'/workspace/out/q4c/grid4_{m}.npz') for m in ('hcn','hco+')}
TK,LN,LNC=G['hcn']['TK'],G['hcn']['LN'],G['hcn']['LNC']
def tau(m,lT,ln,lN,j=2):
    c=np.array([np.interp(lT,np.log10(TK),np.arange(len(TK))),np.interp(ln,LN,np.arange(len(LN))),np.interp(lN,LNC,np.arange(len(LNC)))])
    return map_coordinates(G[m]['tau'][...,j].astype(float),c,order=1,mode='nearest')
mu=2.8*1.6605e-24; pc=3.086e18; Msun=1.989e33
out=[]
for reg,gs in (('S',0.25),('N',0.27)):
  for ds in ('r8a','r8'):
    for tag,Xab in (('X=1e-8',1e-8),('X=1e-9',1e-9)):
        f=sorted(glob.glob(f'/workspace/out/q20/runs/lq0_1z_{reg}_{ds}_tbg2.73_s*.npz'))
        X=np.vstack([np.load(x,allow_pickle=True)['X'] for x in f]); lT,n,NH,NO,lfc=X.T
        # RADEX grid N is the column of a uniform sphere at dv=150 km/s; N_H2 = N(HCN)/X(HCN), X(HCO+)=X(HCN)*N_O/N_H unconstrained -> use HCN
        NH2=10**NH/Xab; L=NH2/10**n/pc; M=10**lfc*NH2*mu*1.133*(gs*PC*pc)**2/Msun; P=10**n*10**lT
        t32h=tau('hcn',lT,n,NH); t32o=tau('hco+',lT,n,NO); t43=tau('hcn',lT,n,NH,3)
        # H13CN(4-3) test: HCN/H13CN flux ratio (1-e^-t)/(1-e^-t/R) must exceed the limit (N>8.3, S>4; Q12) for 12C/13C=R
        lim=8.3 if reg=='N' else 4.0
        pH={R:float(np.mean((1-np.exp(-t43))/(1-np.exp(-t43/R))>lim)) for R in (30,50,100)}
        q=lambda a:np.nanpercentile(a,[16,50,84])
        row=dict(region=reg,data=ds,abund=tag,T=q(10**lT),n=q(n),logNHCN=q(NH),logNHCO=q(NO),logNH2=q(np.log10(NH2)),L_pc=q(L),Mcore_Msun_log=q(np.log10(M)),logP=q(np.log10(P)),eta=q(10**lfc),tau32HCN=q(t32h),tau32HCO=q(t32o),Xratio_HCN_HCO=q(NH-NO),tau43HCN=q(t43),P_H13CN_ok=pH)
        out.append(row); print(reg,ds,tag,{k:(np.round(v,2) if not isinstance(v,(str,dict)) else v) for k,v in row.items() if k not in('region','data','abund')})
pd.DataFrame(out).to_csv('/workspace/out/q20/q20g_phys.csv',index=False)
print('\n== (4) dust continuum cross-check N_H2 (tau=Tb/[J_258(Td)-J(Tbg)], kappa_1.2mm=1.0 cm2/g dust, GDR 100, mu=2.8)')
h=12.4
Jd=lambda T:h/np.expm1(h/T)-h/np.expm1(h/2.73)
for reg,Tb in (('N',0.18),('S',0.59)):
    for Td in (20,30,50,80):
        tau_=Tb/Jd(Td); print(reg,'Td',Td,'tau1.2mm %.3f'%tau_,'N_H2=%.1e cm-2 (beam-averaged, 0.24")'%(tau_*100/(1.0*mu)))
