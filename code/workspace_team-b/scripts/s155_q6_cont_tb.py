"""Q6: continuum Tb at 0.24" (257.76 GHz) and a dust-only order-of-magnitude surface density (upper bound if all continuum is dust). python3 -I s155_q6_cont_tb.py"""
import json,numpy as np
d=json.load(open('/workspace/out/q6/s150.json')); bm=d['cont']['beam']; nu=257.76
J=lambda T: (4.799e-11*nu*1e9)/np.expm1(4.799e-11*nu*1e9/T)
Jbg=J(2.73); res={}
for n,Td in (('S',50),('S',100),('N',50),('N',100),('SE',50)):
    f=d['cont']['fits'][n]['1.0']; A=f['Apk_mJy']; eA=f['eA']
    Tb=1.222e3*A/(nu**2*bm[0]*bm[1]); 
    dd=(f['dec_maj']*f['dec_min'])**.5; bg=(bm[0]*bm[1])**.5
    tau=-np.log(1-Tb/(J(Td)-Jbg)); 
    # kappa_gas(nu)=0.1 (nu/1000 GHz)^1.5 cm2/g(gas)  (Beckwith+1990 convention; factor ~3 uncertain)
    kap=0.1*(nu/1000.)**1.5; Sig_g=tau/kap*4785.   # Msun/pc2 beam-averaged (1 g/cm2=4785 Msun/pc2)
    res[f'{n}_Td{Td}']=dict(Apk_mJy=A,Tb_RJ=Tb,tau=tau,Sigma_beam_Msun_pc2=Sig_g,dec=dd,beam_geo=bg)
    print(n,Td,'Apk %.3f mJy Tb %.3f K  tau %.4f  Sigma(beam-avg, all dust) %.0f Msun/pc2  dec %.2f beam %.2f'%(A,Tb,tau,Sig_g,dd,bg))
json.dump(res,open('/workspace/out/q6/s155.json','w'),indent=1)
