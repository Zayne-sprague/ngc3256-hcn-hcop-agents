# Q2 reviewer fix 2: bound on the dust radiation field from the observed 0.25" 258 GHz continuum (Xb4b spw25_27_29_31.cont; out/continuum_table.csv, pbcor peak, Gaussian beam).
# RJ brightness temperature Tb = 1.222e6 S[Jy/beam] / (nu[GHz]^2 theta_maj theta_min["]) K; extrapolated to the HCN/HCO+ 4-3 / 3-2 frequencies with S~nu^3.5 (optically thin dust, beta~1.5) i.e. Tb~nu^1.5.
# The extra field seen by the gas is <= this (beam-averaged) value if the dust fills the beam; a smaller dust source (filling f) could give up to Tb/f, capped by Td (a thick dust screen at Td~50-100 K would have RJ Tb(258 GHz) ~ 44-90 K).
import pandas as pd, numpy as np
c=pd.read_csv('/workspace/out/continuum_table.csv'); c=c[c.file.str.contains('Xb4b.spw25_27_29_31.cont')&c.src.isin(['N','S'])]
rows=[]
for _,r in c.iterrows():
    bmaj,bmin=[float(x) for x in r.beam.split('x')]; nu=r.nu_GHz; Tb=1.222e6*r.peak_Jybeam/(nu**2*bmaj*bmin)
    # Planck-inverse brightness temperature for completeness
    h=4.799e-11; I=Tb*0+0; x=h*nu*1e9; Tp=x/np.log(1+x/(Tb))  # T_planck: from RJ Tb: B=2k nu^2 Tb/c^2 => Tp = x/ln(1+x/Tb)
    for lab,f in (('3-2 (265 GHz HCN... 1.4x)',(265/nu)**1.5),('HCN/HCO+ 3-2 band-6 (~265)',(265/nu)**1.5),('HCN/HCO+ 4-3 (355 GHz)',(355/nu)**1.5)):
        pass
    rows.append(dict(src=r.src,nu_GHz=nu,beam=r.beam,peak_mJy=r.peak_Jybeam*1e3,Tb_RJ_K=Tb,Tb_Planck_K=Tp,Tb_RJ_at_265GHz_K=Tb*(265/nu)**1.5,Tb_RJ_at_355GHz_K=Tb*(355/nu)**1.5,
                     Tb_RJ_355_if_filling_0p1=Tb*(355/nu)**1.5/0.1,Tbg_equiv_355_beam_avg=2.73+Tb*(355/nu)**1.5))
R=pd.DataFrame(rows); R.round(3).to_csv('/workspace/out/q2/q2_dust_field.csv',index=False); print(R.round(3).to_string())
