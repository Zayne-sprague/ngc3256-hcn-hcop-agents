# Q9: E-CH3OH column from the 254.015 GHz line (Xb4f smoothed to the Xb55 9.36x6.28" beam, beam-averaged) for T=6,10,18 K (Trot range from s120c).
# N_u=8 pi k nu^2 W/(h c^3 A); N(E)=N_u/g_u*Q_E(T)*exp(Eu/T); Q_E from the rigid-rotor-like T^1.5 scaling Q_E~0.5*0.235*T^1.5 (A+E total Q=0.235 T^1.5, half each; +-30%). Beam-averaged, no beam-filling correction.
import json,numpy as np
k=1.380649e-16;h=6.62607e-27;c=2.998e10
S=json.load(open('/workspace/out/s21_ch3oh254.json'))['smoothed']
nu=254.015377e9; A=10**-4.72074; gu=5; Eu=20.09
th=9.36*6.28
out={}
for reg in ('S','Bridge'):
    I=S[reg]['I254_55']; Wk=1.222e6*I/(251.68**2*th)      # K km/s, observed-frame GHz
    Nu=8*np.pi*k*nu**2*(Wk*1e5)/(h*c**3*A)
    out[reg]={'W_K_kms':Wk,'Nu':Nu}
    for T in (6,10,18):
        QE=0.5*0.235*T**1.5; out[reg][f'NE_T{T}']=Nu/gu*QE*np.exp(Eu/T)
    print(reg,{kk:float('%.3g'%v) for kk,v in out[reg].items()})
json.dump(out,open('/workspace/out/q9/s120e.json','w'),indent=1)
