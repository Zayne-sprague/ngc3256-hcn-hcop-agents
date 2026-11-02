# Q13h: LTE column of HDO 2(2,0)-3(1,3) (266.16107 GHz; Splatalogue/JPL: Eu=157.21 K, log10 A=-4.756 -> A=1.75e-5 s^-1, S mu^2=0.399 D^2, g_u=5)
# needed if the whole E (0.115+-0.027 Jy/beam km/s, 1.05" beam, r=1" aperture; Q13a) were HDO. Optically thin LTE, CMB background correction, source fills the 1.05" beam.
# Compared with N(HCO+) from the parent-gas window (HCO+ 2685-2885, same aperture, 0.502 Jy/beam km/s) -> N(H2) via X(HCO+)=1e-9..1e-8; HDO/H2 required.
import numpy as np
h=6.62607e-27;k=1.380649e-16;c=2.99792458e10
Q300=588.   # Q_rot(HDO,300 K) from the asymmetric-rotor high-T formula (A,B,C = 23.1,9.1,6.4 cm-1), no spin; +-10%
def Qh(T): return Q300*(T/300.)**1.5
def Qhco(T): return k*T/(h*44.594e9)+1/3.
def Jy2K(nuG,th): return 1.222e6/(nuG**2*th**2)
def N_tot(Tdv,T,nu,A,gu,Eu,Qf):
    x=h*nu/k; f=1-(np.exp(x/T)-1)/(np.exp(x/2.73)-1)          # CMB-background correction, thin line
    Nu=8*np.pi*k*nu**2*Tdv*1e5/(h*c**3*A)/f
    return Nu*Qf(T)/gu*np.exp(Eu/T)
E=0.115; eE=0.027; th=1.05
Tdv=E*Jy2K(266.161,th); eT=eE*Jy2K(266.161,th); Tp=0.502*Jy2K(267.5576,th)
print(f'E={E}+-{eE} Jy/beam km/s -> {Tdv:.2f}+-{eT:.2f} K km/s; HCO+ parent window {Tp:.2f} K km/s (beam-filling)')
print(' T[K]  N_HDO[cm-2] (1sig range)   N_HCO+thin  N_HDO/N_HCO+  N_H2(X_HCO+=1e-8..1e-9)   X_HDO=N_HDO/N_H2   [beam-filling; multiply by (beam/size)^2 if smaller]')
for T in (50,100,150,200,300):
    Nh=N_tot(Tdv,T,266.16107e9,1.753e-5,5,157.21,Qh); Nl=N_tot(Tdv-eT,T,266.16107e9,1.753e-5,5,157.21,Qh); Nu=N_tot(Tdv+eT,T,266.16107e9,1.753e-5,5,157.21,Qh)
    Nc=N_tot(Tp,T,267.5576e9,1.447e-3,7,25.69,Qhco)
    print(f'{T:4d}  {Nh:.2e} ({Nl:.1e}-{Nu:.1e})  {Nc:.2e}  {Nh/Nc:6.2f}   {Nc/1e-8:.1e}..{Nc/1e-9:.1e}   {Nh/(Nc/1e-8):.1e}..{Nh/(Nc/1e-9):.1e}')
print('\nForward: HDO line for X_HDO=N_HDO/N_H2 = 1e-9, 1e-8 with N_H2 = N_HCO+/X_HCO+ (X_HCO+=1e-9..1e-8), T=150 K, beam-filling:')
T=150; Nc=N_tot(Tp,T,267.5576e9,1.447e-3,7,25.69,Qhco); N1=N_tot(1.0,T,266.16107e9,1.753e-5,5,157.21,Qh)  # column per 1 K km/s
for XH in (1e-8,1e-9):
    for XD in (1e-9,1e-8):
        NH2=Nc/XH; Nd=XD*NH2; ITdv=Nd/N1; print(f' X_HCO+={XH:.0e} X_HDO={XD:.0e}: N_H2={NH2:.1e} N_HDO={Nd:.1e} -> HDO {ITdv:.4f} K km/s = {ITdv/Jy2K(266.161,th)*1e3:.2f} mJy/beam km/s (vs E 115)')
# optical depth check of the HDO line at the required column (T=150, dv=150 km/s FWHM)
T=150; N=N_tot(Tdv,T,266.16107e9,1.753e-5,5,157.21,Qh)
x=h*266.16107e9/k; Nu=N*5/Qh(T)*np.exp(-157.21/T); dv=150e5/1.0645*1.0645  # FWHM cm/s
tau=c**3*1.753e-5*Nu/(8*np.pi*(266.16107e9)**2*dv)*1.0645*(np.exp(x/T)-1)/1*1 if False else None
