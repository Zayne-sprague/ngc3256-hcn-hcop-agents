# Q38: (a) MRS of 12m continuum, (b) global FIR-SED dust prediction at 258/346/171 GHz and 10.5" Tb for disks of 10-30" FWHM,
# (c) missing-fraction vs archive 10.5" matched Tb (q36_cont), (d) J_bg(nu) upper bound for envelope + equivalent Tbg.
import numpy as np, json, sys
sys.path.insert(0,'/workspace/scripts')
H=6.62607e-27;KB=1.380649e-16;C=2.998e10
def B(nu,T): return 2*H*nu**3/C**2/np.expm1(H*nu/(KB*T))
def RJ(nu,T): return (H*nu/KB)/np.expm1(H*nu/(KB*T))
def Tbg_eq(nu,J): return (H*nu/KB)/np.log1p((H*nu/KB)/J)
L=[]
# (a) MRS: theta_MRS ~ 0.983 lambda/L5 (ALMA Tech Handbook), L5 = 5th-percentile baseline. Beam -> config via beam~ 0.574 lambda/Lmax-ish; use beam to infer config.
# Xb4f beam 5.9x3.6 @258 GHz -> C43-1/2 (L5~ 15-20 m... handbook C43-1: MRS 28.5" @100GHz scales: 11" @ 258; C43-2: 26.6"@100 -> 10.3")
# Handbook Cycle 5 values at 100 GHz: C43-1 MRS 28.5", C43-2 26.6", C43-3 16.2", C43-4 14.4"; scale 1/nu. Beam at 100GHz: C43-1 3.4",C43-2 2.3"(handbook: 3.38,2.30?)
# Xb4f beam geometric mean 4.6" at 258 -> 11.9" at 100 GHz => ~C43-1 (beam 3.4-4"?) conservative; Xb47 beam 3.4" at 346 -> 11.7"@100 same; Xb55 beam 7.3" at 170 -> 12.4"@100
for name,nu,beam,mrs100 in (('Xb55 B5',170.7,np.sqrt(9.0*6.1),(28.5,16.2)),('Xb4f B6',257.7,np.sqrt(5.9*3.6),(28.5,16.2)),('Xb47 B7',346.5,np.sqrt(4.15*2.7),(28.5,16.2))):
    L.append(f'{name}: nu={nu} GHz, beam(gm)={beam:.1f}" -> MRS ~ {mrs100[0]*100/nu:.1f}" (C43-1-like) .. {mrs100[1]*100/nu:.1f}" (C43-3-like); handbook 0.983 lambda/L5')
# (b) global FIR dust: S160=93.48 Jy,S100=145.4,S70=120.3 (Michiyama+20). Fit Td (Q36: 29-34 K); scale to mm with beta
res={}
nu160=C/160e-4
for Td in (25,30,34):
  for beta in (1.5,2.0):
    for nu in (170.7e9,257.7e9,346.5e9):
        S=93.48*(nu/nu160)**beta*B(nu,Td)/B(nu160,Td)   # Jy, optically thin extrapolation from 160um
        res[(Td,beta,nu)]=S
    L.append(f'global S_dust(Td={Td},beta={beta}): '+', '.join(f'{res[(Td,beta,n)]*1e3:.0f} mJy @{n/1e9:.0f}' for n in (170.7e9,257.7e9,346.5e9)))
# tau(258 GHz) for 10-30" gaussian disk: tau = S/(Omega B(Td))
beam_area=lambda fwhm: np.pi*(fwhm)**2/(4*np.log(2))*(np.pi/180/3600)**2
arc=np.pi/180/3600
L.append('\nDisk model: gaussian FWHM F carrying total S (Td=30,beta=1.8 central), Tb_RJ at beam centre of 10.5" beam, and recovered fraction if flux on scales >MRS is lost')
Sg=lambda nu,Td=30,beta=1.8: 93.48*(nu/nu160)**beta*B(nu,Td)/B(nu160,Td)
out={}
for nu in (170.7e9,257.7e9,346.5e9):
  for F in (10,15,20,30):
    S=Sg(nu); Om=beam_area(F); Ipk=S*1e-23/Om*4*np.log(2)/(4*np.log(2))  # peak intensity of gaussian = S/Omega_gauss
    Ipk=S*1e-23/Om
    # convolve with 10.5" beam: peak decreases by F^2/(F^2+b^2); Tb_RJ=I c^2/(2 k nu^2)
    b=10.5; f=F**2/(F**2+b**2)
    TbRJ=Ipk*f*C**2/(2*KB*nu**2); tau=Ipk*1e0/B(nu,30)
    out[(nu,F)]=TbRJ
    L.append(f'  nu={nu/1e9:.0f} F={F}": S_tot={S*1e3:.0f} mJy  Tb_RJ(10.5" beam)={TbRJ*1e3:.1f} mK  tau_pk={tau:.1e}  J_CMB={RJ(nu,2.725)*1e3:.0f} mK  -> dust/J_CMB={TbRJ/RJ(nu,2.725):.2f}')
# gas along line of sight: J_bg from a thin disk of tau seen by gas embedded at mid-plane: J = 0.5*I_dust,los(unconvolved, intrinsic), i.e. use unconvolved peak
L.append('\nIntrinsic (unresolved by 10.5" beam) disk mean intensity at centre; J_bg(env) = J_CMB + eta*Tb_RJ_intrinsic (eta=0.5 embedded, 1 behind)')
Jub={}
for nu in (88.63e9,177.26e9,265.89e9,354.51e9):
    S=Sg(nu*1.0, 30,1.8); 
    row=[]
    for F in (10,20,30):
        Ipk=S*1e-23/beam_area(F); Tb=Ipk*C**2/(2*KB*nu**2); row.append(Tb)
    Jub[nu]=row
    L.append(f'  nu={nu/1e9:.1f}: Tb_RJ,pk(F=10/20/30")='+'/'.join(f'{t*1e3:.1f}' for t in row)+' mK (dust only; ff would add at 88 GHz)')
json.dump({str(k):v for k,v in Jub.items()},open('/workspace/out/q38/jub.json','w'))
print('\n'.join(L)); open('/workspace/out/q38/q38_mrs.log','w').write('\n'.join(L))
