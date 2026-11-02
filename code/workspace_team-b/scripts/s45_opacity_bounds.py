# Q4e: turn the H13CN limits (s41) and RADEX table (s44) into opacity bounds.
# (1) Beam-averaged tau(HCN 4-3) upper limit from I(HCN)/I(H13CN) > L for 12C/13C = X (same Tex, same filling): L = (1-exp(-tau12))/(1-exp(-tau12/X)).
# (2) What fraction of the HCN 4-3 flux in the Xb47 N beam-aperture could come from an optically thick compact core (ratio r_c) given the beam-averaged limit L?
#     Core flux fraction f_core from the 0.25" HCN 3-2 annuli (s43; r<0.5" vs the r<1.7" Xb47 aperture, flux = mean * npix / beam_area_pix) assuming the 4-3 flux fraction equals the 3-2 one.
# (3) RADEX regime separation of R_K(J): predicted R43/R32 for thin (tau32<1) and thick (tau32>10) cells with R_K(3-2)=0.85-1.1 (s44 table) vs observed R(J) ladder.
import numpy as np, pandas as pd, json
from scipy.optimize import brentq
def tau_ul(L,X):
    f=lambda t:(1-np.exp(-t))/(1-np.exp(-t/X))-L
    return brentq(f,1e-3,1e4) if f(1e-3)>0 else np.nan   # ratio->X as tau->0, so L<X required; f(1e-3)=X-L>0
out={}
print('(1) tau(HCN 4-3) upper limit from HCN/H13CN(4-3) > L')
tab=[]
for L in (4,6.7,8.3,9.5,12):
    for X in (30,50,100):
        tab.append(dict(L=L,X=X,tau12_max=tau_ul(L,X) if L<X else np.nan))
t=pd.DataFrame(tab).pivot(index='L',columns='X',values='tau12_max'); print(t.round(1).to_string()); out['tau_ul']=t.round(2).to_dict()
# (2) core fraction
df=pd.read_csv('/workspace/out/s43_profiles.csv'); n=df[df.region=='N'].reset_index(drop=True)
ps=0.047; bm=1.1331*0.2456*0.2373/ps**2
n['flux']=n.I_HCN*n.npix/bm
tot17=n.flux[:3].sum()+n.flux[3]*(1.7**2-1.0**2)/(2.0**2-1.0**2)    # annulus 1-2" scaled to 1-1.7" by area (flux per area const)
f05=n.flux[:2].sum()/tot17; f025=n.flux[0]/tot17
print('\n(2) N: HCN 3-2 flux in r<0.25, 0.25-0.5, 0.5-1, 1-2" (Jy km/s):',n.flux.round(2).tolist(),' total r<1.7":',round(tot17,2),' core fractions f(r<0.25)=%.2f f(r<0.5)=%.2f'%(f025,f05))
out['f_core']=dict(f025=f025,f05=f05)
print('  max core flux fraction f allowed by beam-averaged ratio L for a core with HCN/H13CN = r_c and envelope at X:')
rows=[]
for L in (8.3,9.5):
    for X in (30,50,100):
        for rc in (1.5,3,5,8):
            # 1/L >= f/rc + (1-f)/X  -> f <= (1/L-1/X)/(1/rc-1/X)
            f=(1/L-1/X)/(1/rc-1/X); rows.append(dict(L=L,X_env=X,r_core=rc,f_max=min(f,1.0)))
fm=pd.DataFrame(rows); print(fm.pivot_table(index=['L','X_env'],columns='r_core',values='f_max').round(2).to_string()); out['fmax']=fm.round(3).to_dict('records')
# (3)
r0=pd.read_csv('/workspace/out/s44_radex_opacity.csv'); r0['L43_32']=r0.R43/r0.R32; r0['R32_21']=r0.R32/r0.R21
r0['regime']=np.where(r0.tau32_HCN<1,'thin',np.where(r0.tau32_HCN>10,'thick','mid'))
for lab,lo,hi,obs in (('N core-like R_K(3-2)=0.85-1.10',0.85,1.10,'no 4-3/2-1 at the core beam'),('blended/outer R_K(3-2)=0.35-0.70',0.35,0.70,'observed 6.5" N: R43/R32=0.78 (0.43/0.55); 10" midpoint R43/R32=0.75, R32/R21=0.76')):
    r=r0[(r0.R32>lo)&(r0.R32<hi)]
    print('\n(3) RADEX cells,',lab,'; observed:',obs); print(r.groupby('regime')[['L43_32','R32_21','tau32_HCN']].agg(['min','median','max','count']).round(2).to_string())
json.dump(out,open('/workspace/out/s45_opacity_bounds.json','w'),indent=1,default=float)
# (4) optically thin LTE column density limits from the beam-averaged H13CN 4-3 3-sigma limit (N aperture, Xb47 beam 4.26x2.77"): N_u = 8 pi k nu^2 W / (h c^3 A); N = N_u Q / (g_u exp(-Eu/Tex)).
h=6.62607e-27;k=1.380649e-16;c=2.99792458e10
nu=345.3398e9; A=1.9e-3; Eu=41.4; gu=9.; B=43.17e9    # H13CN 4-3: A_ul (LAMDA h13cn@xpol; check) ; B rotational constant H13CN
S_HCN=10.923   # Jy/beam km/s, N beam aperture (s41 I_HCN 10923 mJy km/s)
sky=342.14e9; th=4.263*2.772
W_HCN=1.222e3*S_HCN*1e3/((sky/1e9)**2*th)          # K km/s (RJ, 1.222e3 S[mJy]/(nu^2 th1 th2))
print('\n(4) HCN 4-3 beam-aperture N: I=%.2f Jy/beam km/s = %.2f K km/s (observed-frame RJ)'%(S_HCN,W_HCN))
res4={}
for al in (0.105,0.12):
    W=al*W_HCN
    for Tex in (20,50,100):
        Nu=8*np.pi*k*nu**2*W*1e5/(h*c**3*A); Q=k*Tex/(h*B)+1/3.; N13=Nu*Q/(gu*np.exp(-Eu/Tex))
        print(f'  alpha<{al}: W(H13CN)<{W:.2f} K km/s, Tex={Tex} K: N(H13CN)<{N13:.2e} cm-2 ; N(HCN)<X*N: X=30 {30*N13:.1e}, 50 {50*N13:.1e}, 100 {100*N13:.1e}')
        res4[f'{al}_{Tex}']=N13
json.dump(out|{'N13_thin':res4,'W_HCN_Kkms':W_HCN},open('/workspace/out/s45_opacity_bounds.json','w'),indent=1,default=float)
