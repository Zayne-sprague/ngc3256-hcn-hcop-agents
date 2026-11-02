# Q20d: (a) eta lower bound from 0.24" T_b (T_ex<=T_k); (b) core share of the 10" aperture flux (is the 10" ladder blind to a cool core?);
# (c) envelope ladder after subtracting a cool-core ladder.  Inputs: out/q4b/dilution.csv, out/q1 deblend ladders, team-b p0028 10" intensities.
import numpy as np, pandas as pd
D=pd.read_csv('/workspace/out/q4b/dilution.csv'); L=pd.read_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv')
hv=12.8  # h nu/k for 3-2 (K), HCN 265.9 / HCO+ 267.6 GHz
def J(T,Tbg=2.73): return hv/np.expm1(hv/T)-hv/np.expm1(hv/Tbg)
BEAM=0.24; FW={'N':150.,'S':117.}   # FWHM km/s (q5: sigma_obs 49.6 S, 56-70 N)
print("(a) eta_min = Tb/J(Tk) for thick lines (Tex<=Tk). Tint = pedestal-inclusive deconvolved (dilution.csv); Tcore = Tb_025/dil_gauss for Q5 Gaussian cores")
GS={'N':0.27,'S':0.25}   # Q5 Gaussian FWHM (arcsec)
rows=[]
for _,r in D[D.region.isin(['N','S'])].iterrows():
    dg=GS[r.region]**2/(GS[r.region]**2+BEAM**2); Tcore=r.Tb_025/dg
    for Tk in (20,25,30,50,80,150):
        rows.append(dict(region=r.region,line=r.line,Tb_025=r.Tb_025,Tint=r.Tint_K,Tcore=Tcore,Tk=Tk,eta_min_Tb=r.Tb_025/J(Tk),eta_min_Tint=r.Tint_K/J(Tk),eta_min_Tcore=Tcore/J(Tk)))
A=pd.DataFrame(rows); A.to_csv('/workspace/out/q20/q20d_eta_min.csv',index=False); print(A.round(3).to_string())
print("\n(b) core share of 10\" 3-2 flux; I_core(10\")=Tint*s^2/(s^2+10^2)*1.0645*FWHM (pedestal-inclusive s = upper bound); Gaussian-core variant uses Tcore, s=Q5")
I10={'N':{'HCN':(12.38,4.81,1.62),'HCO+':(20.27,10.04,4.51)},'S':{'HCN':(11.66,4.26,1.35),'HCO+':(20.89,9.95,4.28)}}
sh=[]
for _,r in D[D.region.isin(['N','S'])].iterrows():
    sp=r.line.split()[0]; dil=r.s_fwhm_arcsec**2/(r.s_fwhm_arcsec**2+100); Ic=r.Tint_K*dil*1.0645*FW[r.region]
    dg=GS[r.region]**2/(GS[r.region]**2+BEAM**2); Tcore=r.Tb_025/dg; Icg=Tcore*GS[r.region]**2/(GS[r.region]**2+100)*1.0645*FW[r.region]
    tot=I10[r.region][sp][1]; sh.append(dict(region=r.region,sp=sp,I10=tot,Icore_ped=Ic,frac_ped=Ic/tot,Icore_gauss=Icg,frac_gauss=Icg/tot))
S=pd.DataFrame(sh); S.to_csv('/workspace/out/q20/q20d_core_share10.csv',index=False); print(S.round(3).to_string())
print("\n(c) envelope ladder = 10\" minus cool core (core 3-2 share f=frac_ped; core ladder from 1.7\" deblend)")
for reg in ('N','S'):
    for sp in ('HCN','HCO+'):
        l=lambda q: L[(L.comp==reg)&(L.species==sp)&(L.quantity==q)].iloc[0].value
        r32=l('3-2 [Xb4f]/2-1 (Tb)'); r43=l('4-3/3-2 [Xb4f 3-2] (Tb)')
        i21,i32,i43=I10[reg][sp]; f=float(S[(S.region==reg)&(S.sp==sp)].frac_ped.iloc[0])
        c32=f*i32; c21=c32/r32; c43=c32*r43
        print(f"{reg} {sp}: 10\" 3-2/2-1={i32/i21:.3f} 4-3/3-2={i43/i32:.3f} | core share {f:.2f}, core ladder {r32:.2f},{r43:.2f} | envelope {(i32-c32)/(i21-c21):.3f}, {(i43-c43)/(i32-c32):.3f}")
