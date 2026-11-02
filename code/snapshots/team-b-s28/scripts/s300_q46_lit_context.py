"""Q46: place NGC 3256 N and S HCN/HCO+ among published values.
Literature numbers are transcribed from lit.py full text (Nishimura+2024 Tables 1,2,4; Izumi+2016 Table 1;
Michiyama+2018 Table 4; Harada+2018 Table 13; Ueda+2021 text quoting Privon+2015). Our numbers: report.json
Q26 (3-2, Xb4b/Xb4d), Q29 (1-0, X145/X961), Xb47 4-3 line properties. All ratios are T_b (K km/s) ratios."""
import numpy as np
from scipy.stats import spearmanr
D=40.0; asec=D*1e6*4.848e-6  # pc per arcsec ~ 194
names=["IRAS17208-0014","IRAS09022-3615","IRAS13120-5453","IRASF14378-3651","IRASF05189-2524","IRASF10565+2448",
"IRAS19542+1110","ESO148-IG002","NGC6240","IRASF17138-1017","IRAS17578-0400","ESO173-G015","UGC11763","NGC3110",
"IC4734","NGC5135","ESO221-IG10","IC5179","UGC2982","NGC2369","ESO286-G035","ESO320-G030","NGC5734"]
DL=[183,253,135,295,176,184,264,185,106,77.3,60.0,32.7,265,74.3,67.8,52.6,59,47.3,70.6,44.6,73.2,36,60.2]
aAGN=[.05,.19,.04,.05,.60,.04,.08,.19,.11,.07,.03,.03,np.nan,.10,.07,.24,.07,.10,.11,.10,.06,.03,.14]
res_pc=[298,382,225,388,343,281,423,344,161,139,102,51,397,142,126,83,95,63,117,59,101,62,89]  # geometric-ish major axis
resolved=[1.738,.793,2.671,1.323,1.630,1.320,1.789,.445,.691,.629,1.340,1.132,.728,1.490,1.076,1.028,.771,.564,np.nan,1.489,.666,2.344,np.nan]
globl=[1.464,.555,1.706,1.140,1.534,.878,1.825,.387,.693,.587,1.283,.865,np.nan,1.044,1.063,1.100,.651,.842,np.nan,1.133,.587,2.254,np.nan]
R=np.array(resolved);G=np.array(globl);A=np.array(aAGN);P=np.array(res_pc)
ok=~np.isnan(R)
print("Nishimura+2024 HCN/HCO+ 3-2 (L' ratio) n_res=%d n_glob=%d"%(ok.sum(),(~np.isnan(G)).sum()))
print("resolved: median %.2f mean %.2f range %.2f-%.2f; global: median %.2f mean %.2f"%(np.nanmedian(R),np.nanmean(R),np.nanmin(R),np.nanmax(R),np.nanmedian(G),np.nanmean(G)))
def frac(a,x): a=a[~np.isnan(a)]; return (a<x).mean(),(a<=x).mean()
ours={"N 55pc (Xb4b)":0.99,"N 250pc":0.69,"S 55pc":0.44,"S 250pc":0.37,"N+S+Bridge 6.5\" (1.3 kpc)":0.44}
for k,v in ours.items(): print("  ours %-28s %.2f  fraction of Nishimura resolved below: %.2f; global below: %.2f"%(k,v,frac(R,v)[0],frac(G,v)[0]))
# matched-resolution subsets
for lo,hi in[(0,110),(110,450)]:
    m=ok&(P>=lo)&(P<hi); print("resolution %d-%d pc: n=%d resolved ratios median %.2f [%.2f-%.2f]"%(lo,hi,m.sum(),np.median(R[m]),R[m].min(),R[m].max()))
m=ok&(P<110)
for n,p,r,a in zip(np.array(names)[m],P[m],R[m],A[m]): print("   %-16s %3d pc  R=%.2f  aAGN=%.2f"%(n,p,r,a))
# AGN fraction split
for lab,mm in[("aAGN<0.10",ok&(A<0.10)),("0.10<=aAGN<0.25",ok&(A>=0.10)&(A<0.25)),("aAGN>=0.25",ok&(A>=0.25))]:
    print("  %-16s n=%2d resolved median %.2f  range %.2f-%.2f"%(lab,mm.sum(),np.median(R[mm]),R[mm].min(),R[mm].max()))
mm=ok&~np.isnan(A)
print("Spearman resolved vs aAGN: rho=%.2f p=%.2f (n=%d); global: rho=%.2f p=%.2f"%(*spearmanr(A[mm],R[mm]),mm.sum(),*spearmanr(A[mm&~np.isnan(G)],G[mm&~np.isnan(G)])))
# resolved/global (central enhancement) vs ours
rg=R/G; print("Nishimura resolved/global: median %.2f; >1 in %d/%d"%(np.nanmedian(rg),np.nansum(rg>1),np.sum(~np.isnan(rg))))
print("ours N core(55pc)/N(6.5\" blend)=%.2f ; S 0.44/0.44=1.0"%(0.99/0.44))
# Izumi 2016 J=4-3
iz={"AGN(Seyfert) hi-res":[1.53,2.84,3.19,3.14,2.58,2.01,1.11],"SB hi-res (M82,NGC253,1097ring,1614,7469 B/C/D)":[.41,1.03,.82,.24,.75,.45,.48],
    "LMC N113,N159":[.21,.11],"low-res buried AGN (NGC4418,IRAS12127)":[1.64,1.58],"low-res AGN (1068,1365,4945)":[1.85,.59,.78],
    "low-res SB (NGC3256,3628,7552,IRAS13242)":[.39,.36,.54,1.23]}
print("Izumi+2016 HCN4-3/HCO+4-3:")
for k,v in iz.items(): print("  %-48s median %.2f range %.2f-%.2f"%(k,np.median(v),min(v),max(v)))
n43=14.679/28.725*(356.734/354.505)**2; s43=5.948/20.969*(356.734/354.505)**2
print("ours 4-3 (Xb47 4.14x2.70\", peak pixel, +-350, Tb ratio): N %.2f S %.2f; N+S combined (Q38r2, 10\") 0.38; Izumi NGC3256 (APEX) 0.39+-0.11"%(n43,s43))
print("1-0: Harada 1.7\" N 0.96 S 0.645; Michiyama 2.1\" N 0.94+-0.13 S 0.60+-0.08; ours 1.5\" N 0.94 S 0.58, 2.2\" N 0.84 S 0.57 ; GOALS(Privon via Ueda) mean 0.92-1.02, 0.77 w/o AGN; merger remnants 0.88+-0.08")
print("scale: 1\"=%.0f pc ; 0.26\"=%.0f pc ; 1.7\"=%.0f pc ; 4.1x2.7\"=%.0fx%.0f pc"%(asec,0.26*asec,1.7*asec,4.1*asec,2.7*asec))
