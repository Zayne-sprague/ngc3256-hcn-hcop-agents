"""Q51: is NGC 3256 N CON-like?  (1) our HCN-vib (v2=1f) limits converted to Falstad+2021 quantities
(L_HCN-vib in Lsun, Sigma_HCN-vib, CON criterion Sigma>1 Lsun/pc2) and to vib/HCN L' ratios of published HCN-vib detections;
(2) 1.1 mm continuum brightness of N, S vs Falstad Table 3 CON / non-CON nuclei; (3) Nishimura+2024 HCN/HCO+ split by Paper I (Falstad+2021) CON flag.
Literature numbers transcribed from lit.py full text: Falstad+2021 (2021A&A...649A.105F) Table 3, Nishimura+2024 (2024A&A...686A..48N) Tables 3,4.
Our limits: out/s42_hcn_vib.csv (scripts/s42_hcn_vib.py)."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
D_L=41.0  # Mpc (z=0.0094); Falstad Table 1 DL for NGC 3256 not extracted; 194 pc/arcsec as project
PC=194.; z=0.0094
def Lsun(S_Jykms,nu_rest=267.199283,DL=D_L): return 1.04e-3*S_Jykms*nu_rest/(1+z)*DL**2   # Falstad eq.1
def Lprime(Ls,nu): return 3.125e10*Ls/nu**3     # K km/s pc2 per Lsun (z-independent form derived from Falstad eq.1 + Solomon&Vanden Bout)
print("conversion 1e3 Lsun at 267.199 GHz = %.3g K km/s pc2"%(Lprime(1e3,267.199283)))
# ---- Falstad Table 3 (LIRGs incl. N,S) : name, L_vib(1e3 Lsun; negative=3sigma UL), Sigma, S_cont(mJy 1.1mm), size a,b (mas)
F=[("NGC1614",-.425,-.002,16.6,1930,1550),("NGC7469",-.121,-.02,3.76,342,229),("NGC3256 S",-.065,-.004,9.28,1247,545),
("NGC3256 N",-.065,-.001,17.6,1790,1290),("IRASF17138",-.247,-.003,3.75,1217,687),("IRAS17578",15.2,4.76,54.57,227,212),
("NGC7130",-.203,-.005,13.03,745,583),("ESO173",1.2,.19,50.6,959,333),("NGC3110",-.157,-.004,1.17,987,358),("IC4734",-.18,-.004,10.64,767,621),
("Zw049.057",11.3,3.37,37.2,294,171),("NGC5135",-.126,-.01,3.07,669,368),("ESO221",-.143,-.004,1.98,1130,472),("IC5179",-.055,-.01,1.93,408,323),
("UGC2982",-.124,-.003,.49,800,457),("ESO286",-.15,-.004,2.15,1213,304),("NGC4418",4.9,9.12,99.3,190,138),("NGC2369",-.075,-.003,11.82,1387,495),
("NGC5734",-.115,-.02,.38,638,166),("ESO320",2.0,1.4,33.0,282,208),
("IRAS17208",62.2,2.15,37.9,273,171),("ArpII W",61.7,11.96,189.3,220,190),("ArpII E",10.1,.95,80.5,380,230),("IRASF14378E",3.9,.019,3.61,350,349),("IRASF22491",27.2,4.84,4.41,66,43),
("IRASF12112NE",22.1,.23,9.29,201,261)]
rows=[]
nu_obs=265.0
for n,L,S,Sc,a,b in F:
    Tb=1.222e3*Sc/(nu_obs**2*(a/1e3)*(b/1e3))      # RJ T_b of the Gaussian-deconvolved continuum source at ~265 GHz
    rows.append((n,L,S,Sc,a,b,Tb,abs(S)>=1 and S>0 and L>0))
df=pd.DataFrame(rows,columns=["src","Lvib_k","Sig","Scont","a","b","Tb_K","CON"])
df["detected"]=df.Lvib_k>0
print(df.round(3).to_string())
con=df[df.CON]; ncon=df[~df.CON&~df.src.str.startswith("NGC3256")]
print("\nFalstad continuum T_b(RJ, deconvolved Gaussian) at ~265 GHz [K]: CON median %.2f (range %.2f-%.2f, n=%d); non-CON median %.2f (%.2f-%.2f, n=%d)"%(
 con.Tb_K.median(),con.Tb_K.min(),con.Tb_K.max(),len(con),ncon.Tb_K.median(),ncon.Tb_K.min(),ncon.Tb_K.max(),len(ncon)))
for n in ("NGC3256 N","NGC3256 S"):
    r=df[df.src==n].iloc[0]; print("  %s: S=%.1f mJy, %.0fx%.0f pc, T_b=%.3f K, frac of non-CON below %.2f, of CON below %.2f"%(n,r.Scont,r.a*PC/1e3,r.b*PC/1e3,r.Tb_K,(ncon.Tb_K<r.Tb_K).mean(),(con.Tb_K<r.Tb_K).mean()))
print("  Falstad N: L_vib<65 Lsun, Sigma<%.4f (check: 65/(pi/4*%.0f*%.0f pc2)=%.4f)"%(.001,1790*PC/1e3,1290*PC/1e3,65/(np.pi/4*1790*PC/1e3*1290*PC/1e3)))
# ---- our HCN-vib limits -> L_vib, Sigma at beam scale (assume source fills the beam, Gaussian beam area pi/4 a b), 3 sigma
v=pd.read_csv('/workspace/out/s42_hcn_vib.csv'); v=v[(v.kind=='pix')&(v.chi2red<1.5)]
beams={'J=3-2 Xb4b':(0.245,0.235,267.199283),'J=3-2 Xb4d':(0.95,0.75,267.199283),'J=3-2 Xb4f':(5.7,3.5,267.199283),'J=4-3 Xb47':(4.2,2.8,356.256)}
print("\nOur 3-sigma HCN-vib limits (peak pixel; widest of 4 wing models):")
out=[]
for setup,(a,b,nu) in beams.items():
    for reg in ('N','S'):
        q=v[(v.setup==setup)&(v.region==reg)]
        if q.empty: continue
        al=q.alpha_3s_ul.max(); I=q.I_HCN_mJykms.iloc[0]/1e3     # Jy/beam km/s
        Lv=Lsun(al*I,nu_rest=nu); A=np.pi/4*a*b*PC**2; Sg=Lv/A
        print("  %-11s %s beam %.2fx%.2f\" (%.0fx%.0f pc): alpha<%.3f  I_HCN=%.2f Jy/beam km/s  L_vib(beam)<%.0f Lsun  Sigma_vib<%.3f Lsun/pc2 (CON: >1; Falstad N/S: <0.001/<0.004)"%(setup,reg,a,b,a*PC,b*PC,al,I,Lv,Sg))
        out.append((setup,reg,al,Lv,Sg))
# ---- literature vib/HCN ratios in detections (K km/s ratio = L' ratio); Nishimura Table 3 HCN L' (1e7 K km/s pc2): resolved (beam), global (10x beam)
Lp=lambda L,nu=267.199283: Lprime(L*1e3,nu)
det={"IRAS17208 (CON)":(62.2,20.5,94.5),"IRAS17578 (CON)":(15.2,3.613,13.82),"ESO320 (CON)":(2.0,.798,5.67),"ESO173 (vib det, not CON)":(1.2,.297,3.38),"IRASF14378 (vib det, not CON)":(3.9,6.25,29.4)}
print("\nPublished HCN-vib / HCN 3-2 (L' ratio) = Falstad L_vib (whole nucleus) / Nishimura HCN L' (resolved = 1 beam, global = 10 beams):")
for k,(L,r,g) in det.items(): print("  %-32s vs resolved %.2f   vs global %.3f"%(k,Lp(L)/(r*1e7),Lp(L)/(g*1e7)))
print("  (F14378 E vib vs HCN of whole system: lower bound on dilution; Falstad L_vib sums nucleus E only)")
print("  ours N: alpha<0.025-0.05 (55 pc), <0.02 (250 pc, 4-3 800 pc), <0.01 (1.1 kpc); S: <0.085-0.10, <0.03-0.05, <0.015-0.02")
# ---- Nishimura split by Paper I flag
names=["IRAS17208","IRAS09022","IRAS13120","IRASF14378","IRASF05189","IRASF10565","IRAS19542","ESO148","NGC6240","IRASF17138","IRAS17578","ESO173","UGC11763","NGC3110","IC4734","NGC5135","ESO221","IC5179","UGC2982","NGC2369","ESO286","ESO320","NGC5734"]
res=[1.738,.793,2.671,1.323,1.630,1.320,1.789,.445,.691,.629,1.340,1.132,.728,1.490,1.076,1.028,.771,.564,np.nan,1.489,.666,2.344,2.5]
glo=[1.464,.555,1.706,1.140,1.534,.878,1.825,.387,.693,.587,1.283,.865,np.nan,1.044,1.063,1.100,.651,.842,np.nan,1.133,.587,2.254,2.4]
notP1={"IRASF05189","IRASF10565","IRAS19542","ESO148","NGC6240","UGC11763"}   # asterisks in Table 1: not in Paper I -> flag unknown
CON={"IRAS17208","IRAS17578","ESO320"}; VIB={"IRASF14378","ESO173"}
cls={n:("unclassified" if n in notP1 else "CON" if n in CON else "vib-detected non-CON" if n in VIB else "vib non-detected") for n in names}
t=pd.DataFrame(dict(src=names,cls=[cls[n] for n in names],res=res,glo=glo))
print("\nNishimura HCN/HCO+ 3-2 by Paper-I class (NGC5734 limits >2.5/>2.4 treated as 2.5/2.4 in medians; UGC2982 no value):")
for c,g in t.groupby("cls"): print("  %-22s n=%2d resolved median %.2f [%s]  global median %.2f"%(c,len(g),g.res.median(),", ".join("%.2f"%x for x in g.res.dropna()),g.glo.median()))
a=t[t.cls=="CON"]; b=t[t.cls=="vib non-detected"]; ab=t[t.cls.isin(["CON","vib-detected non-CON"])]
for lab,x in(("CON",a),("any vib det",ab)):
    for col in("res","glo"):
        u=mannwhitneyu(x[col].dropna(),b[col].dropna()); print("  Mann-Whitney %s vs non-detected, %s: U=%.1f p=%.3f (n=%d vs %d)"%(lab,col,u.statistic,u.pvalue,x[col].notna().sum(),b[col].notna().sum()))
print("  excluding NGC5734 limit: non-detected resolved median %.2f"%b[b.src!="NGC5734"].res.median())
print("  ours: N 0.99 (55 pc) / 0.69 (250 pc), S 0.44 / 0.37")
t.to_csv('/workspace/out/q51/nishimura_split.csv',index=False); df.to_csv('/workspace/out/q51/falstad_table3.csv',index=False)
# ---- continuum-derived gas surface density at N (our Q6 numbers) vs CON threshold N_H2>1e24
NH=1e24; Sig_thr=NH*2.8*1.6726e-24/2.0e-4*1.0   # Msun/pc2: 1 Msun/pc2 = 2.0e-4 g/cm2
print("\nCON column NH2=1e24 cm-2 -> Sigma_gas = %.1e Msun/pc2 ; our N (0.24\", dust upper bound, Td 50-100K) 7e2-1.5e3 ; S 2.3e3-5.0e3 (Q6)"%Sig_thr)
