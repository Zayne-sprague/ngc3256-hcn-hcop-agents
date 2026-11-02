# Q5: from line ratios to abundance ratios with RADEX (LVG, dv=150 km/s, Tbg 2.73 K, H2 collisions), assuming both lines of a ratio share the same Tkin, n, filling and velocity width.
# Beam-averaged brightness temperatures are only 0.01-0.1 K (compact emission in a 4" beam), so absolute columns need an unknown filling factor eta; the RATIO of T_R at fixed N(ref) does not.
#   required X = N(CS)/N(HCN) to give R_K(CS 7-6/HCN 4-3) for N(HCN) in 1e14,1e15,1e16 per 150 km/s (and the implied filling eta = T_b,obs/T_R(HCN 4-3) of the reference line).
#   limits   N(SO)/N(HCO+) < from SO 8_8-7_7/HCO+ 4-3 (Xb47) and SO 6_5-5_4/HCO+ 3-2 (Xb4f/Xb4d), 3 sigma.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
R=pd.read_csv('/workspace/out/q5/ratios.csv')
g3=np.load('/workspace/out/q3/grid_lvg.npz'); g5=np.load('/workspace/out/q5/radex_cs_so.npz'); TK,LN,LNC=g3['TK'],g3['LN'],g3['LNC']
def solveN(trarr,target):     # trarr over LNC (T_R), return logN where T_R=target (monotone part), nan otherwise
    ok=np.isfinite(trarr)&(trarr>0); x=LNC[ok]; y=np.log10(trarr[ok]); 
    if len(x)<3: return np.nan
    imax=np.argmax(y); x,y=x[:imax+1],y[:imax+1]
    t=np.log10(target)
    return float(np.interp(t,y,x)) if y[0]<=t<=y[-1] else np.nan
def TR(arr,iT,iN,lnN,k): 
    ok=np.isfinite(arr[iT,iN,:,k,0]); return float(10**np.interp(lnN,LNC[ok],np.log10(np.clip(arr[iT,iN,:,k,0][ok],1e-30,None))))
sel=lambda t,reg,ous='Xb47': R[(R.ratio==t)&(R.region==reg)&(R.ous==ous)&R.note.str.contains('shape=profile')].iloc[0]
out=[]
for reg in ('N','S','Bridge'):
    r1=sel('CS 7-6 / HCN 4-3',reg); r2=sel('SO 8_8-7_7 / HCO+ 4-3',reg)
    for tk in (28,56,80,160):
        iT=list(TK).index(tk)
        for ln in (4.0,4.5,5.0,5.25,5.5,6.0,7.0):
            iN=int(np.argmin(abs(LN-ln)))
            for lNr in (14,15,16):
                T_hcn=TR(g3['hcn'][...,2,:][...,None,:] if False else g3['hcn'],iT,iN,lNr,2) if False else float(10**np.interp(lNr,LNC,np.log10(np.clip(g3['hcn'][iT,iN,:,2,0],1e-30,None))))
                T_hco=float(10**np.interp(lNr,LNC,np.log10(np.clip(g3['hcop'][iT,iN,:,2,0],1e-30,None))))
                lNcs=solveN(g5['cs'][iT,iN,:,0,0],r1.R_K*T_hcn)
                lNso=solveN(g5['so'][iT,iN,:,1,0],r2.R_K_UL3*T_hco)
                out.append(dict(region=reg,Tk=tk,logn=ln,logN_ref=lNr,T_R_HCN43=T_hcn,logN_CS=lNcs,logX_CS_HCN=lNcs-lNr,logN_SO_UL3=lNso,logX_SO_HCOp_UL3=lNso-lNr,R_K_CS_HCN=r1.R_K,R_K_SO_HCOp_UL3=r2.R_K_UL3))
D=pd.DataFrame(out); D.to_csv('/workspace/out/q5/abundance_ratios_b7.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',700)
piv=D.pivot_table(index=['region','Tk','logn'],columns='logN_ref',values='logX_CS_HCN'); print('log10 X=N(CS)/N(HCN) required\n',piv.round(2).to_string())
piv=D.pivot_table(index=['region','Tk','logn'],columns='logN_ref',values='logX_SO_HCOp_UL3'); print('log10 N(SO)/N(HCO+) 3sig upper limit\n',piv.round(2).to_string())

# ---- Band 6: SO 6_5-5_4 / HCO+ 3-2 (SO-only fit = upper bound if CH3OH ladder contributes): limit = (W+3 sigma)/W_HCO+ in K units
out=[]
for ous in ('Xb4f','Xb4d'):
  for reg in ('N','S','Bridge','E1'):
    r=R[(R.ratio=='SO 6_5-5_4 / HCO+ 3-2')&(R.region==reg)&(R.ous==ous)&R.note.str.startswith('SO only')&R.note.str.contains('shape=profile')].iloc[0]
    kf=(267.55763/251.82577)**2; ul=(r.W_num+3*r.sig_used)/r.W_den/kf**-1 if False else (r.W_num+3*r.sig_used)/r.W_den*(267.55763/251.82577)**2
    for tk in (40,80,160):
      iT=list(TK).index(tk) if tk in TK else int(np.argmin(abs(TK-tk)))
      for ln in (4.5,5.0,5.25,6.0):
        iN=int(np.argmin(abs(LN-ln)))
        for lNr in (14,15):
            T_hco=float(10**np.interp(lNr,LNC,np.log10(np.clip(g3['hcop'][iT,iN,:,1,0],1e-30,None))))
            lNso=solveN(g5['so'][iT,iN,:,0,0],ul*T_hco)
            out.append(dict(ous=ous,region=reg,R_K_meas=r.R_K,R_K_UL=ul,Tk=tk,logn=ln,logN_ref=lNr,logN_SO_over_HCOp_UL=lNso-lNr))
D2=pd.DataFrame(out); D2.to_csv('/workspace/out/q5/abundance_limits_b6.csv',index=False)
print(D2.pivot_table(index=['ous','region','Tk','logn'],columns='logN_ref',values='logN_SO_over_HCOp_UL').round(2).to_string())
print(D2.groupby(['ous','region'])[['R_K_meas','R_K_UL']].first())
