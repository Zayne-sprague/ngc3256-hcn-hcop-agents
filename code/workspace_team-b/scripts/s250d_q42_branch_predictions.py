# Q42 diagnostic: for the best-fit cell of each Tk branch (cold<=30, mid, hot>=100) of the N one-zone fit (Harada 12C + H13CO+, variant 'H', band priors 10/15%),
# list the model N(HCN), tau(HCN 1-0, 4-3), the predicted alpha = I(H13CN 4-3)/I(HCN 4-3) [Jy km/s] and H13CN/HCN 1-0 and HCO+ ratios, versus the data.
# usage: python3 -I s250d_q42_branch_predictions.py -> out/q42/branch_predictions.txt
import sys,os,numpy as np
sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']
import s250b_q42_ratio_fit as q
s=q.s; m=q.m; d=q.d; S=q.S235; TK=m.TK; LN=m.LN; NF=m.NF
L=[]
def P(x): print(x); L.append(x)
k=0; sg=S.SIGS[k]
for lx in (1.5,1.65,1.8):
    q.CUR['terms']=[]; X=10**lx
    f=S.O+f'fit_N_H_s{k}_lX{int(round(lx*100))}.npz'; r=np.load(f)
    qq=s.iso('H')['N']; s.X_iso={'H':qq['H'],'O':qq['O'],'C':[]}
    s.MOD13['H']=s.lnI13(d['hcn'],lx); s.MOD13['O']=s.lnI13(d['hcop'],lx)
    m.setup('N',*sg,['H','O','C'])
    Rh4=q.ratio_model('H',2,lx); Rh1=q.ratio_model('H',0,lx); Ro1=q.ratio_model('O',0,lx); Ro3=q.ratio_model('O',1,lx)
    P(f'--- X={X:.0f} (band priors {sg})')
    for lab,mk in (('cold Tk<=30',TK<=30),('mid',(TK>30)&(TK<100)),('hot >=100',TK>=100)):
        cm=r['chimin'].copy(); cm[~mk,:]=np.inf; it,jn=np.unravel_index(cm.argmin(),cm.shape)
        tot=m.X.pi.copy(); F={sp:m.chi_s(it,jn,sp) for sp in 'HOC'}
        for sp in 'HOC': tot=tot+F[sp].min(1)
        ib=int(np.argmin(tot)); iN=int(np.argmin(F['H'][ib])); iO=int(np.argmin(F['O'][ib]))
        tau4=d['hcn'][it,jn,:,2,2]; tau1=d['hcn'][it,jn,:,0,2]
        from scipy.interpolate import CubicSpline
        t4=CubicSpline(m.LNC,np.nan_to_num(tau4),axis=0)(NF[iN]); t1=CubicSpline(m.LNC,np.nan_to_num(tau1),axis=0)(NF[iN])
        P(f' {lab:12s}: Tk={TK[it]:.0f} K log n={LN[jn]:.1f} chi2={tot[ib]:.2f} logN(HCN)={NF[iN]:.2f} tau(HCN 1-0)={t1:.1f} tau(HCN 4-3)={t4:.2f} | alpha_mod(H13CN4-3/HCN4-3 flux)={Rh4[it,jn,iN]:.3f} (obs <0.12 3sig; 0.045+-0.04) H13CN/HCN 1-0={Rh1[it,jn,iN]:.3f} (obs 0.0207+-0.0034 @2.2\") H13CO+/HCO+ 1-0={Ro1[it,jn,iO]:.3f} (obs 0.0231+-0.0031 @2.2\", Harada 0.046) 3-2={Ro3[it,jn,iO]:.3f} (Harada 0.0284)')
open(q.O+'branch_predictions.txt','w').write('\n'.join(L))
# ---- part 2: per-datum ratio terms of the 'T' (common-taper 2.2" ratios) variant at N for the best cell of each branch (cell chosen from the T fit itself)
L2=[]
def P2(x): print(x); L2.append(x)
for var in ('T','T32'):
  for k in (0,1):
    for lx in (1.65,1.8,1.95):
        f=q.O+f'fit_N_{var}_s{k}_lX{int(round(lx*100))}.npz'
        if not os.path.exists(f): continue
        r=np.load(f); q.CUR['terms']=q.terms('N',var); q.RAT.clear()
        for (spc,li,Ro,sg_) in q.CUR['terms']: q.RAT[(spc,li)]=q.ratio_model(spc,li,lx)
        m.setup('N',*S.SIGS[k],['H','O','C']); s.X_iso={'H':[],'O':[],'C':[]}
        P2(f'--- N variant {var}, band priors {S.SIGS[k]}, X={10**lx:.0f}')
        for lab,mk in (('cold Tk<=30',TK<=30),('mid',(TK>30)&(TK<100)),('hot >=100',TK>=100)):
            cm=r['chimin'].copy(); cm[~mk,:]=np.inf; it,jn=np.unravel_index(cm.argmin(),cm.shape)
            tot=m.X.pi.copy(); F={sp:m.chi_s(it,jn,sp) for sp in 'HOC'}
            for sp in 'HOC': tot=tot+F[sp].min(1)
            ib=int(np.argmin(tot)); iN=int(np.argmin(F['H'][ib])); iO=int(np.argmin(F['O'][ib]))
            parts=[]
            for (spc,li,Ro,sg_) in q.CUR['terms']:
                i_=iN if spc=='H' else iO; Rm=q.RAT[(spc,li)][it,jn,i_]; s2=sg_**2+(0.1*Rm)**2
                parts.append(f"{'H13CN' if spc=='H' else 'H13CO+'}({['1-0','3-2','4-3'][li]}) mod {Rm:.4f} obs {Ro:.4f} chi2 {(Rm-Ro)**2/s2:.2f}")
            P2(f' {lab:12s} Tk={TK[it]:.0f} logn={LN[jn]:.1f} total chi2 {tot[ib]:.2f}: '+' | '.join(parts))
open(q.O+'branch_ratio_terms.txt','w').write('\n'.join(L2))
