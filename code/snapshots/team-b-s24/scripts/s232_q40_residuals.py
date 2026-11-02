# Q40: per-datum residuals at the best-fit cell of the N one-gas fit with isotopologue data (variant A/H, X=40,60,90)
import sys,os,numpy as np
sys.path.insert(0,'/workspace/scripts'); sys.argv=['x']; import s230_q40_isotopologue_fit as s
m=s.m
for v in ('A','H'):
  for X in (40,60,90):
    pos='N'; q=s.iso(v)[pos]
    s.X_iso={'H':q['H'],'O':q['O'],'C':[]}; s.MOD13['H']=s.lnI13(s.d['hcn'],np.log10(X)); s.MOD13['O']=s.lnI13(s.d['hcop'],np.log10(X))
    m.setup(pos,0.10,0.15,['H','O','C'])
    o=np.load(f'/workspace/out/q40/fit_N_{v}_X{X}.npz'); cm=o['chimin']; best=np.unravel_index(cm.argmin(),cm.shape)
    # branch-wise best cells
    for lab,sel in (('best',None),('cold(Tk<=30)',m.TK<=30)):
        c=cm.copy()
        if sel is not None: c[~sel]=np.inf
        it,jn=np.unravel_index(c.argmin(),c.shape); tot=m.X.pi.copy(); info={}
        arg={}
        for sp in 'HOC':
            f=s.chi_s(it,jn,sp); tot=tot+f.min(1)
        ib=int(np.argmin(tot)); out=[]
        for sp in 'HOC':
            f=s.chi_s(it,jn,sp); iN=int(np.argmin(f[ib])); Mo=m.MOD[sp][it,jn][iN]; M3=s.MOD13[sp][it,jn][iN] if sp in s.MOD13 else Mo
            i0,i1=m.SPL[sp]; off=(m.X.o3[ib],m.X.o6[ib])
            mod12=np.exp(Mo+np.array(off)); mod13=np.exp(M3+np.array(off))
            out.append(f'{sp}: logN={m.NF[iN]:.2f} 12C model {mod12[0]:.1f}/{mod12[1]:.1f} (obs {np.exp(m.X.y[i0]):.1f}/{np.exp(m.X.y[i1]):.1f}); 13C model {mod13[0]:.2f}/{mod13[1]:.2f}; obs13 {[(round(I,2)) for _,I,_ in s.X_iso[sp]]}')
        print(f'{v} X={X} {lab}: Tk={m.TK[it]:.0f} logn={m.LN[jn]:.1f} chi2={tot[ib]:.1f} eta-offs {m.X.o3[ib]:.2f},{m.X.o6[ib]:.2f}\n   '+'\n   '.join(out))
