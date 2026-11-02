# Q3 round 2, fix 1: effect of a fixed Xb4f (3-2) flux loss of 0/30/50% (Xb55 and Xb47 unaffected, or Xb47 losing the same) on the ladder fit and on T(2-1)/T(3-2), T(4-3)/T(3-2).
# usage: python3 -I s55e_xb4f_loss.py -> out/q3b/s55e_xb4f_loss.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s55_lib as S, s55b_recovery_fit as B
L=S.L; F=B.F; log=[]; res={}
def P(s=''): print(s,flush=True); log.append(s)
for pos in ('N','S','E1','E2','W'):
    I,e=L.DATA[pos]; extra=B.extra_for(pos); rel=e.copy()
    if pos=='W': rel[2]=5.0
    En=S.Engine(F,S.covmat(rel,extra=extra))
    for lab,(k47,k55) in (('Xb4f only',(0.,0.)),('Xb4f and Xb47 equal',(1.,0.))):
        for Lf in (0,0.3,0.5):
            f=S.recovery_vec(Lf,0,k47,k55); y=np.log(I)-np.log(f); d=S.Post(F,En.lnlike(y)[0]).summary(); res[f'{pos}|{lab}|{Lf}']=d
            Ic=I/f; P(f'[{pos}] loss {int(Lf*100)}% in {lab}: T(2-1)/T(3-2) HCN {Ic[0]/Ic[1]:.2f} HCO+ {Ic[3]/Ic[4]:.2f}; T(4-3)/T(3-2) HCN {Ic[2]/Ic[1]:.2f} HCO+ {Ic[5]/Ic[4]:.2f} | log n {d["log n"][2]:.2f} [{d["log n"][1]:.2f},{d["log n"][3]:.2f}] Tk {d["Tk"][2]:.0f} [{d["Tk"][1]:.0f},{d["Tk"][3]:.0f}] log A {d["log A"][2]:.2f} | P(n>1e5)={d["P"]["n>1e5"]:.2f} P(Tk<30)={d["P"]["Tk<30"]:.2f}')
json.dump(res,open(S.OUT+'s55e_xb4f_loss.json','w'),indent=1,default=float); open(S.OUT+'s55e_xb4f_loss.txt','w').write('\n'.join(log))
