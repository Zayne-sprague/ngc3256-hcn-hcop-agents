# Q57: best-fit (max marginal-likelihood cell, g=0) observables and per-line residuals for each hypothesis of s353, plus the contrast each reproduces; uses the s353 cubes.
# usage: python3 -I s354_q57_bestfits.py [XLO XHI] -> out/q57/s354_<XLO>-<XHI>.txt
import sys,numpy as np; sys.path.insert(0,'/workspace/scripts')
a=sys.argv[1:]; sys.argv=['x']+a
import s353_q57_joint_bracket as S
xlo,xhi=(float(a[0]),float(a[1])) if a else (1.6,2.0)
out=[]
def P(t): print(t,flush=True); out.append(t)
CN=S.build('N')[1].astype(float); CS=S.build('S')[1].astype(float)   # g=0: [X,T,n,NO,A]
xs=(S.LXg>=xlo-1e-6)&(S.LXg<=xhi+1e-6); CN=CN[xs]; CS=CS[xs]; LX=S.LXg[xs]
M={}
def best(pos,idx,X):
    I,e=S.data(pos); m=S.model_cube(X); mm=m[:,idx[0],idx[1],idx[2],idx[3]].astype(float)
    # best ln eta
    bl=None
    for le in S.LETA:
        mod=mm*np.exp(le); s2=e**2+(0.05*I)**2+(0.1*mod)**2; c=((I-mod)**2/s2).sum()
        if bl is None or c<bl[0]: bl=(c,le,mod)
    return I,e,bl
def report(name,iN,iS):
    # iN,iS = (ix,iT,in,iO,iA) cells
    rows=[]
    for pos,ii in (('N',iN),('S',iS)):
        I,e,bl=best(pos,ii[1:],LX[ii[0]])
        rows.append((pos,ii,I,e,bl))
    P(f'-- {name}')
    for pos,ii,I,e,bl in rows:
        P(f'  {pos}: Tk={S.TK[ii[1]]:.0f} log n={S.LNs[ii[2]]:.1f} logN(HCO+)={S.NOg[ii[3]]:.1f} A={10**S.LAg[ii[4]]:.2f} X={10**LX[ii[0]]:.0f} eta={np.exp(bl[1]):.2f} chi2={bl[0]:.1f}')
        P('     obs  HCN1,HCN3,HCO1,HCO3,H13CN1,H13CO1 = '+' '.join(f'{x:6.2f}' for x in I)+'\n     mod  '+' '.join(f'{x:6.2f}' for x in bl[2]))
    mN,mS=rows[0][4][2],rows[1][4][2]
    P(f'  model contrast N/S of HCN/HCO+: 1-0 {(mN[0]/mN[2])/(mS[0]/mS[2]):.2f}, 3-2 {(mN[1]/mN[3])/(mS[1]/mS[3]):.2f}  (obs 1.45, 1.50)')
def am(C,axes): return None
# full per-nucleus best cells
bN=np.unravel_index(CN.argmax(),CN.shape); bS=np.unravel_index(CS.argmax(),CS.shape)
report('Hfull best',bN,bS)
# Hexc: shared (NO,A,X) ; T,n free
LN_=CN.max((1,2)); LS_=CS.max((1,2)); J=LN_*LS_; sh=np.unravel_index(J.argmax(),J.shape)
iN=(sh[0],)+tuple(np.unravel_index(CN[sh[0],:,:,sh[1],sh[2]].argmax(),CN.shape[1:3]))+(sh[1],sh[2])
iS=(sh[0],)+tuple(np.unravel_index(CS[sh[0],:,:,sh[1],sh[2]].argmax(),CS.shape[1:3]))+(sh[1],sh[2])
report('Hexc best (same A, N(HCO+), X; Tk, n free)',iN,iS)
# Hcol: shared (A,X)
LN_=CN.max((1,2,3)); LS_=CS.max((1,2,3)); J=LN_*LS_; sh=np.unravel_index(J.argmax(),J.shape)
def cell(C,x,a): 
    sub=C[x,:,:,:,a]; t=np.unravel_index(sub.argmax(),sub.shape); return (x,t[0],t[1],t[2],a)
report('Hcol best (same A, X; N(HCO+), Tk, n free)',cell(CN,*sh),cell(CS,*sh))
# HcolX: shared A only
LN_=CN.max((0,1,2,3)); LS_=CS.max((0,1,2,3)); a_=(LN_*LS_).argmax()
def cellA(C,a): 
    sub=C[...,a]; t=np.unravel_index(sub.argmax(),sub.shape); return (t[0],t[1],t[2],t[3],a)
report('HcolX best (same A; N(HCO+), X, Tk, n free)',cellA(CN,a_),cellA(CS,a_))
# Habun: shared X,T,n,NO ; A free
J=CN.max(4)*CS.max(4); sh=np.unravel_index(J.argmax(),J.shape)
report('Habun best (same Tk,n,N(HCO+),X; A free)',sh+(CN[sh].argmax(),),sh+(CS[sh].argmax(),))
# Posterior predictive: contrast reproduced in the excitation-only family, fraction of Hexc posterior mass with 1-0 contrast >=1.35 ?
open(S.O+f's354_{xlo}-{xhi}.txt','w').write('\n'.join(out))
