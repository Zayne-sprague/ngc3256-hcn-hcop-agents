# Q11: translate the HCN/HCO+ annulus-ratio change (inner r<=4" -> outer 8-12") into RADEX (n, X=N(HCN)/N(HCO+)) changes, using the existing grids out/q1/radex_{hcn,hco+}.npz
# (LAMDA, H2, Tbg 2.73 K, dv=150 km/s, T_R of 2-1,3-2,4-3). Tb ratio = flux ratio *1.0126 (nu^2). Inner: R21=0.587,R32=0.447 (Xb4f); outer: R21=0.498 (2-1 8-12") and R32=0.447*q with q = true 3-2 fall 0.58-0.77 (forward model).
import numpy as np, itertools
TK=np.array([10,15,20,30,50,80,150.]); LN=np.arange(3,8.01,0.25); LNC=np.arange(13,18.51,0.25)
H=np.load('/workspace/out/q1/radex_hcn.npz')['TR']; O=np.load('/workspace/out/q1/radex_hco+.npz')['TR']   # [T,n,N,J]
def Rj(it,iN,jH,jO,J): return H[it,iN,jH,J]/O[it,iN,jO,J]
f=1.0126; inner=(0.579*f,0.442*f)
def match(target,tol,it,iN):
    t21,t32=target; m=np.zeros((len(LNC),len(LNC)),bool)
    R21=H[it,iN,:,0][:,None]/O[it,iN,:,0][None,:]; R32=H[it,iN,:,1][:,None]/O[it,iN,:,1][None,:]
    return (abs(R21/t21-1)<tol)&(abs(R32/t32-1)<tol)
print('Inner target R21,R32 (Tb)',np.round(inner,3))
for q in (0.58,0.68,0.77,0.85):
    outer=(0.492*f,0.442*f*q)
    print(f'\n== outer: R21={outer[0]:.3f} (2-1 8-12") R32={outer[1]:.3f} (3-2 fall q={q} of inner)')
    for T in (20,30,50,80):
        it=list(TK).index(T); rows=[]
        for iN,n in enumerate(LN):
            mi=match(inner,0.07,it,iN); mo=match(outer,0.07,it,iN)
            if mi.any() and mo.any():
                # same N(HCO+) allowed to differ; summarise X=NH/NO range = 10^(lNH-lNO)
                xi=[LNC[a]-LNC[b] for a,b in zip(*np.where(mi))]; xo=[LNC[a]-LNC[b] for a,b in zip(*np.where(mo))]
                rows.append((n,min(xi),max(xi),min(xo),max(xo)))
        print(f' Tk={T:3d} K: log n (cm-3) where both zones match: '+(', '.join(f'{r[0]:.2f}[inner logX {r[1]:+.2f}..{r[2]:+.2f}; outer {r[3]:+.2f}..{r[4]:+.2f}]' for r in rows[:6]) if rows else 'none'))
# direct sensitivity: at Tk=30 K, N(HCO+)=1e14.5, X=logN diff 0.0: R21,R32 vs n
it=3
print('\nSensitivity at Tk=30 K, log N(HCO+)=14.5, X=1 and 0.5 : n -> R21,R32 (Tb)')
for iN,n in enumerate(LN):
    if 3.5<=n<=7: 
        a=LNC.tolist().index(14.5)
        print(f' logn={n:.2f}',[ (round(H[it,iN,a+d,0]/O[it,iN,a,0],3),round(H[it,iN,a+d,1]/O[it,iN,a,1],3)) for d in (0,-2)])
