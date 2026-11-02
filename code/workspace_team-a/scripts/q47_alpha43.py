# Q47(a'): predicted alpha(4-3)=I(H13CN 4-3)/I(HCN 4-3) [Jy km/s; Jy/K ~ nu^2 so alpha=TR ratio*(345.34/354.51)^2] from our RADEX grids (out/q4c/grid4_*.npz, Tbg 2.73 K) for cold-thick, warm and hot cells at X=45, to test p0050's claim alpha=0.016-0.035 in every branch.
import numpy as np
H=np.load('/workspace/out/q4c/grid4_hcn.npz'); I=np.load('/workspace/out/q4c/grid4_h13cn.npz'); TK=H['TK']; LN=H['LN']; LNC=H['LNC']; F=(345.3398/354.5055)**2
print('T   logn  logN12  tau43(HCN)  alpha(X=32,45,63)')
for T in (15,20,30,50,100):
  for n in (4.5,5.0,5.5,6.0):
    for N in (16.0,16.75,17.5):
      it=list(TK).index(T); jn=int(round((n-LN[0])/.25)); a=[]
      iN=int(round((N-LNC[0])/.25)); t43=H['tau'][it,jn,iN,3]; trh=H['TR'][it,jn,iN,3]
      for X in (32,45,63):
        k=int(round((N-np.log10(X)-LNC[0])/.25)); a.append(I['TR'][it,jn,k,0]/trh*F if k>=0 else np.nan)
      if N==16.75: print(f'{T:3d} {n:4.1f} {N:6.2f} {t43:8.1f}   '+' '.join(f'{x:.3f}' for x in a))
