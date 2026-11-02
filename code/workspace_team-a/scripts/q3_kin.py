# Q3 step 5: kinematic context around S: E-W velocity cut (disk rotation axis) vs N-S cut; Xb4d vs Xb4b flux of the Soutflow blue lobe (clumpiness / recovered flux).
import sys; sys.path.insert(0,'/workspace/scripts')
from q3_pos import *
import pandas as pd
c={sp:cube(sp,fw=1.05) for sp in ('HCN','HCO+')}
print('E-W cut through S (y=-2.55), HCO+ 3-2 1.05": centroid v (2300-3150, I>0) and flux')
rows=[]
for x in np.arange(-4.35,4.2,0.6):
    cq=c['HCO+']; s=spec(cq,x,-2.55); v=cq['v']; m=(v>2300)&(v<3150); I=s[m].sum()*cq['dv']
    rows.append(dict(x=x,I=I,vmean=(s[m]*v[m]).sum()/s[m].sum() if I>0.2 else np.nan,vpeak=v[np.argmax(s)])); 
K=pd.DataFrame(rows);print(K.round(1).to_string());K.round(2).to_csv('/workspace/out/q3/q3_EW_cut.csv',index=False)
# lobe flux: sum in box (x -1.7..1.3, y -6.0..-3.3) of blue window 2300-2575, Xb4d native vs Xb4b native
for ous in ('Xb4d','Xb4b'):
    for sp in ('HCN','HCO+'):
        cc=cube(sp,ous,crop=8.,fw=None); m0,n=wsum(cc,2300,2575); sg=wnoise(cc,2300,2575)[0]
        box=(cc['dx']>-1.7)&(cc['dx']<1.3)&(cc['dy']>-6.0)&(cc['dy']<-3.3); ppb=OM(cc['beam'])/cc['pix']**2
        F=m0[box].sum()/ppb; eF=np.sqrt(box.sum()/ppb)*np.median(sg[box])   # approx. (noise correlated over beam): sqrt(Nbeam)*sigma_beam
        print(f'{ous} {sp} blue 2300-2575 flux in box (no S/N clip) = {F:.3f} +- {eF:.3f} Jy km/s; peak {m0[box].max():.3f} Jy/beam km/s; beam {cc["beam"][0]:.2f}x{cc["beam"][1]:.2f}')
