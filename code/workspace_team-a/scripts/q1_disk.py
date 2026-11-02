# Q1: 'disk' candidate = Xb55 HCO+ 2-1 S/N>5 pixels at r>8" (10" smoothed map). Mean brightness (Jy/beam km/s, 10" beam) and ratios over the mask with different zero-level treatments.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
FW=10.0
keys=['HCN 2-1','HCO+ 2-1','HCN 3-2@Xb4f','HCO+ 3-2@Xb4f','HCN 4-3','HCO+ 4-3']
D={k:load(k) for k in keys}; S={k:smooth(D[k],FW) for k in keys}; P={k:noise_profile(D[k],FW,rbins=(0,4,8,12,16,20,24,28)) for k in keys}
ref='HCO+ 2-1'; d0=D[ref]; r0=np.hypot(d0['dx'],d0['dy'])
sig=np.vectorize(lambda q:sigma_at(P[ref],q))(r0)
mask=(S[ref]/sig>5)&(r0>8)&(r0<24)
print('mask pixels',mask.sum(),'(independent 10" beams ~',round(mask.sum()*1.2**2/(1.13*FW**2),1),')')
px,py=d0['dx'][mask],d0['dy'][mask]; pr=np.hypot(px,py)
print('noise profile Xb55 HCO+ 2-1 at 10" (r_in,r_out,sigma,null median):',[(a,b,round(s,3),round(m,3)) for a,b,s,m,_ in P[ref]])
print('Band-5 PB gain (HWHM 17") at r=8,12,16,20 :',[round(float(np.exp(-np.log(2)*(r/17)**2)),2) for r in (8,12,16,20)])
rows=[]
vals={}
for k in keys:
    d=D[k]; sm=S[k]; v=np.array([sample(sm,d,x,y) for x,y in zip(px,py)])
    # zero levels: (a) none; (b) median of the same 10" map in the 14-22" ring (used in p0001); (c) median of the null-window maps (baseline offset of the line-free channels) in the mask
    rr=np.hypot(d['dx'],d['dy']); ring=np.nanmedian(sm[(rr>14)&(rr<22)&(sm!=0)])
    nullmed=np.nanmedian(np.concatenate([smooth(d,FW,n)[(rr>8)&(rr<24)] for _,n in d['nulls']])) if d['nulls'] else np.nan
    e=np.sqrt(np.mean([sigma_at(P[k],q)**2 for q in pr]))/np.sqrt(max(1,mask.sum()*1.2**2/(1.13*FW**2)))
    vals[k]=(v.mean(),ring,nullmed,e)
    rows.append(dict(line=k,mean_mask=v.mean(),ring_level=ring,null_level=nullmed,e_mean=e,mean_ring_corr=v.mean()-ring,mean_null_corr=v.mean()-nullmed))
Tm=pd.DataFrame(rows); print(Tm.round(3).to_string()); Tm.to_csv('/workspace/out/q1/q1_disk_levels.csv',index=False)
def Rt(nk,dk,mode):
    f=lambda k:{'none':vals[k][0],'ring':vals[k][0]-vals[k][1],'null':vals[k][0]-vals[k][2]}[mode]
    kT=lambda k:Tb_factor(LINES[k][2]/(1+VSYS/C),FW)
    a,b=f(nk)*kT(nk),f(dk)*kT(dk); ea,eb=vals[nk][3]*kT(nk),vals[dk][3]*kT(dk)
    R=a/b; return R,abs(R)*np.hypot(ea/a,eb/b),a/ea,b/eb
out=[]
for mode in ('none','ring','null'):
    for nm,nk,dk in [('HCN/HCO+ 2-1','HCN 2-1','HCO+ 2-1'),('HCN/HCO+ 3-2','HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'),('HCN/HCO+ 4-3','HCN 4-3','HCO+ 4-3'),
                     ('HCO+ 3-2/2-1','HCO+ 3-2@Xb4f','HCO+ 2-1'),('HCN 3-2/2-1','HCN 3-2@Xb4f','HCN 2-1')]:
        R,eR,sa,sb=Rt(nk,dk,mode); out.append(dict(zero_level=mode,ratio=nm,R=R,eR_stat=eR,snr_num=sa,snr_den=sb))
O=pd.DataFrame(out); O.to_csv('/workspace/out/q1/q1_disk_ratios.csv',index=False); print(O.round(3).to_string())
