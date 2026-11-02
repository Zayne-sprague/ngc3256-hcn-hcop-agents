# Q3 round 2: probability that the disk (E1+E2+W joint) is colder / lower-n than the centre (N+S+Bridge joint), from the marginal posteriors of s55b (independent fits).
# usage: python3 -I s55f_tk_contrast.py -> out/q3b/s55f_tk_contrast.txt
import numpy as np
d=np.load('/workspace/out/q3b/s55b_posteriors.npz'); LT=d['LT']; LN=d['LN']; out=[]
keys=[k for k in d.files if 'joint' in k]
for sc in ('nominal','one-sided','stress','ratios-only'):
    kc=[k for k in keys if k.startswith(sc) and 'centre' in k][0]; kd=[k for k in keys if k.startswith(sc) and 'disk' in k][0]
    pc=d[kc]; pd_=d[kd]; Tc=pc.sum(1); Td=pd_.sum(1); nc=pc.sum(0); nd=pd_.sum(0)
    P_T=sum(Td[i]*Tc[:i].sum() for i in range(len(LT)))+0.5*np.sum(Td*Tc)     # P(Tdisk>Tcentre)
    P_n=sum(nd[i]*nc[:i].sum() for i in range(len(LN)))+0.5*np.sum(nd*nc)
    out.append(f'{sc:12s}: P(Tk_disk < Tk_centre)={P_T:.3f}; P(n_disk > n_centre)={P_n:.3f}')
print('\n'.join(out)); open('/workspace/out/q3b/s55f_tk_contrast.txt','w').write('\n'.join(out))
