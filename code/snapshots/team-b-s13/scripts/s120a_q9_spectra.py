# Q9: Xb55 spw22 (163.24-165.11 GHz obs) spectra at N, S, Bridge (peak pixel + aperture) in the CH3OH J_-1-J_0 E (Q-branch) region; plot.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
z=np.load('/workspace/out/regspec/Xb55.spw22.cube.pbcor.npz'); print([k for k in z.files])
nu=z['nu']; regs=sorted({k.rsplit('_',1)[0] for k in z.files if k!='nu' and k.endswith(('_pix','_ap'))}); print(regs)
fig,ax=plt.subplots(len(regs),1,figsize=(10,2.2*len(regs)),sharex=True)
f0=165.061130
for a,r in zip(np.atleast_1d(ax),regs):
    for st,c in (('pix','k'),('ap','r')):
        if f'{r}_{st}' in z.files:
            y=z[f'{r}_{st}']; v=C*(1-nu/(f0*1e9)); a.step(v,y*1e3,c,lw=.6,label=st)
    a.set_ylabel(r+' mJy/b'); a.axvline(2780,ls=':'); a.set_ylim(-12,25)
plt.xlabel('v (radio, rel. 165.061 GHz) km/s'); plt.savefig('/workspace/out/q9/s120a_spectra.png',dpi=90)
for r in regs:
  y=z[f'{r}_pix']; v=C*(1-nu/(f0*1e9)); m=np.isfinite(y)
  print(r,'n',m.sum(),'rms(|v-2780|>1200)',np.nanstd(y[abs(v-2780)>1200])*1e3,'median far',np.nanmedian(y[abs(v-2780)>1200])*1e3)
