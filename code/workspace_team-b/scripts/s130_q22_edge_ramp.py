# Q22: edge ramp in Xb55 spw16/18/20/22 pbcor spectra (peak pixel, continuum-subtracted by pipeline): profile vs channel from each edge, vs continuum.
import numpy as np, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
regs=['N','S','Bridge','E1','E2','W','Sblue']; spws=[16,18,20,22]; B=10
out={}; fig,ax=plt.subplots(2,4,figsize=(16,7),sharey='row')
for j,s in enumerate(spws):
    z=np.load(f'/workspace/out/regspec/Xb55.spw{s}.cube.pbcor.npz'); n=len(z['nu'])
    for r in regs:
        y=z[r+'_pix']*1e3
        lo=y[:B*30].reshape(30,B).mean(1); hi=y[::-1][:B*30].reshape(30,B).mean(1)
        mid=np.nanmedian(y[300:650]); rms=np.nanstd(y[300:650]-np.nanmedian(y[300:650]))
        out[f'{s}_{r}']=dict(lo=lo.tolist(),hi=hi.tolist(),mid=float(mid),rms=float(rms))
        ax[0,j].plot(np.arange(30)*B+B/2,lo,label=r); ax[1,j].plot(np.arange(30)*B+B/2,hi,label=r)
    ax[0,j].set_title(f'spw{s} low-freq edge (mJy/b)'); ax[1,j].set_title(f'spw{s} high-freq edge'); ax[0,j].axhline(0,c='k',lw=.5); ax[1,j].axhline(0,c='k',lw=.5)
ax[0,0].legend(fontsize=6); [a.set_xlabel('channels from edge (1.95 MHz)') for a in ax[1]]
plt.tight_layout(); plt.savefig('/workspace/out/q22/s130_edges.png',dpi=80)
json.dump(out,open('/workspace/out/q22/s130.json','w'))
for s in spws:
    print('spw',s)
    for r in regs:
        o=out[f'{s}_{r}']; print(f'  {r:7s} lo(ch0-10,30-40,100-110,200-210)',[round(o["lo"][k],2) for k in (0,3,10,20)],' hi',[round(o["hi"][k],2) for k in (0,3,10,20)],'rms',round(o['rms'],2))
