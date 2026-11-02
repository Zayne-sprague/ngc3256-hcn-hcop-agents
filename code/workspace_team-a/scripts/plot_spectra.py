import sys,warnings; warnings.filterwarnings('ignore')
import numpy as np, glob, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ks=sorted(glob.glob('/workspace/out/spec_*.npz'))
for grp in ['Xb55','Xb47','Xb4f','Xb4d']:
    fs=[k for k in ks if grp+'.' in k]
    fig,ax=plt.subplots(len(fs),1,figsize=(14,3*len(fs)))
    for a,k in zip(np.atleast_1d(ax),fs):
        d=np.load(k)
        for n,c in [('N','C0'),('S','C3')]: a.plot(d['freq']/1e9,d[n]*1e3,c,lw=.7,label=n)
        a.set_title(k.split('spec_')[1]); a.legend(); a.set_ylabel('mJy/beam')
    plt.tight_layout(); plt.savefig(f'/workspace/out/spec_{grp}.png',dpi=70); plt.close()
