import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
d=np.load('/workspace/out/region_spectra.npz'); C=299792.458
fig,axs=plt.subplots(4,2,figsize=(16,14))
for r,(ous,spw) in enumerate([('Xb4f','spw20'),('Xb4f','spw22'),('Xb4d','spw29'),('Xb4d','spw31')]):
    fr=d[f'{ous}__{spw}__freq']
    for c,reg in enumerate(['N','S']):
        s=d[f'{ous}__{spw}__{reg}__pk']; a=axs[r,c]
        a.plot(fr,s*1e3,lw=.4,c='0.6'); a.plot(fr,uniform_filter1d(s,12)*1e3,'k',lw=1)
        a.set_title(f'{ous} {spw} {reg}'); a.axhline(0,c='r',lw=.3)
        for nm,f0 in [('HC3N',254.6995),('34SO',253.208),('NSe',253.570),('NSf',253.968),('SiS',254.103),('H2COH+',252.870),('SO',251.8258),('CH3OH',251.906),('lC3H',251.434),('CH3CHO',251.364)]:
            if f0==0: continue
            fs=f0*(1-2780/C)
            if fr.min()<fs<fr.max(): a.axvline(fs,c='b',ls=':'); a.text(fs,a.get_ylim()[1]*0.8,nm,fontsize=7,rotation=90)
plt.tight_layout(); plt.savefig('/workspace/out/weak_b6.png',dpi=55)
