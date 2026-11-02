# Xb55 spw20/spw22 frequency overlap (obs 164.536-165.110 GHz): independent-baseband consistency of the continuum-subtracted-free spectra (N, S beam apertures, and r=8" sums)
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
z20=np.load('/workspace/out/regspec/Xb55.spw20.cube.pbcor.npz'); z22=np.load('/workspace/out/regspec/Xb55.spw22.cube.pbcor.npz')
lo,hi=max(z20['nu'].min(),z22['nu'].min()),min(z20['nu'].max(),z22['nu'].max()); print(f'overlap obs {lo/1e9:.3f}-{hi/1e9:.3f} GHz = rest {lo/1e9/(1-VSYS/C):.3f}-{hi/1e9/(1-VSYS/C):.3f} GHz at {VSYS}')
for reg in ('N','S'):
    nu=z20['nu']; m=(nu>=lo)&(nu<=hi); a=z20[reg+'_ap'][m]
    b=np.interp(nu[m],z22['nu'][np.argsort(z22['nu'])],z22[reg+'_ap'][np.argsort(z22['nu'])])
    d=a-b; print(f'  {reg}: n={m.sum()} ch; mean(spw20)={np.nanmean(a)*1e3:.2f} mean(spw22)={np.nanmean(b)*1e3:.2f} mJy/beam; diff mean {np.nanmean(d)*1e3:.2f} std {np.nanstd(d)*1e3:.2f} mJy (per-channel rms ~{np.nanstd(a)*1e3:.2f}); sum over overlap: {np.nansum(a)*abs(np.median(np.diff(nu)))/1e9*C/177.0:.3f} vs {np.nansum(b)*abs(np.median(np.diff(nu)))/1e9*C/177.0:.3f} Jy/beam km/s')
print('CH3OH 165.050/.061/.099 rest = obs',[round(f*(1-VSYS/C),3) for f in (165.050,165.061,165.099)],'GHz: inside spw22 only (spw20 starts at 164.536 obs? no: CH3OH obs is 163.5 < 164.536) -> no independent overlap check possible')
