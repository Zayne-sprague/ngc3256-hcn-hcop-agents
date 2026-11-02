# Q20 step a: native-resolution spectra (peak pixel, Jy/beam) at S, N, ridge, lobe for every window that could hold the ridge/lobe wings:
# Xb47 HCN 4-3 (spw16), HCO+ 4-3 (spw18), H13CN 4-3 (+CO 3-2 edge, spw20), CS 7-6 (spw22); Xb4d SO 6_5-5_4 (spw31), HCN/HCO+ 3-2 for reference.
# Output: per-channel noise (line-free, native Jy/beam), spectra csv and figure.   usage: python3 -I s110a_q20_spectra.py
import sys; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q20'
WIN={'HCN43':('Xb47','spw16',354.50547e9),'HCO+43':('Xb47','spw18',356.73424e9),'H13CN43':('Xb47','spw20',345.33977e9),'CS76':('Xb47','spw22',342.88285e9),
     'SO65':('Xb4d','spw31',251.82577e9),'HCN32':('Xb4d','spw27',265.88618e9),'HCO+32':('Xb4d','spw25',267.55762e9)}
POS={'S':S0,'N':N0,'ridge':RED,'lobe':LOBE}
txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
fig,ax=plt.subplots(len(WIN),4,figsize=(16,2.4*len(WIN)),sharex=False)
for i,(nm,(o,sp,f0)) in enumerate(WIN.items()):
    c,v,ps,wc,x0,y0,b0=load_cube(o,sp,f0,-1500,1500,half=24.0)
    px=pxy_fn(wc,x0,y0); dv=abs(np.median(np.diff(v)))
    # robust noise per channel: MAD of off-source pixels (>8" from N-S axis midpoint) in line-free channels
    mid=offset_pos(S0,0.1,2.55); mx,my=px(*mid); Y,X=np.mgrid[:c.shape[1],:c.shape[2]]
    off=(np.hypot(X-mx,Y-my)*ps>9)&(np.hypot(X-mx,Y-my)*ps<14)
    P(f'== {nm} ({o} {sp}) beam {b0[0]:.2f}x{b0[1]:.2f} PA {b0[2]:.0f}, dv {dv:.2f} km/s, v range {v.min():.0f}..{v.max():.0f}, pix {ps:.3f}"')
    for j,(pn,pos) in enumerate(POS.items()):
        x,y=px(*pos); ix,iy=int(round(x)),int(round(y)); s=c[:,iy,ix]
        ax[i,j].step(v,s*1e3,lw=.7); ax[i,j].set_title(f'{nm} @ {pn}',fontsize=8); ax[i,j].axvspan(180,320,color='r',alpha=.15); ax[i,j].axvspan(-320,-180,color='b',alpha=.15)
        np.savetxt(f'{OUT}/s110a_{nm}_{pn}.csv',np.c_[v,s],header='v_kms Jy/beam',delimiter=',')
        P(f'  {pn}: peak pixel ({ix},{iy}); peak {s.max()*1e3:.1f} mJy/beam at v={v[np.argmax(s)]:.0f}')
    # noise per channel in all-channel-from-off pixels
    nz=np.array([np.std(c[k][off]) for k in range(len(v))]); P(f'  off-source rms/channel (9-14" annulus; median over chans) {np.median(nz)*1e3:.2f} mJy/beam')
    ax[i,3].plot(v,nz*1e3); ax[i,3].set_title('off-source rms (mJy/beam)',fontsize=8)
plt.tight_layout(); plt.savefig(f'{OUT}/s110a_spectra.png',dpi=70)
open(OUT+'/s110a.txt','w').write('\n'.join(txt))
