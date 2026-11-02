# Q17 step d: figure (PV along PA185 through S, HCN/HCO+ 3-2 Xb4d@1.05", r=0.75" slit-averaged) and mirrored spectra of the ridge (+3.3" N of S) and the blue lobe (-1.85"), plus the wing sub-window shape test.
# Reads out/q17/s104a_pv_*.npz, s104b_spectra.npz.  usage: python3 -I s104d_q17_figs.py
import numpy as np, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q17'; txt=[]; P=lambda s:(txt.append(s),print(s))
fig,ax=plt.subplots(2,2,figsize=(13,9))
for j,mol in enumerate(('HCN','HCO+')):
    d=np.load(f'{OUT}/s104a_pv_{mol}.npz'); v=d['v']; o=d['offs']; pv=d['pv'].T   # (v, off)
    e=np.std(pv[((v>-900)&(v<-450))])
    a=ax[0,j]; a.imshow(pv,origin='lower',aspect='auto',extent=[o[0]-.25,o[-1]+.25,v[0],v[-1]],cmap='viridis',vmin=-3*e,vmax=12*e)
    a.contour(o,v,pv,levels=[3*e,5*e,10*e,20*e],colors='w',linewidths=.5); a.set_ylim(-450,450); a.set_xlabel('offset N of S along PA185 (")'); a.set_ylabel('v - VSYS (km/s)'); a.set_title(f'{mol} 3-2 PV (1.05" beam, r=0.75" slit); contours 3,5,10,20 sigma ({e*1e3:.2f} mJy/b)')
    for x in (0,5.2): a.axvline(x,color='r',ls=':',lw=.7)
    s=np.load(f'{OUT}/s104b_spectra.npz'); b=ax[1,j]; vv=s[f'v_{mol}']
    b.step(vv,s[f'{mol}_ridge']*1e3,where='mid',label='ridge (+3.3" N of S)'); b.step(-vv,s[f'{mol}_lobe']*1e3,where='mid',label='blue lobe (-1.85"), velocity sign flipped'); b.set_xlim(-100,450); b.set_ylim(-1,5); b.set_xlabel('|v - VSYS| (km/s) (red side for the ridge; sign flipped for the lobe)'); b.set_ylabel('mJy/beam (mean over r=1")'); b.legend(); b.axhline(0,color='k',lw=.5)
    # sub-window shape test
    dv=abs(np.median(np.diff(vv)))
    def F(sp,a_,b_): k=(vv>=a_)&(vv<=b_); return sp[k].sum()*dv, s[f'{mol}_ridge_e']*np.sqrt(k.sum())*dv
    for nm,sp,sg in (('ridge',s[f'{mol}_ridge'],1),('lobe',s[f'{mol}_lobe'],-1)):
        w=[]
        for a_,b_ in ((120,200),(200,280),(280,360)):
            k=(vv*sg>=a_)&(vv*sg<=b_); w.append((sp[k].sum()*dv,s[f'{mol}_{nm}_e']*np.sqrt(k.sum())*dv))
        P(f'{mol} {nm}: wing flux (mean-over-aperture Jy/b km/s) 120-200: {w[0][0]:.3f}+-{w[0][1]:.3f}, 200-280: {w[1][0]:.3f}+-{w[1][1]:.3f}, 280-360: {w[2][0]:.3f}+-{w[2][1]:.3f}; (120-200)/(200-280)={w[0][0]/w[1][0]:.2f}+-{w[0][0]/w[1][0]*np.hypot(w[0][1]/w[0][0],w[1][1]/w[1][0]):.2f}')
        if nm=='ridge': wr=w
        else: wl=w
    ra=wr[0][0]/wr[1][0]; rl=wl[0][0]/wl[1][0]; er=np.hypot(wr[0][1]/wr[0][0],wr[1][1]/wr[1][0])*ra; el=np.hypot(wl[0][1]/wl[0][0],wl[1][1]/wl[1][0])*rl
    P(f'{mol} shape difference (120-200)/(200-280): ridge {ra:.2f}+-{er:.2f} vs lobe {rl:.2f}+-{el:.2f}: {abs(ra-rl)/np.hypot(er,el):.1f} sigma (white noise)')
plt.tight_layout(); plt.savefig(OUT+'/s104_pv_spectra.png',dpi=80); open(OUT+'/s104d.txt','w').write('\n'.join(txt))
