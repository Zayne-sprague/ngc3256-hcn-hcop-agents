# Fig. 3: N/S contrast of HCN/HCO+ versus beam. Values: Tables B1, B1b, B1c of the paper. Usage: python3 -I fig3_contrast_beam.py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
B,G,O='#0072B2','#009E73','#D55E00'
fig,ax=plt.subplots(figsize=(7.2,4.8))
def pts(x,y,e,**k): ax.errorbar(x,y,yerr=e,capsize=2.5,linestyle='none',**k)
# 3-2, Xb4b convolved in the image plane (12-m only)
bx=[0.3,0.6,1.05,1.5,2.2]
ax.errorbar(bx,[2.26,2.21,2.09,2.01,1.89],yerr=[.11,.11,.10,.08,.07],fmt='o-',mfc='w',color=B,capsize=2.5,label='3-2, 12-m, annulus r<0.25" (convolved)')
ax.errorbar([x*1.03 for x in bx],[2.25,2.11,1.98,1.85,1.66],yerr=[.13,.11,.08,.07,.07],fmt='v--',mfc='w',color=B,alpha=.6,capsize=2.5,label='3-2, 12-m, aperture r = beam (convolved)')
# 3-2 flux-complete 12m+7m, uv-tapered
pts([1.5*0.97,2.2*0.97],[1.73,1.67],[.03,.03],fmt='o',color=B,label='3-2, 12-m+7-m, aperture r=0.5" (uv taper)')
pts([1.5*0.94,2.2*0.94],[1.53,1.44],[.03,.02],fmt='v',color=B,alpha=.6,label='3-2, 12-m+7-m, aperture r = beam (uv taper)')
# 1-0 flux-complete
pts([1.5*1.06,2.2*1.06],[1.77,1.62],[.10,.06],fmt='s',color=G,label='1-0, flux-complete, r=0.5"')
pts([1.5*1.09,2.2*1.09],[1.62,1.48],[.06,.04],fmt='v',color=G,alpha=.6,label='1-0, flux-complete, r = beam')
# ACA only, matched recipe
pts([5.0,6.5],[1.36,1.23],[.18,.16],fmt='o',color=B,mfc='w',label='3-2, ACA only')
pts([5.0*1.04,6.5*1.04],[1.45,1.25],[.13,.14],fmt='o',color=O,label='4-3, ACA only')
ax.errorbar([5.0*1.1],[1.8],yerr=[[.23],[.32]],xerr=[[1.0],[2.0]],fmt='s',color=O,mfc='w',capsize=2.5,label='4-3, ACA visibility fit (no beam)')
ax.axhspan(1.46,1.52,color='0.8',zorder=0); ax.text(0.27,1.54,'Harada et al. 1.7", 1-0: 1.49',fontsize=7.5,color='0.35')
ax.axhline(1,color='k',lw=.8)
ax.set_xscale('log'); ax.set_xlim(0.22,9); ax.set_ylim(0.9,2.55)
ax.set_xticks([.25,.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8']); ax.minorticks_off()
ax.set_xlabel('beam FWHM (arcsec; 1" = 194 pc)'); ax.set_ylabel('N/S contrast of HCN/HCO$^+$')
ax.legend(fontsize=7,frameon=False,loc='upper right',ncol=1)
fig.tight_layout(); fig.savefig('../figs/fig3_contrast_beam.png',dpi=160)
