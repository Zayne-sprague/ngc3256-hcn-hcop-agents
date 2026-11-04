# Fig. 2: HCN/HCO+ 3-2 ratio vs radius (values: Table C2 / r09_profile.json). Usage: python3 -I fig2_profile.py
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
r=np.array([0.125,0.375,0.75,1.5,2.5]); 
N=np.array([.998,.877,.648,.432,.377]); Ne=np.array([.019,.017,.011,.017,.054])
S=np.array([.443,.386,.331,.340,.446]); Se=np.array([.014,.017,.021,.033,.027])
sysx=np.array([.03,.03,.03,.05,.05])   # array-to-array (Xb4b vs Xb4d) and missing-flux allowance
fig,ax=plt.subplots(figsize=(6.2,4.3))
for y,e,c,l,dx in ((N,Ne,'#D55E00','N',0.97),(S,Se,'#0072B2','S',1.03)):
    x=r*dx
    ax.errorbar(x,y,yerr=np.hypot(e,sysx),fmt='none',ecolor=c,alpha=.35,capsize=0,lw=5)
    ax.errorbar(x,y,yerr=e,fmt='o-',color=c,capsize=3,label=l)
ax.set_xscale('log'); ax.set_xlim(0.09,3.6); ax.set_ylim(0,1.2)
ax.set_xticks([.125,.375,.75,1.5,2.5]); ax.set_xticklabels(['0-0.25','0.25-0.5','0.5-1','1-2','2-3'])
ax.minorticks_off()
ax.set_xlabel('annulus radius (arcsec; 1" = 194 pc)'); ax.set_ylabel('HCN/HCO$^+$ J=3-2 $T_b$ ratio')
ax.legend(frameon=False,loc='upper right')
ax.set_title('Thin bars: statistical; thick bars: statistical + array/flux systematic',fontsize=8)
fig.tight_layout(); fig.savefig('../figs/fig2_profile.png',dpi=160)
