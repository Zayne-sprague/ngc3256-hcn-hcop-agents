"""Figure 6: HCN v2=1f 3-2 non-detection at N and S (Xb4b 48 pc beam). Reuses the s42_hcn_vib.py method
(/snapshots/team-b-s29/scripts) on the cached region spectra /snapshots/team-b-s15/out/regspec/ (pbcor, peak pixel)."""
import sys,warnings; warnings.filterwarnings('ignore')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
C=299792.458; VSYS=2780.; R='/snapshots/team-b-s15/out/regspec/'
Fv,Fh,Fo=267.199283,265.886434,267.557626
def vel(nu,f0): return (1-nu/f0)*C-VSYS
zv=np.load(R+'Xb4b.spw25.cube.pbcor.npz'); zh=np.load(R+'Xb4b.spw27.cube.pbcor.npz')
nuv=zv['nu']/1e9; nuh=zh['nu']/1e9
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(11,4),sharey=False); res={}
for ax,reg in zip(axs,('N','S')):
    yv=zv[f'{reg}_pix']*1e3; yh=zh[f'{reg}_pix']*1e3
    vh=vel(nuh,Fh); o=np.argsort(vh); vh,yhs=vh[o],yh[o]
    yhs=yhs-np.median(yhs[(abs(vh)>500)&(abs(vh)<1000)]); msk=abs(vh)<600
    I=np.trapezoid(yhs[msk],vh[msk]); Ti=interp1d(vh,yhs,bounds_error=False,fill_value=0.)
    vv=vel(nuv,Fv); vo=vel(nuv,Fo); u=vo/100.
    m=(vv>-250)&(vv<250)&np.isfinite(yv); lf=(abs(vo)>600)&(abs(vv)>450)&np.isfinite(yv); sg=np.std(yv[lf])
    t=Ti(vv[m]); y=yv[m]; um=u[m]
    fn=lambda p:np.exp(p[0]+p[1]*um)
    def prof(al,starts=([0,-.5],[-1,-1],[1,-.2])):
        b=None
        for p0 in starts:
            r=least_squares(lambda p:(y-al*t-fn(p))/sg,p0,max_nfev=500)
            if b is None or r.cost<b[0]: b=(r.cost,r.x)
        return 2*b[0],b[1]
    als=np.linspace(0,0.2,81); ch=np.array([prof(a)[0] for a in als]); i=ch.argmin(); dof=m.sum()-3; s2=max(ch[i]/dof,1.)
    d=(ch-ch[i])/s2; ok=np.where((d>=9)&(als>=als[i]))[0]; ul=als[ok[0]]; res[reg]=(als[i],ul,I)
    _,p=prof(als[i])
    # HCO+ wing model over full window, from best-fit
    f=nuv/(1-VSYS/C); ax.step(f,yv,where='mid',c='0.3',lw=1,label='observed (peak pixel)')
    mm=np.isfinite(yv)&(vo>-900)&(vo<2000)
    ax.plot(f[m],np.exp(p[0]+p[1]*um),c='#0072B2',lw=1.8,label=f'HCO$^+$ wing model (best, exp-linear)')
    tt=Ti(vel(nuv,Fv)); ax.plot(f[m],np.exp(p[0]+p[1]*um)+ul*tt[m],c='#D55E00',lw=1.5,label=rf'+ HCN $v=0$ template × $\alpha_{{3\sigma}}$={ul:.3f}')
    ax.axvline(Fv,c='#009E73',ls=':'); ax.text(Fv-0.005,-0.95,'HCN v$_2$=1f 3-2\n267.199 GHz',color='#009E73',fontsize=8,ha='right')
    ax.axvline(Fo,c='gray',ls=':'); ax.set_xlim(266.85,267.9)
    # HCO+ peak is far above; clip y to wing level
    ax.set_ylim(-1.2,2.4); ax.set_xlabel('rest-frame frequency (GHz), v$_{sys}$=2780 km/s'); ax.set_ylabel('Jy/beam ×10$^{-3}$ (mJy/beam), pbcor')
    ax.set_title(f'{reg}: 0.24" (48 pc) beam; $\\sigma$={sg:.2f} mJy/beam / 3.9 MHz ch',fontsize=9,loc='left'); ax.text(Fo-0.01,-0.95,'HCO$^+$ 3-2 (267.558)',ha='right',fontsize=8,color='gray')
    ax.legend(frameon=False,fontsize=7,loc='upper left',bbox_to_anchor=(0.0,1.0))
fig.tight_layout(); fig.savefig('hcn_vib_limits.png',dpi=180)
print(res)
