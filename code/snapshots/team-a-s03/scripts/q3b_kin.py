# Q3 round 2, fix 3: outflow vs disk/tidal feature from kinematics.  PV diagrams along the N-S axis (x=-0.2") and parallel cuts, E-W rotation curve at the S and Soutflow declinations, centroid/dispersion/red-edge vs position,
# beam-smearing prediction of a rotating disk, Xb4b (0.5") moment-1 gradient in the lobe.  Outputs: out/q3b/kin_*.csv, kin.log, q3b_pv.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
import pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
log=open(OUT+'kin.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
CB={sp:vcube(sp,'Xb4d',1.05,9.) for sp in ('HCN','HCO+')}
def strip_spec(c,x0,y0,halfw_x=None,halfw_y=None):
    m=np.ones(c['dx'].shape,bool)
    if halfw_x is not None: m&=abs(c['dx']-x0)<=halfw_x
    if halfw_y is not None: m&=abs(c['dy']-y0)<=halfw_y
    return c['a'][:,m].mean(1)*1e3
def moments(v,s,e,lo=2150,hi=3150,thr=3.):
    """first moment/dispersion of the channels in lo..hi after 30 km/s smoothing, using only channels >thr*e_smooth; red/blue edge = velocity where the smoothed profile falls to 25% of its peak on the red/blue side"""
    n=7; sm=np.convolve(s,np.ones(n)/n,'same'); es=e/np.sqrt(n)
    k=(v>=lo)&(v<=hi); w=np.where(k&(sm>thr*es),sm,0.)
    if w.sum()<=0: return dict(vm=np.nan,sd=np.nan,vred=np.nan,vblue=np.nan,I=np.nan,peak=np.nan)
    vm=(w*v).sum()/w.sum(); sd=np.sqrt((w*(v-vm)**2).sum()/w.sum()); pk=sm[k].max(); vp=v[k][np.argmax(sm[k])]
    kk=k&(v>vp); red=v[kk][np.where(sm[kk]<0.25*pk)[0][0]] if (sm[kk]<0.25*pk).any() else np.nan   # v descending: first channel below 25% going to larger v is wrong -> handled below
    return dict(vm=vm,sd=sd,vpeak=vp,pk=pk,I=w.sum()*abs(np.median(np.diff(v)))/1e3)
def edges(v,s,e):
    n=7; sm=np.convolve(s,np.ones(n)/n,'same'); k=(v>=2100)&(v<=3200); o=np.argsort(v); vv=v[o]; ss=sm[o]; kk=(vv>=2100)&(vv<=3200)
    vv=vv[kk]; ss=ss[kk]; ip=np.argmax(ss); pk=ss[ip]
    if pk<4*e/np.sqrt(n): return np.nan,np.nan,np.nan
    r=np.where(ss[ip:]<0.25*pk)[0]; b=np.where(ss[:ip+1][::-1]<0.25*pk)[0]
    return (vv[ip+r[0]] if len(r) else np.nan),(vv[ip-b[0]] if len(b) else np.nan),vv[ip]
rows=[]
ys=np.arange(-7.2,3.01,0.45)
# ---- cuts along y at several x, strips 1" wide in x, 0.45" in y
for sp in ('HCN','HCO+'):
    c=CB[sp]
    for x0 in (-4.0,-2.0,-0.2,2.0,4.0):
        for y0 in ys:
            s=strip_spec(c,x0,y0,0.5,0.225); s,e=baseline(c['v'],s,FREE[sp],1); m=moments(c['v'],s,e); red,blue,vp=edges(c['v'],s,e)
            rows.append(dict(line=sp,x=x0,y=y0,vm=m['vm'],sd=m['sd'],vred25=red,vblue25=blue,vpeak=vp,I=m['I'],noise=e))
K=pd.DataFrame(rows); K.round(2).to_csv(OUT+'kin_ns_cuts.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_rows',500)
pr('== N-S cut at x=-0.2" (1" wide), per 0.45" step: centroid, dispersion (flux-weighted, 3sigma-clipped, 30 km/s smoothed), red edge (25% of peak), blue edge, peak velocity')
pr(K[(K.x==-0.2)].pivot_table(index='y',columns='line',values=['vm','sd','vred25','vblue25','vpeak']).round(0).to_string())
# ---- E-W cuts at S and Soutflow dec (1" wide in y)
rows=[]
for sp in ('HCN','HCO+'):
    c=CB[sp]
    for y0,nm in ((-2.55,'S dec'),(-4.5,'Soutflow dec'),(-3.45,'mid')):
        for x0 in np.arange(-7,7.01,0.5):
            s=strip_spec(c,x0,y0,0.25,0.5); s,e=baseline(c['v'],s,FREE[sp],1); m=moments(c['v'],s,e); red,blue,vp=edges(c['v'],s,e)
            rows.append(dict(line=sp,cut=nm,y=y0,x=x0,vm=m['vm'],sd=m['sd'],vred25=red,vblue25=blue,vpeak=vp,I=m['I']))
EW=pd.DataFrame(rows); EW.round(1).to_csv(OUT+'kin_ew_cuts.csv',index=False)
pr('== E-W cuts (HCN): centroid / dispersion / peak velocity vs x (E positive)')
pr(EW[(EW.line=='HCN')].pivot_table(index='x',columns='cut',values=['vm','sd','vpeak']).round(0).dropna(how='all').to_string())
# ---- disk-rotation model from the E-W cut at the S dec: linear gradient of the HCN peak velocity within |x|<=2"
h=EW[(EW.line=='HCN')&(EW.cut=='S dec')&(abs(EW.x)<=2.0)&EW.vpeak.notna()]
p=np.polyfit(h.x,h.vpeak,1); pr(f'rotation gradient at the S dec (HCN peak velocity, |x|<=2"): {p[0]:.0f} km/s/arcsec, v(x=0)={p[1]:.0f} km/s')
h2=EW[(EW.line=='HCO+')&(EW.cut=='S dec')&(abs(EW.x)<=2.0)&EW.vpeak.notna()]; p2=np.polyfit(h2.x,h2.vpeak,1); pr(f'   HCO+: {p2[0]:.0f} km/s/arcsec, v(0)={p2[1]:.0f}')
g=abs(p[0]); sig_smear=g*1.05/2.355
pr(f'beam smearing of that gradient over the 1.05" beam adds sigma ~ {sig_smear:.0f} km/s in quadrature to the intrinsic ~60 km/s: expected sigma {np.hypot(60,sig_smear):.0f} km/s at Soutflow if it were disk gas on the minor axis; observed sigma 140-150 km/s (HCN, N-S cut table)')
# ---- Soutflow spectrum vs disk-model expectation: centroid offset from rotation-curve value at that x, red edge
for sp in ('HCN','HCO+'):
    c=CB[sp]; 
    for nm,(x,y,w) in {'Soutflow centroid-ap (-0.7,-4.8)':(-0.72,-4.79,1.0),'Soutflow peak (-0.15,-4.5)':(-0.15,-4.5,1.0),'S nucleus':(-0.15,-2.55,1.0),'W control (-4.4,-4.5)':(-4.4,-4.5,1.0),'E control (+4.0,-4.5)':(4.0,-4.5,1.0),'mirror N of S':(-0.09,-0.62,1.0)}.items():
        s=aper(c,x,y,w)*1e3; s,e=baseline(c['v'],s,FREE[sp],1); m=moments(c['v'],s,e); red,blue,vp=edges(c['v'],s,e)
        pr(f'{sp:5s} {nm:34s} centroid {m["vm"]:.0f} sd {m["sd"]:.0f} peak v {vp:.0f} red25 {red:.0f} blue25 {blue:.0f}  (rotation-curve value at that x: {p[1]+p[0]*x:.0f})')
# ---- Xb4b (0.5") HCN and HCO+ moment-1 in the S lobe: gradients
from scipy.stats import linregress
for sp in ('HCN','HCO+'):
    c=vcube(sp,'Xb4b',0.5,crop=7.5); v=c['v']
    m0=((v>=2300)&(v<=2900)); I=c['a'][m0].sum(0)*c['dv']; V=(c['a'][m0]*v[m0,None,None]).sum(0)*c['dv']/np.where(I>0,I,np.nan)
    # noise of m0 from the line-free channels
    nn=int(m0.sum()); fm=np.zeros(len(v),bool)
    for a,b in FREE[sp]: fm|=(v>=a)&(v<=b)
    cs=np.cumsum(np.concatenate([np.zeros((1,)+c['a'].shape[1:]),c['a']]),0); S=[(cs[i+nn]-cs[i])*c['dv'] for i in range(0,len(v)-nn,10) if fm[i] and fm[i+nn-1]]
    sg=1.4826*np.median(np.abs(np.array(S)-np.median(S,0)),0)
    sel=(I>5*sg)&(c['dy']<-3.2)&(c['dy']>-6.5)&(abs(c['dx'])<3)&np.isfinite(V)
    ly=linregress(c['dy'][sel],V[sel]); lx=linregress(c['dx'][sel],V[sel])
    pr(f'Xb4b {sp} 0.5" beam, S-lobe pixels (>5sigma in 2300-2900 m0, {sel.sum()} px = {sel.sum()*c["pix"]**2:.1f} arcsec^2): dV/dDec = {ly.slope:.1f}+-{ly.stderr:.1f} km/s/arcsec (r={ly.rvalue:.2f}); dV/dRA(E+) = {lx.slope:.1f}+-{lx.stderr:.1f} km/s/arcsec (r={lx.rvalue:.2f}); mean V {np.mean(V[sel]):.0f}, rms scatter {np.std(V[sel]):.0f} km/s')
# ---- figure
fig,axs=plt.subplots(2,5,figsize=(17,7),sharey=True)
for i,sp in enumerate(('HCN','HCO+')):
    c=CB[sp]
    for j,x0 in enumerate((-4.0,-2.0,-0.2,2.0,4.0)):
        pv=[]
        for y0 in ys:
            s=strip_spec(c,x0,y0,0.5,0.225); s,e=baseline(c['v'],s,FREE[sp],1); pv.append(s)
        pv=np.array(pv); ax=axs[i,j]; sel=(c['v']>2150)&(c['v']<3150)
        im=ax.imshow(pv[:,sel][:,::-1].T,origin='lower',aspect='auto',extent=[ys[0]-0.225,ys[-1]+0.225,c['v'][sel].min(),c['v'][sel].max()],vmin=0,vmax=np.percentile(pv[:,sel],99.5) if j==2 else 6,cmap='viridis')
        ax.axvline(-2.55,color='w',lw=.6,ls=':'); ax.axvline(2.55,color='w',lw=.6,ls=':'); ax.axvline(-4.5,color='r',lw=.6,ls=':')
        ax.set_title(f'{sp}  x={x0:+.1f}"'); 
        if j==0: ax.set_ylabel('v (km/s)')
        ax.set_xlabel('Dec offset (")')
plt.tight_layout(); plt.savefig(OUT+'q3b_pv.png',dpi=110)
