"""Q26 analysis: tables, N/S contrast significance, scale dependence, estimator/baseline/position systematics, wing vs core, figure.
usage: python3 -I s171_q26_analysis.py   (needs out/q26/s170_*.json, maps_*.npz)"""
import json, numpy as np, sys
O='/workspace/out/q26/'
A=json.load(open(O+'s170_Xb4b.json')); B=json.load(open(O+'s170_Xb4b_blueonly.json')); D=json.load(open(O+'s170_Xb4d.json')); F=json.load(open(O+'s170_Xb4f.json'))
rng=np.random.default_rng(26)
txt=[]
def P(s=''): txt.append(s); print(s)
SCS=list(A.keys()); TH={'0.26"(50pc)':0.26,'55pc(0.2835")':0.2835,'120pc(0.6186")':0.6186,'250pc(1.2887")':1.2887}
SYS=0.10   # flux scale per line
# ---- Table 1
P('TABLE 1  HCN/HCO+ 3-2 T_b ratio (Xb4b). R_ap = aperture (r = target beam FWHM) mean of m0 at HCO+ peak [window N,S,Bridge +-350 km/s about 2780; others vc+-250]; R_pix = same pixel; R_pk = peak m0 / peak m0 (own peaks); R_ch = 3-channel-smoothed peak T_b ratio')
P('stat errors from line-free channels; total adds 14%% flux-scale (2 x %d%%) in quadrature; spread = max-min of 4 estimators' % (SYS*100))
P('%-16s %-4s %-6s %-6s %-10s %-10s %-8s %-8s %-6s %-8s | TbHCN TbHCO+ (K, peak chan) | snr_ap(HCN,HCO+)'%('scale','reg','R_ap','e_st','R_pix','e_st','R_pk','R_ch','spread','e_tot'))
for sn in SCS:
    for rn,o in A[sn].items():
        ests=[o['R_ap'][0],o['R_pix'][0],o['R_peakpeak'][0],o['R_chanpeak'][0]]
        e=o['R_ap'][1]; etot=np.hypot(e,0.14*o['R_ap'][0]) if np.isfinite(e) else float('nan')
        P('%-16s %-6s %-6.2f %-6.2f %-6.2f %-6.2f %-8.2f %-8.2f %-6.2f %-8.2f | %.2f %.2f | %.0f %.0f'%(sn,rn,*o['R_ap'],*o['R_pix'],ests[2],ests[3],max(ests)-min(ests),etot,o['HCN']['Tchan_peak'],o['HCO+']['Tchan_peak'],*o['snr_ap']))
# ---- N/S contrast
P('\nTABLE 2  N/S contrast (flux scale cancels in the contrast; stat errors only, MC 20000)')
def mc(x,ex,y,ey,n=20000):
    a=rng.normal(x,ex,n); b=rng.normal(y,ey,n); r=a/b; return np.median(r),np.percentile(r,16),np.percentile(r,84)
CON={}
for sn in SCS:
    N,S,Br=A[sn]['N'],A[sn]['S'],A[sn]['Bridge']
    for est,key in (('aperture','R_ap'),('pixel','R_pix')):
        m,lo,hi=mc(*N[key],*S[key]); P('%-16s %-8s R_N/R_S = %.2f [%.2f,%.2f]   (R_N %.2f R_S %.2f; N/Bridge %.2f)'%(sn,est,m,lo,hi,N[key][0],S[key][0],N[key][0]/Br[key][0]))
        CON[(sn,est)]=(m,lo,hi)
    # peak T_b per line N/S
    P('%-16s Tb_pk(chan) N/S: HCO+ %.2f HCN %.2f ; peak-chan ratio contrast %.2f'%(sn,N['HCO+']['Tchan_peak']/S['HCO+']['Tchan_peak'],N['HCN']['Tchan_peak']/S['HCN']['Tchan_peak'],N['R_chanpeak'][0]/S['R_chanpeak'][0]))
# Xb4d / Xb4f
P('\nTABLE 3  consistency with coarser arrays (same convolution code; circular beams)')
for nm,X,sc in (('Xb4b',A,'250pc(1.2887")'),('Xb4d',D,'250pc(1.2887")'),('Xb4d',D,'6.5"'),('Xb4f',F,'6.5"')):
    P('%s @ %s: '%(nm,sc)+'  '.join('%s %.2f+-%.2f'%(r,X[sc][r]['R_ap'][0],X[sc][r]['R_ap'][1]) for r in ('N','S','Bridge','E1','E2','Ridge','Slobe') if r in X[sc]))
# ---- blue-only baseline
P('\nTABLE 4  baseline systematics: line-free baseline blue-only (v_rel<-450) vs both sides: R_ap')
for sn in SCS:
    P('%-16s '%sn+'  '.join('%s %.2f/%.2f'%(r,A[sn][r]['R_ap'][0],B[sn][r]['R_ap'][0]) for r in ('N','S','Bridge','E1','E2','Sblue')))
# ---- scale dependence
P('\nTABLE 5  scale dependence: power-law R ~ theta^alpha over 0.26-1.29" (aperture, 4 scales; 0.26 and 0.2835 near-duplicates -> use 3: 0.2835, 0.62, 1.29), alpha +- from stat err')
for rn in ('N','S','Bridge','E1','E2','Sblue'):
    xs=np.log([TH[s] for s in ('55pc(0.2835")','120pc(0.6186")','250pc(1.2887")')]); ys=[A[s][rn]['R_ap'] for s in ('55pc(0.2835")','120pc(0.6186")','250pc(1.2887")')]
    y=np.log([a for a,b in ys]); w=np.array([ (b/a) for a,b in ys]); w=np.where(np.isfinite(w),w,0.2)
    c=np.polyfit(xs,y,1,w=1/np.hypot(w,0.03)); cov=np.linalg.inv(np.vstack([xs,np.ones(3)])@np.diag(1/np.hypot(w,0.03)**2)@np.vstack([xs,np.ones(3)]).T)
    P('%-7s alpha = %+.2f +- %.2f   R(55)=%.2f R(250)=%.2f  (R55/R250 = %.2f)'%(rn,c[0],np.sqrt(cov[0,0]),ys[0][0],ys[2][0],ys[0][0]/ys[2][0]))
# ---- wing vs core
P('\nTABLE 6  velocity-resolved: wing/core ratio (aperture, windows relative to the local HCO+ centroid vc: blue -300..-100, core +-100, red +100..+300 km/s)')
for sn in ('55pc(0.2835")','250pc(1.2887")'):
    for rn in ('N','S','Bridge','E1','E2','Sblue'):
        cw=A[sn][rn]['cw']; s=[]
        for k in ('blue','core','red'):
            c=cw[k]; s.append('%s %s'%(k,'%.2f+-%.2f'%(c['R'],c['eR']) if c['R'] is not None and c['eR']==c['eR'] else 'n/a (HCN %.2f+-%.2f, HCO+ %.2f+-%.2f)'%(c['HCN'],c['eHCN'],c['HCO'],c['eHCO'])))
        P('%-16s %-7s vc %+4.0f | '%(sn,rn,A[sn][rn]['vc'])+' | '.join(s))
P('\nTABLE 7  50 km/s channel-bin ratio (HCN/HCO+ T_b ratio, bins centred rel. to vc, only bins with both >3 sigma) at 55 pc and 250 pc')
for sn in ('55pc(0.2835")','250pc(1.2887")'):
    for rn in ('N','S','Bridge'):
        r=[]
        for c,a,ea,b,eb in A[sn][rn]['bins']:
            r.append('%+d:%s'%(c,'%.2f'%(a/b) if a>3*ea and b>3*eb else '--'))
        P('%-16s %-7s '%(sn,rn)+' '.join(r))
# special windows
P('\nTABLE 8  fixed windows rel. to 2780: Slobe (-320..-180) and Ridge (+180..+320) in Xb4b and Xb4d at 250 pc')
for nm,X in (('Xb4b',A),('Xb4d',D)):
    for sn in ('250pc(1.2887")',):
        for rn in ('Slobe','Ridge'):
            o=X[sn][rn]; P('%s %s %s R_ap %.2f+-%.2f R_pix %.2f+-%.2f snr_ap HCN %.1f HCO+ %.1f'%(nm,sn,rn,*o['R_ap'],*o['R_pix'],*o['snr_ap']))
# ---- map-based position robustness: R in 3x3 pixel neighbourhoods
P('\nTABLE 9  map-based (5-sigma both) R in N and S at the finest scales: median and 16-84 range over the pixels within 0.15" of the HCO+ peak position')
def pixpos(sn,rn): return A[sn][rn]['pos_pix']
for sn in ('0.26"(50pc)','250pc(1.2887")'):
    tag=sn.replace('"','').replace('(','_').replace(')','')
    m=np.load(O+'maps_%s.npz'%tag); R=m['R']; pix=float(m['pix']); ny,nx=R.shape; Y,X=np.mgrid[:ny,:nx]
    for rn in ('N','S'):
        px,py=pixpos(sn,rn); sel=np.hypot(X-px,Y-py)*pix<=0.15; v=R[sel]; v=v[np.isfinite(v)]
        P('%-14s %s: n=%d median %.2f [%.2f,%.2f]'%(sn,rn,len(v),np.median(v),*np.percentile(v,[16,84])))
open(O+'s171_summary.txt','w').write('\n'.join(txt))
# ---- figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
cols={'N':'#d55e00','S':'#0072b2','Bridge':'#009e73','E1':'#cc79a7','E2':'#999999','Sblue':'#e69f00'}
for rn,c in cols.items():
    xs=[TH[s] for s in TH]; ys=[A[s][rn]['R_ap'][0] for s in TH]; es=[A[s][rn]['R_ap'][1] for s in TH]
    ax[0].errorbar(np.array(xs)*194,ys,es,color=c,marker='o',label=rn)
ax[0].set_xscale('log'); ax[0].set_xlabel('beam FWHM (pc)'); ax[0].set_ylabel('T_b(HCN 3-2)/T_b(HCO+ 3-2)'); ax[0].legend(fontsize=7); ax[0].set_title('aperture ratio vs scale (Xb4b)')
for k,(sn,tt) in enumerate((('0.26"(50pc)','0.26" (50 pc)'),('250pc(1.2887")','1.29" (250 pc)'))):
    tag=sn.replace('"','').replace('(','_').replace(')','')
    m=np.load(O+'maps_%s.npz'%tag); R=m['R']; pix=float(m['pix']); h=R.shape[0]//2
    sl=(slice(None),slice(None)); a=ax[1+k]; ext=None
    im=a.imshow(R,origin='lower',vmin=0,vmax=1.4,cmap='viridis'); a.set_title('HCN/HCO+ T_b ratio, %s, >5 sigma both'%tt)
    for rn in ('N','S','Bridge','E1','E2'):
        px,py=pixpos(sn,rn); a.plot(px,py,'r+'); a.text(px,py,rn,color='w',fontsize=7)
    plt.colorbar(im,ax=a,fraction=0.04)
plt.tight_layout(); plt.savefig(O+'q26_ratio_figure.png',dpi=110)
