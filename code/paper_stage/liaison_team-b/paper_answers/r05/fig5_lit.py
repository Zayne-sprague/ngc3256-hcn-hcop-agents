"""R05 / Figure 5: HCN/HCO+ vs scale compared with Nishimura+2024 (3-2) and Izumi+2016 (4-3).
Literature values are transcribed in /workspace/scripts/s300_q46_lit_context.py (from lit.py full text: Nishimura+2024 A&A 686, A48 Tables 1,2,4; Izumi+2016 Table 1).
Our 3-2 values: out/q26/s171_summary.txt (Xb4b aperture ratios; 6.5" from Xb4d s170_Xb4d.json). Our 4-3: Xb47 values as in s300 (4.14x2.70", peak pixel, +-350 km/s).
usage: python3 -I fig5_lit.py"""
import numpy as np, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
exec(open('/workspace/scripts/s300_q46_lit_context.py').read().split('R=np.array(resolved)')[0].split('import numpy as np')[1].replace('from scipy.stats import spearmanr',''))
R=np.array(resolved); G=np.array(globl); A=np.array(aAGN); P=np.array(res_pc); ok=~np.isnan(R)
# our measurements (aperture R, stat err ~0.01; flux scale 14% not shown)
ours={'N':{55:.99,120:.84,250:.69},'S':{55:.43,120:.40,250:.37}}
oursx={'N':{55:55*0+0.2835*194,120:0.6186*194,250:1.2887*194}}
xs={55:0.2835*194,120:0.6186*194,250:1.2887*194}
d=json.load(open('/workspace/out/q26/s170_Xb4d.json'))['6.5"']; r65={k:d[k]['R_ap'][0] for k in ('N','S')}
def frac(x,sel=ok): a=R[sel]; return (a<x).mean(),(a<=x).mean()
lines=['Percentile = fraction of Nishimura+2024 resolved nuclear ratios BELOW the value (n=%d)'%ok.sum()]
for k,x in (('N',.99),('N',.84),('N',.69),('S',.43),('S',.40),('S',.37)):
    pass
rows=[]
for k in ('N','S'):
    for sc,v in ours[k].items():
        lo,hi=frac(v); m=ok&(abs(P-xs[sc])<0.5*xs[sc]) if True else ok
        mlo,mhi=frac(v,m); rows.append('%s %3d pc R=%.2f: all resolved below %.0f%% (<= %.0f%%); within +-50%% in resolution (n=%d) below %.0f%%'%(k,sc,v,100*lo,100*hi,m.sum(),100*mlo))
lines+=rows
for k,v in r65.items(): lines.append('%s 6.5" (1.26 kpc, Xb4d) R=%.2f: below %.0f%% of resolved; global-ratio sample (n=%d) below %.0f%%'%(k,v,100*frac(v)[0],(~np.isnan(G)).sum(),100*(G[~np.isnan(G)]<v).mean()))
lines.append('Nishimura resolved: median %.2f, range %.2f-%.2f; AGN-dominated (aAGN>=0.5) only IRAS F05189-2524 (aAGN 0.60, R=1.63)'%(np.nanmedian(R),np.nanmin(R),np.nanmax(R)))
open('percentiles.txt','w').write('\n'.join(lines)+'\n'); print('\n'.join(lines))
plt.rcParams.update({'font.size':10})
fig,(a,b)=plt.subplots(1,2,figsize=(12,5),gridspec_kw=dict(width_ratios=[1.35,1]))
cl={'lo':'#999999','mid':'#56B4E9','hi':'#CC79A7'}
def grp(x): return 'hi' if x>=0.25 else ('mid' if x>=0.10 else 'lo')
for i in np.where(ok)[0]:
    g='nan' if np.isnan(A[i]) else grp(A[i])
    a.scatter(P[i],R[i],s=34,marker='o' if g!='hi' else 'D',facecolor=cl.get(g,'none'),edgecolor='k' if g!='nan' else '0.4',lw=0.6,zorder=2)
a.axhspan(1.1,3.8,color='0.93',zorder=0); a.axhline(1,color='0.6',lw=.6,ls=':')
for k,c,m in (('N','#D55E00','o'),('S','#0072B2','o')):
    x=[xs[s] for s in ours[k]]; y=list(ours[k].values()); a.plot(x,y,'-',color=c,lw=1.5,zorder=4)
    a.scatter(x,y,s=70,marker='*',color=c,edgecolor='k',lw=.6,zorder=5)
    a.scatter([6.5*194],[r65[k]],s=70,marker='*',facecolor='w',edgecolor=c,lw=1.5,zorder=5)
    a.annotate(f'NGC 3256 {k}',(xs[55],ours[k][55]),xytext=(-6,-16 if k=='N' else -14),textcoords='offset points',color=c,ha='right',fontsize=9,weight='bold',annotation_clip=False)
a.set_xscale('log'); a.set_xlim(35,1700); a.set_ylim(0,3.8); a.set_xlabel('spatial scale (pc): beam FWHM (ours) / nuclear size, Nishimura+24 Table 4'); a.set_ylabel(r'HCN/HCO$^+$ J=3–2 brightness-temperature ratio')
a.text(1650,1.12,'HCN/HCO$^+$ > 1 (AGN-like)',ha='right',fontsize=8,color='0.4')
h=[Line2D([],[],marker='o',ls='',mfc=cl['lo'],mec='k',label=r'Nishimura+24 resolved, $\alpha_{\rm AGN}<0.10$'),Line2D([],[],marker='o',ls='',mfc=cl['mid'],mec='k',label=r'0.10–0.25'),Line2D([],[],marker='D',ls='',mfc=cl['hi'],mec='k',label=r'$\geq0.25$ (IRAS F05189: 0.60)'),
   Line2D([],[],marker='o',ls='',mfc='none',mec='0.4',label=r'$\alpha_{\rm AGN}$ not listed (UGC 2982, NGC 5734 have no resolved ratio)'),
   Line2D([],[],marker='*',ls='-',color='k',ms=10,label='NGC 3256 (this work): 55, 120, 250 pc'),Line2D([],[],marker='*',ls='',mfc='w',mec='k',ms=10,label='6.5" (1.3 kpc, Xb4d, flux-loss affected)')]
a.legend(handles=h,fontsize=7.5,loc='upper left',frameon=False); a.set_title('(a) vs spatial scale',fontsize=10)
# panel b
iz={"AGN (Seyfert), hi-res":[1.53,2.84,3.19,3.14,2.58,2.01,1.11],"starburst, hi-res":[.41,1.03,.82,.24,.75,.45,.48]}
b.axvspan(0.24,1.03,color='#009E73',alpha=0.25,lw=0); b.axvspan(1.11,3.19,color='#D55E00',alpha=0.2,lw=0)
b.text(0.63,4.55,'SB\n0.24–1.03',ha='center',fontsize=9,color='#006e50'); b.text(2.15,4.55,'AGN\n1.1–3.2',ha='center',fontsize=9,color='#a03d00')
rng=np.random.default_rng(1)
for k,(lab,v) in enumerate(iz.items()): b.scatter(v,[3.6+0.0*np.array(v)+rng.uniform(-.12,.12,len(v))],s=18,color='0.35',zorder=3,label='Izumi+16 hi-res individual' if k==0 else None)
b.scatter([0.52],[2.2],s=130,marker='*',color='#D55E00',edgecolor='k',zorder=5); b.scatter([0.29],[1.6],s=130,marker='*',color='#0072B2',edgecolor='k',zorder=5)
b.scatter([0.38],[0.9],s=60,marker='s',color='k',zorder=5); b.errorbar([0.39],[0.4],xerr=[0.11],fmt='o',color='0.3',zorder=5)
b.text(0.60,2.2,'N: 0.52',va='center',fontsize=9); b.text(0.37,1.6,'S: 0.29',va='center',fontsize=9)
b.text(0.47,0.9,'N+S, 10" (Q38r2): 0.38',va='center',fontsize=8); b.text(0.55,0.4,'Izumi+16 APEX NGC 3256: 0.39±0.11',va='center',fontsize=8)
b.text(3.45,1.3,'ours: ALMA Xb47 4-3,\nbeam 4.1"×2.7"\n(800×520 pc),\npeak pixel, ±350 km/s',ha='right',fontsize=8)
b.set_yticks([]); b.set_ylim(0,5.3); b.set_xlim(0,3.5); b.set_xlabel(r'HCN/HCO$^+$ J=4–3 brightness-temperature ratio'); b.set_title('(b) on the Izumi+16 J=4-3 diagnostic bands',fontsize=10)
b.text(1.8,3.15,'Izumi+16 high-resolution sources (individual)',ha='center',fontsize=8,color='0.3')
fig.tight_layout(); fig.savefig('fig5_ratio_vs_scale_and_izumi.png',dpi=180)
