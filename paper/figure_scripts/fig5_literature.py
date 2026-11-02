# Fig. 5: literature context. Literature values: Nishimura+2024 Tables 1,2,4; Izumi+2016 Table 1 (transcribed in team script s300_q46_lit_context.py).
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=np.array([1.738,.793,2.671,1.323,1.630,1.320,1.789,.445,.691,.629,1.340,1.132,.728,1.490,1.076,1.028,.771,.564,np.nan,1.489,.666,2.344,np.nan])
P=np.array([298,382,225,388,343,281,423,344,161,139,102,51,397,142,126,83,95,63,117,59,101,62,89])
A=np.array([.05,.19,.04,.05,.60,.04,.08,.19,.11,.07,.03,.03,np.nan,.10,.07,.24,.07,.10,.11,.10,.06,.03,.14])
ok=~np.isnan(R)
# ours: ring median within 0.13" radius (aperture diameter = one beam, as in Nishimura) at 0.26"; aperture means at 120, 250 pc
xs=[50,120,250]; ours={'N':[.99,.84,.69],'S':[.48,.40,.37]}
fig,(a,b)=plt.subplots(1,2,figsize=(12,4.8),gridspec_kw=dict(width_ratios=[1.3,1]))
cl={'lo':'#999999','mid':'#56B4E9','hi':'#CC79A7'}
for i in np.where(ok)[0]:
    g=None if np.isnan(A[i]) else ('hi' if A[i]>=.25 else 'mid' if A[i]>=.10 else 'lo')
    a.scatter(P[i],R[i],s=34,marker='D' if g=='hi' else 'o',facecolor=cl.get(g,'none'),edgecolor='k' if g else '0.4',lw=.6)
for k,c in (('N','#D55E00'),('S','#0072B2')):
    a.plot(xs,ours[k],'-',color=c,lw=1.5); a.scatter(xs,ours[k],s=80,marker='*',color=c,edgecolor='k',zorder=5)
    a.annotate(f'NGC 3256 {k}',(xs[0],ours[k][0]),xytext=(-6,-14),textcoords='offset points',color=c,ha='right',fontsize=9,weight='bold')
a.set_xscale('log'); a.set_xlim(30,700); a.set_xticks([50,100,200,400]); a.set_xticklabels(['50','100','200','400']); a.minorticks_off(); a.set_ylim(0,3); a.set_xlabel('spatial scale (pc)'); a.set_ylabel('HCN/HCO$^+$ J=3-2 ratio')
from matplotlib.lines import Line2D
h=[Line2D([],[],marker='o',ls='',mfc=cl['lo'],mec='k',label=r'Nishimura+24, AGN fraction <0.10'),Line2D([],[],marker='o',ls='',mfc=cl['mid'],mec='k',label='0.10-0.25'),
   Line2D([],[],marker='D',ls='',mfc=cl['hi'],mec='k',label=r'$\geq$0.25 (one nucleus)'),Line2D([],[],marker='o',ls='',mfc='none',mec='0.4',label='AGN fraction not listed'),
   Line2D([],[],marker='*',ls='-',color='k',ms=10,label='NGC 3256 (this work)')]
a.legend(handles=h,fontsize=8,loc='upper left',frameon=False)
a.set_title('(a) J=3-2, indicative only: definitions differ (see text)',fontsize=9)
# (b) J=4-3 vs Izumi+16 ranges (hi-res SB 0.24-1.03, AGN 1.11-3.19)
b.axvspan(0.24,1.03,color='#009E73',alpha=.25,lw=0); b.axvspan(1.11,3.19,color='#D55E00',alpha=.2,lw=0)
b.text(.63,4.6,'starburst\n0.24-1.03',ha='center',fontsize=9); b.text(2.15,4.6,'AGN\n1.11-3.19',ha='center',fontsize=9)
b.errorbar([.515],[3.2],fmt='*',ms=13,color='#D55E00',mec='k'); b.text(.58,3.2,'N, 4.5": 0.515',va='center',fontsize=9)
b.errorbar([.317],[2.4],xerr=[[.317-.22],[.0]],fmt='*',ms=13,color='#0072B2',mec='k',capsize=3,ecolor='#0072B2'); b.text(.12,2.0,'S, 4.5": 0.317; 0.22-0.32\nrange from the 1.0" high-resolution\ncubes, which lose 39% of S HCN flux',va='top',fontsize=8)
b.errorbar([.39],[0.8],xerr=[.11],fmt='o',color='0.3'); b.text(.55,0.8,'Izumi+16, NGC 3256 single dish: 0.39$\\pm$0.11',va='center',fontsize=8)
b.set_yticks([]); b.set_ylim(0,5.4); b.set_xlim(0,3.5); b.set_xlabel('HCN/HCO$^+$ J=4-3 brightness-temperature ratio')
b.set_title('(b) J=4-3 against the Izumi+16 ranges',fontsize=9)
fig.tight_layout(); fig.savefig('../figs/fig5_literature.png',dpi=160)
