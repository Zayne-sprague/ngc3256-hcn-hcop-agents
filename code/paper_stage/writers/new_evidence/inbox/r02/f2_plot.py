"""F2: N/S HCN/HCO+ contrast vs effective beam FWHM. All numbers are transcribed from the cited posts/outputs (paths in 'src');
4-3 12-m/ACA image-plane points and their stat errors are recomputed from /snapshots/team-a-s34/out/q60/imgplane.csv."""
import csv,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# columns: id, J, beam_arcsec, val, err_lo_stat, err_hi_stat, err_sys_lo, err_sys_hi, complete(1/0), plane(img/uv), src, estimator/flag
rows=[]
def add(*a): rows.append(a)
# ---- 3-2, 12-m only, Xb4b/Xb4d image plane
add('p0064_Xb4b_r<0.25','3-2',0.30,2.26,0.09,0.09,0.12,0.063,0,'img','p0064 out/q64r2/final_budget.txt','annulus sums r<0.25", med350; published +0.11/-0.15 total; stat 0.09, baseline 0.05, array 0.03, resid/restored 0/-0.08')
add('p0064_Xb4d_r<1','3-2',1.05,1.98,0.06,0.06,0.08,0.08,0,'img','p0064 out/q64r2/final_budget.txt','annulus sums r<1", 12-m Xb4d/Xb4b 1.05" beam; published 2.0+-0.1')
for b,v in [(0.26,2.25),(0.2835,2.32),(0.6186,2.11),(1.2887,1.9)]:
    add(f'p0031_Xb4b_ap_{b}','3-2',b,v,0.0,0.0,0.2,0.2,0,'img','p0031 out/q26/s171_summary.txt (post text)','FLAG: aperture r=beam FWHM at HCO+ peak, RJ T_b; error = adopted +-0.2 estimator spread (stat 0.01-0.03 on each ratio)')
for b,v,lo,hi in [(0.246,1.94,1.64,2.07),(0.3,2.13,1.92,2.25),(0.4,2.21,2.06,2.31),(0.6,2.26,2.19,2.29),(1.0,2.14,2.09,2.18),(1.5,2.03,1.97,2.08)]:
    add(f'p0065_teamA_Xb4b_{b}','3-2',b,v,0,0,v-lo,hi-v,0,'img','p0065 snapshots/team-a-s33/out/q64r2/summary.txt','FLAG team-a estimator: median of HCN/HCO+ over baseline windows (A,B,C...), N/S of window medians; bar = min-max of window variants, not stat')
add('p0036_c12_1.5','3-2',1.5,1.578,0.027,0.027,0.05,0.05,0,'img','p0036 snapshots/team-b-s21/out/q34/s192.txt (c12 w350)','12-m only (1 of 2 12-m EBs) uvtaper to 1.5"; sys 0.05 = 12m+7m vs 12m-only')
add('p0036_c12_2.2','3-2',2.2,1.451,0.027,0.027,0.05,0.05,0,'img','p0036 s192.txt (c12 w350)','12-m only uvtaper to 2.2"')
# ---- 3-2 12m+7m / ACA
add('p0036_comb_1.5','3-2',1.5,1.532,0.026,0.026,0.05,0.05,1,'img','p0036 s192.txt (comb w350)','12m+7m (1/2 12-m EB, 3/14 7-m EBs), uvtaper; image-plane beam matching, not strict uv matching')
add('p0036_comb_2.2','3-2',2.2,1.441,0.023,0.023,0.05,0.05,1,'img','p0036 s192.txt (comb w350)','12m+7m uvtaper to 2.2"')
add('p0031_Xb4f_6.5','3-2',6.5,1.07,0.0,0.0,0.10,0.10,1,'img','p0031 (N 0.46, S 0.43 at 6.5")','FLAG: ACA 7-m only; 6.5" circular beam; err <=0.10 propagated from quoted max stat 0.03 per ratio; N+S blended')
add('p0031_Xb4d_6.5','3-2',6.5,1.10,0.0,0.0,0.10,0.10,0,'img','p0031 (N 0.43, S 0.39 at 6.5")','12-m Xb4d convolved to 6.5"; same error as above')
# ---- 1-0 (short-spacing-complete archive cubes)
add('p0036_1-0_1.5','1-0',1.5,1.616,0.059,0.059,0,0,1,'img','p0036 s192.txt (X145 1.5")','2015.1.00412 X145 (Harada project) short-spacing-complete cube; stat only')
add('p0036_1-0_2.2','1-0',2.2,1.467,0.022,0.022,0,0,1,'img','p0036 s192.txt / p0034 (X961 2.2")','2015.1.00993 X961; stat only')
add('p0058_1-0_2.5','1-0',2.5,1.56,0.03,0.03,0,0,1,'img','p0058 out/q52 (HCN 2.21/HCO+ 1.42 N/S line ratios; 0.912/0.585)','FLAG: pixel at fixed N,S positions, linear baseline, +-350 km/s; X961 convolved to 2.5"; stat only')
# ---- 4-3
add('p0066_uvfit','4-3',5.0,1.8,0.12,0.12,0.2,0.3,1,'uv','p0066 snapshots/team-a-s34/out/q65/tables.txt','FLAG: uv-plane two-Gaussian fit, uv>8 kl, ACA-visible flux; NO defined beam, x=5" is nominal ACA resolution (bar 4-7"); also not same estimator')
D=list(csv.DictReader(open('/snapshots/team-a-s34/out/q60/imgplane.csv')))
def get(b,l,s,r):
    for d in D:
        if float(d['beam'])==b and d['line']==l and d['src']==s and d['reg']==r: return float(d['I']),float(d['e'])
for b in (4.5,5.5,7.0):
    for arr,sh,hs,comp in [('ACA','Xb47','Xb47',1),('12m','X1df1','X858',0)]:
        q=[get(b,'HCN',sh,'N'),get(b,'HCN',sh,'S'),get(b,'HCO+',hs,'N'),get(b,'HCO+',hs,'S')]
        v=q[0][0]/q[1][0]/(q[2][0]/q[3][0]); e=v*np.sqrt(sum((x[1]/x[0])**2 for x in q))
        add(f'p0060_4-3_{arr}_{b}','4-3',b,round(v,3),round(e,3),round(e,3),0,0,comp,'img','p0060/p0066 snapshots/team-a-s34/out/q60/imgplane.csv',
            'FLAG: pixel values; '+('ACA Xb47 both lines' if comp else '12-m: HCN 2018.1.00493 X1df1, HCO+ 2018.1.00223 X858 (different projects/EBs); not re-imaged')+'; stat only from file errors')
OFF={'p0031_Xb4b_ap_0.26':0.93,'p0031_Xb4b_ap_0.2835':1.0,'p0064_Xb4b_r<0.25':1.07,'p0036_c12_1.5':0.95,'p0036_1-0_1.5':1.05,'p0036_c12_2.2':0.95,'p0036_1-0_2.2':1.05,'p0036_comb_2.2':1.0,'p0031_Xb4d_6.5':0.95}
cols={'1-0':'#009E73','3-2':'#0072B2','4-3':'#D55E00'}
fig,ax=plt.subplots(figsize=(8,5.2))
ax.axhspan(1.48-0.03,1.48+0.03,color='0.8',zorder=0,label='_nolegend_')
ax.axhline(1,color='k',lw=0.8)
for r in rows:
    id_,J,b,v,sl,sh,yl,yh,comp,pl,src,note=r
    lo=np.hypot(sl,yl); hi=np.hypot(sh,yh)
    team_a=id_.startswith('p0065'); m='s' if pl=='uv' else 'o'
    ms=4 if team_a else 7
    x=b*(1.04 if team_a else OFF.get(id_,1.0))
    ax.errorbar(x,v,[[lo],[hi]],fmt=m,color=cols[J],mfc=cols[J] if comp else 'white',mec=cols[J],ms=ms,mew=1.3 if not team_a else 0.8,
        elinewidth=1.0 if not team_a else 0.6,capsize=2,alpha=0.6 if team_a else 1,zorder=3)
    if r[0]=='p0066_uvfit': ax.plot([4,7],[v,v],color=cols[J],lw=1,alpha=.5)
    if 'FLAG' in note and not team_a: ax.annotate('†',(x,v),xytext=(5,5),textcoords='offset points',fontsize=10,color=cols[J])
ax.set_xscale('log'); ax.minorticks_off(); ax.set_xlim(0.2,9); ax.set_ylim(0.8,3.2)
ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8'])
ax.set_xlabel('effective beam FWHM  (arcsec)'); ax.set_ylabel('N/S contrast of HCN/HCO$^+$ ratio')
top=ax.secondary_xaxis('top',functions=(lambda r:r*194,lambda p:p/194)); top.set_xlabel('beam FWHM (pc), 1" = 194 pc'); top.minorticks_off()
top.set_xticks([50,100,200,400,800,1600]); top.set_xticklabels(['50','100','200','400','800','1600'])
h=[Line2D([],[],marker='o',ls='',color=cols[j],label=f'J={j}') for j in cols]
h+=[Line2D([],[],marker='o',ls='',color='k',label='flux-complete (filled)'),Line2D([],[],marker='o',ls='',color='k',mfc='white',label='12-m only (open)'),
    Line2D([],[],marker='o',ls='',color='k',label='image plane (circle)'),Line2D([],[],marker='s',ls='',color='k',label='uv plane (square)'),
    Line2D([],[],marker='o',ls='',color='0.4',ms=4,mfc='white',alpha=.6,label='small: team-a p0065 estimator'),
    Line2D([],[],marker='$†$',ls='',color='k',label='† beam/estimator differs'),
    plt.Rectangle((0,0),1,1,color='0.8',label='Harada+2018 1.7" 1-0: 1.48')]
ax.legend(handles=h,fontsize=7.5,ncol=2,loc='upper right',frameon=False)
fig.savefig('f2_contrast_vs_beam.png',dpi=180,bbox_inches='tight'); fig.savefig('f2_contrast_vs_beam.pdf',bbox_inches='tight')
with open('f2_values.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['id','J','beam_fwhm_arcsec','beam_pc','N_over_S','stat_lo','stat_hi','sys_lo','sys_hi','total_lo_plotted','total_hi_plotted','flux_complete','plane','source','estimator_note'])
    for r in rows:
        id_,J,b,v,sl,sh,yl,yh,comp,pl,src,note=r
        w.writerow([id_,J,b,round(b*194,1),v,sl,sh,yl,yh,round(np.hypot(sl,yl),3),round(np.hypot(sh,yh),3),comp,pl,src,note])
