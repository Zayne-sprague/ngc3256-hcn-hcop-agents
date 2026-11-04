"""Figure 3 (R2): N/S contrast of HCN/HCO+ vs J, and N/S per unit 1.3 mm continuum. Numbers transcribed from
published team snapshots (paths in fig3_numbers.csv / caption.md); run: python3 fig3.py"""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, csv
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
Cb,Co,Cg,Cp='#0072B2','#E69F00','#009E73','#CC79A7'
def spread(vals,stat):
    return np.mean(vals), float(np.hypot(stat,(max(vals)-min(vals))/2))
rows=[]  # (panel,label,value,err,src)
# (a) 2.2": p0036 s192 (comb/c12 x w350/h150), 1-0 from 2015.1.00993
v32=[1.441,1.503,1.451,1.518]; m,e=spread(v32,0.023); a22_32=(m,e)
v10=[1.467,1.448]; a22_10=spread(v10,0.022)
v32=[1.532,1.634,1.578,1.666]; a15_32=spread(v32,0.026)
a15_10=spread([1.616,1.616],0.059)
# 2.5" p0060 s345 (pixel/aperture estimators, 1-0 X961, 3-2 Xb4d 12m+7m)
a25_10=(1.49,0.07); a25_32=(1.64,0.13)   # means/half-range over 4 estimators: 1-0 1.40-1.53, 3-2 1.48-1.74
fig=plt.figure(figsize=(12,5.2)); gs=fig.add_gridspec(1,2,width_ratios=[1.35,1],wspace=0.28)
ax=fig.add_subplot(gs[0]); 
groups=[('2.2" (≈400 pc)',a22_10,a22_32),('1.5"',a15_10,a15_32),('2.5" (estimator spread)',a25_10,a25_32)]
x=0; ticks=[];tl=[]
for g,(lab,u,v) in enumerate(groups):
    ax.bar(x,u[0],0.8,yerr=u[1],color=Cb,capsize=3); ax.bar(x+1,v[0],0.8,yerr=v[1],color=Co,capsize=3)
    for xx,q in ((x,u),(x+1,v)): ax.text(xx,q[0]+q[1]+0.06,f'{q[0]:.2f}',ha='center',fontsize=8)
    ticks+=[x,x+1]; tl+=['1-0','3-2']; ax.text(x+0.5,-0.42,lab.replace(' (estimator spread)',' (est.)'),ha='center',va='top',fontsize=8); x+=3
# ~1" flux-loss-biased
ax.bar(x,1.45,0.8,yerr=0.22,fill=False,hatch='//',edgecolor=Cb,capsize=3); ax.bar(x+1,2.08,0.8,fill=False,hatch='//',edgecolor=Co)
ax.text(x,1.45+0.22+0.06,'1.45',ha='center',fontsize=8); ax.text(x+1,2.14,'2.08',ha='center',fontsize=8)
ticks+=[x,x+1]; tl+=['1-0','3-2']; ax.text(x+0.5,-0.42,'1.0-1.1"\n(flux-loss\nbiased)',ha='center',va='top',fontsize=8); x+=3
# 4-3
ax.bar(x,1.62,0.8,fill=False,hatch='xx',edgecolor=Cg); ax.bar(x+1,3.67,0.8,yerr=0.26,fill=False,hatch='xx',edgecolor=Cg,capsize=3)
ax.text(x,1.68,'1.6',ha='center',fontsize=8); ax.text(x+1,3.67+0.26+0.06,'3.7',ha='center',fontsize=8)
ax.annotate('',xy=(x+1.9,3.67),xytext=(x+1.9,1.62),arrowprops=dict(arrowstyle='<->',color='k'))
ax.text(x+1.0,4.35,'4-3 unresolved\n(4.5" Xb47 vs 1.0" 00493+00223)',ha='center',fontsize=8)
ticks+=[x,x+1]; tl+=['4.5"','1.0"']; ax.text(x+0.5,-0.42,'4-3',ha='center',va='top',fontsize=8)
ax.set_xticks(ticks); ax.set_xticklabels(tl,fontsize=8); ax.set_ylim(0,4.9); ax.axhline(1,c='gray',lw=.6,ls=':')
ax.set_ylabel('N/S contrast of HCN/HCO+',fontsize=9); ax.set_title('(a) N/S contrast by J, matched beams',loc='left',fontsize=10)
# D inset
# separate D strip as inset axes top-left
d=ax.inset_axes([0.05,0.60,0.30,0.30]); d.set_xlim(0.6,1.25); d.set_ylim(-0.6,3.6)
d.errorbar([1.00],[3],xerr=0.05,fmt='o',c=Cb,capsize=2); d.errorbar([1.07],[2],xerr=0.05,fmt='s',c=Co,capsize=2)
d.errorbar([1.10],[1],xerr=0.06,fmt='o',mfc='w',c='gray',capsize=2); d.errorbar([0.71],[0],fmt='s',mfc='w',c='gray')
d.yaxis.tick_right(); d.spines['right'].set_visible(True); d.spines['left'].set_visible(False); d.set_yticks([3,2,1,0]); d.set_yticklabels(['p0036 1.00±0.05','p0060 1.07±0.05','p0034 1.10±0.06 (sup.)','0.71 (HCN R32/10 S/N, sup.)'],fontsize=6.5)
d.axvline(1,c='gray',lw=.6,ls=':'); d.set_xlabel('D = contrast(3-2)/contrast(1-0)',fontsize=7); d.tick_params(labelsize=7)
d.patch.set_alpha(0.9)
# (b)
bx=fig.add_subplot(gs[1]); cont={'1-0 2.5"':1.541,'3-2 1.1"':1.0746,'3-2 1.4"':1.2551}
data={ # N/S line ratio: 1-0 2.5" (team-b s322, 2015.1.00993), 3-2 (team-a norm_cont.csv)
 'HCN':[2.21,3.307,3.212],'HNC':[2.69,3.234,2.969],'CN':[2.98,3.419,3.740],'HCO+':[1.42,1.548,1.581]}
cols={'HCN':Cb,'HNC':Co,'CN':Cg,'HCO+':Cp}
labs=['1-0 (2.5")','3-2 (1.1")','3-2 (1.4")']; w=0.2; out=[]
for i,(sp,vals) in enumerate(data.items()):
    for j,k in enumerate(cont):
        y=vals[j]/cont[k]; err=0.10*y
        xx=j+(i-1.5)*w
        bx.bar(xx,y,w*0.92,color=cols[sp],label=sp if j==0 else None)
        bx.errorbar(xx,y,yerr=err,c='k',capsize=2,lw=.8)
        bx.plot([xx,xx],[y*0.5,y*2],c='k',lw=3,alpha=.25,solid_capstyle='butt')
        out.append((sp,k,vals[j],cont[k],y,err))
bx.axhline(1,c='gray',lw=.6,ls=':'); bx.set_xticks(range(3)); bx.set_xticklabels(labs)
bx.set_ylabel('N/S of line ÷ 1.3 mm continuum'); bx.set_title('(b) N/S per unit continuum',loc='left',fontsize=10)
bx.legend(frameon=False,ncol=4,fontsize=8,loc='upper center',bbox_to_anchor=(0.5,-0.08)); bx.set_ylim(0,6.5)
bx.text(0.02,0.82,'grey band: ×0.5-2 for Td$_N$/Td$_S$\n(per-dust-mass value = bar × Td$_N$/Td$_S$);\nerror bars 10% (stat+sampling)',transform=bx.transAxes,fontsize=7.5)
fig.subplots_adjust(bottom=0.2); fig.savefig('contrast_NS.png',dpi=180,bbox_inches='tight')
with open('fig3_numbers.csv','w',newline='') as f:
    w_=csv.writer(f); w_.writerow(['panel','beam','J','value','err','source'])
    for lab,u,v in groups:
        w_.writerow(['a',lab,'1-0',round(u[0],3),round(u[1],3),'p0036 s192 / p0060 s345']); w_.writerow(['a',lab,'3-2',round(v[0],3),round(v[1],3),''])
    w_.writerow(['a','1.0"','1-0',1.45,0.22,'team-a q55f r10_core.csv (00965, flux-loss)']); w_.writerow(['a','1.1"','3-2',2.08,'', 'team-a Q55a Xb4d 12m only'])
    w_.writerow(['a','4.5" Xb47','4-3',1.62,'','0.515/0.317']); w_.writerow(['a','1.0" 00493+00223','4-3',3.67,0.26,'0.814/0.222'])
    w_.writerow(['b:species','beam','N/S line','N/S cont','N/S per cont','err(10%)'])
    for r in out: w_.writerow(['b']+[r[0],r[1]]+[round(q,3) for q in r[2:]])
for a,b in zip(['a22_10','a22_32','a15_10','a15_32'],[a22_10,a22_32,a15_10,a15_32]): print(a,b)
for r in out: print(r)
