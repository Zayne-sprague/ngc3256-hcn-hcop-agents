"""F5 (R1 artefact check): Xb4b 0.30" spectra of HCN/HCO+ 3-2 at N (r<0.15") and S, baseline windows/variants, N core ratio per variant.
Reuses q64r2c_images.py machinery (same crop cubes, FFT convolution to a 0.30" circular beam, per-pixel polynomial baseline, +-350 km/s moment).
python3 -I f5_baseline_check.py  -> f5_baseline_check.png, f5_spectra.csv, f5_variants.csv"""
import sys, csv, numpy as np
sys.path.insert(0,'/workspace/scripts')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import q64r2c_images as Q
O='/workspace/paper_answers/r04/'; VSYS=2780.; VC=Q.VC; TH=0.30; HW=350; RAP=0.15
LINES={'HCO+':'25','HCN':'27'}
WIN={'A (adopted)':([(-850,-450),(450,600)],1),'order 0':([(-850,-450),(450,600)],0),'C (+far blue)':([(-1400,-1150),(-850,-450),(450,600)],1),
     'blue-only (-850..-450)':([(-850,-450)],0),'red-only':([(450,600)],0),'pipeline only':(None,1)}
cols={'A (adopted)':'#000000','order 0':'#56B4E9','C (+far blue)':'#009E73','blue-only (-850..-450)':'#D55E00','red-only':'#CC79A7','pipeline only':'#999999'}
# 0.30" beam cubes, aperture means
cube={};
for ln,s in LINES.items():
    c,v,info=Q.load('Xb4b',s); c=np.nan_to_num(c); X,Y=Q.xy(c.shape[1:],info)
    nat=Q.kernel_ft(c.shape[1:],info['pix'],info['bmaj'],info['bmin'],info['bpa']); tgt=Q.kernel_ft(c.shape[1:],info['pix'],TH,TH,0)
    F=np.fft.fft2(c,axes=(1,2))*tgt/np.maximum(nat,1e-6); F[:,nat<1e-4]=0
    cc=np.fft.ifft2(F,axes=(1,2)).real*(TH*TH)/(info['bmaj']*info['bmin'])
    cube[ln]=(cc,v,X,Y)
def aperture(ln,key):
    cc,v,X,Y=cube[ln]; px,py=Q.POS[key]; m=np.hypot(X-px,Y-py)<RAP
    return cc[:,m].mean(1),v
def baseline(spec,v,fw,od):
    if fw is None: return np.zeros_like(spec)
    m=np.zeros(len(v),bool)
    for a,b in fw: m|=(v>=a)&(v<=b)
    p=np.polyfit(v[m]/1000,spec[m],od); return np.polyval(p,v/1000)
rows=[];res={}
S={(ln,k):aperture(ln,k) for ln in LINES for k in 'NS'}
for k in 'NS':
    for wn,(fw,od) in WIN.items():
        I={}
        for ln in LINES:
            sp,v=S[(ln,k)]; b=baseline(sp,v,fw,od); dv=abs(np.median(np.diff(v))); sel=(v>=-HW)&(v<=HW)
            I[ln]=((sp-b)[sel].sum()*dv, b[np.argmin(abs(v))]*1e3)
        r=I['HCN'][0]/I['HCO+'][0]; res[(k,wn)]=r
        rows.append(dict(core=k,variant=wn,HCN_Jybm_kms=I['HCN'][0],HCOp_Jybm_kms=I['HCO+'][0],ratio_HCN_HCOp=r,base_at_line_HCN_mJy=I['HCN'][1],base_at_line_HCOp_mJy=I['HCO+'][1]))
with open(O+'f5_variants.csv','w') as f:
    w=csv.DictWriter(f,list(rows[0])); w.writeheader(); w.writerows(rows)
# spectra csv (pipeline cont-sub, 0.30" beam) + adopted baseline
with open(O+'f5_spectra.csv','w') as f:
    w=csv.writer(f); w.writerow(['v_minus_2780_kms']+[f'{ln}_{k}_mJybm' for ln in LINES for k in 'NS']+[f'{ln}_{k}_baselineA_mJybm' for ln in LINES for k in 'NS'])
    v=S[('HCN','N')][1]+VC-VSYS
    cols_=[S[(ln,k)][0]*1e3 for ln in LINES for k in 'NS']+[baseline(S[(ln,k)][0],S[(ln,k)][1],*WIN['A (adopted)'])*1e3 for ln in LINES for k in 'NS']
    for i in range(len(v)): w.writerow([f'{v[i]:.2f}']+[f'{c[i]:.4f}' for c in cols_])
# rms of line-free channels
for ln in LINES:
    for k in 'NS':
        sp,v=S[(ln,k)]; m=((v>=-850)&(v<=-450))|((v>=450)&(v<=600)); print(ln,k,'rms in baseline windows after A fit (mJy/b): %.3f'%np.std((sp-baseline(sp,v,*WIN['A (adopted)']))[m]*1e3),
          ' A baseline at line centre %.3f ; blue-only %.3f'%(baseline(sp,v,*WIN['A (adopted)'])[np.argmin(abs(v))]*1e3,baseline(sp,v,*WIN['blue-only (-850..-450)'])[np.argmin(abs(v))]*1e3))
for k,v_ in res.items(): print(k,'%.3f'%v_)
# figure
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(11,8.2)); gs=fig.add_gridspec(2,4,height_ratios=[1.55,1],hspace=0.38,wspace=0.35)
x0=VC-VSYS
def shade(ax,ln):
    ax.axvspan(x0-350,x0+350,color='#f0e442',alpha=.25,lw=0)
    for a,b in [(-850,-450),(450,600)]: ax.axvspan(x0+a,x0+b,color='#0072B2',alpha=.18,lw=0)
    ax.axvspan(x0-1400,x0-1150,color='#009E73',alpha=.12,lw=0)
    ax.axvspan(x0-1300,x0-700,fc='none',ec='#D55E00',hatch='////',lw=0,alpha=.5)
    if ln=='HCO+': ax.axvspan(x0-1100,x0-900,color='#D55E00',alpha=.22,lw=0)
for j,(ln,k) in enumerate([('HCN','N'),('HCO+','N'),('HCN','S'),('HCO+','S')]):
    ax=fig.add_subplot(gs[0,j]); sp,v=S[(ln,k)]; vv=v+VC-VSYS
    shade(ax,ln); ax.step(vv,sp*1e3,color='0.35',lw=.7,where='mid')
    for wn in ('A (adopted)','blue-only (-850..-450)','C (+far blue)'):
        fw,od=WIN[wn]; ax.plot(vv,baseline(sp,v,fw,od)*1e3,color=cols[wn],lw=1.4,ls={'A (adopted)':'-','blue-only (-850..-450)':'--','C (+far blue)':':'}[wn],label=wn)
    ax.set_xlim(-1500,1500); ax.set_title(f'{ln} 3-2, {k} core',fontsize=10)
    top=np.nanmax(sp)*1e3; ax.set_ylim(-0.12*top,top*1.12)
    if ln=='HCO+': ax.annotate('weak extra\nline',xy=(x0-1000,np.nanmax(sp[(v<-900)&(v>-1100)])*1e3),xytext=(-1450,top*.45),fontsize=8,color='#D55E00',arrowprops=dict(arrowstyle='->',color='#D55E00'))
    ax.set_xlabel('v - 2780 km/s'); ax.set_ylabel('mJy/beam' if j in(0,) else '')
    if j==0: ax.legend(fontsize=7,loc='upper left',frameon=False)
# inset zoom of baselines is in caption numbers; add a legend patch row
import matplotlib.patches as mp
fig.legend(handles=[mp.Patch(color='#f0e442',alpha=.4,label='line window +-350 (about 2753)'),mp.Patch(color='#0072B2',alpha=.3,label='adopted baseline windows'),mp.Patch(color='#009E73',alpha=.25,label='extra window, variant C'),mp.Patch(color='#D55E00',alpha=.35,label='extra HCO+ line (-1100..-900)'),mp.Patch(fc='none',ec='#D55E00',hatch='////',label='uv variant B window (-1300..-700), rejected')],loc='upper center',ncol=3,fontsize=8,frameon=False,bbox_to_anchor=(.5,.975))
# beam ellipse
axb=fig.add_axes([0.905,0.80,0.07,0.07]); axb.add_patch(Ellipse((.5,.5),.8,.8,fc='0.4')); axb.set_xlim(0,1);axb.set_ylim(0,1);axb.axis('off'); axb.set_title('0.30" beam',fontsize=7)
# variant panel
ax=fig.add_subplot(gs[1,:2]); names=list(WIN); y=np.arange(len(names))[::-1]
for i,wn in enumerate(names): ax.barh(y[i],res[('N',wn)],color=cols[wn],height=.6); ax.text(res[('N',wn)]+.01,y[i],'%.2f'%res[('N',wn)],va='center',fontsize=9)
ax.set_yticks(y); ax.set_yticklabels(names); ax.set_xlim(0,1.45); ax.set_xlabel('N core HCN/HCO+ 3-2 (0.30" beam, r<0.15", +-350 km/s)')
ax.axvline(res[('N','A (adopted)')],color='k',lw=.6,ls=':')
ax=fig.add_subplot(gs[1,2:]); ax.axis('off')
ax.text(0,.95,'uv-plane check (Xb4d+ACA, 0.4" model):\nvariant A (blue -1300..-700 + red >+560): reduced chi2 5.5, N = 0.96\nvariant B (blue -1300..-700 only): reduced chi2 47.7, N = 1.12\n=> B rejected. Its window contains the extra HCO+ line.\n\nImage plane (left, Xb4b): windows -850..-450 / +450..+600\nexclude that line, so the variants move N by <=0.05.',va='top',fontsize=9)
fig.suptitle('F5: continuum-baseline artefact check, Xb4b 0.30" beam',y=1.02,fontsize=11)
fig.savefig(O+'f5_baseline_check.png',dpi=170,bbox_inches='tight')
