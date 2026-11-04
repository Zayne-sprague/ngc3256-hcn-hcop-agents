import json,re,csv,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D=json.load(open('f4_profiles.json'))
# --- Bayes factors Habun/Hexc from s356 outputs (3 error models x 2 X priors x 8 physical priors = 48 cells)
ERR={'s0.025_m0.05':'errors halved','s0.05_m0.1':'default errors','s0.1_m0.2':'errors doubled'}
BF=[]
for tag,lab in ERR.items():
    for l in open(f'/workspace/out/q57/s356_{tag}.txt'):
        m=re.match(r'## (.*): X (\d+-\d+);',l)
        if m: cur=(m.group(1),m.group(2))
        m=re.search(r'BF\(Habun/Hexc\)=([\d.e+]+)',l)
        if m: BF.append((lab,cur[1],cur[0],float(m.group(1))))
print(len(BF),min(b[3] for b in BF),max(b[3] for b in BF))
fig=plt.figure(figsize=(7.5,9.5)); gs=fig.add_gridspec(3,1,height_ratios=[3,0.9,2.8],hspace=0.38)
ax=fig.add_subplot(gs[0]); st=fig.add_subplot(gs[1],sharex=ax); bx=fig.add_subplot(gs[2])
col={'Hexc':'#D55E00','Hcol':'#0072B2'}
sty={'X40-100 unconstrained':'-','X40-100 log n>=3.8':':','X20-200 unconstrained':'--'}
rows=[]
for k,v in D.items():
    b=10**np.array(v['bins'])
    for m in ('Hexc','Hcol'):
        y=np.array(v[m]); y0=y-y.min()
        ax.plot(b,y0,sty[k],color=col[m],lw=1.8 if k.startswith('X40-100 unc') else 1.1,label=None)
        for c,yy,yy0 in zip(b,y,y0): rows.append([k,m,round(c,3),round(yy,2),round(yy0,2)])
for lv in (1,4): ax.axhline(lv,color='0.4',lw=.8,ls='-.'); ax.text(1.74,lv+0.25,f'$\\Delta\\chi^2$ = {lv}',fontsize=8,ha='right',color='0.3')
ax.axvspan(1.09,1.13,color='#CC79A7',alpha=.35); ax.text(1.11,22.5,'best-fit\nequal-A\n1.09-1.13',fontsize=7.5,ha='center',va='top',color='#8a2d63')
ax.axvspan(1.448-0.015,1.448+0.015,color='k',alpha=.6); ax.text(1.46,22.5,'observed 1-0\n1.448$\\pm$0.015\n(stat; 2.2")',fontsize=7.5,ha='left',va='top')
ax.set_xlim(0.8,1.75); ax.set_ylim(0,23); ax.set_ylabel('$\\Delta\\chi^2$ = $\\Delta(-2\\ln L)$ above own best fit')
ax.set_title('Equal [HCN]/[HCO$^+$] at N and S: profile likelihood of the model 1-0 N/S contrast',fontsize=9.5)
h=[Line2D([],[],color=col['Hexc'],label='H$_{exc}$: excitation only (Tk, n, filling free)'),Line2D([],[],color=col['Hcol'],label='H$_{col}$: + N(HCO$^+$) free per nucleus'),
   Line2D([],[],color='k',lw=1.8,label='12C/13C 40-100, any n'),Line2D([],[],color='k',ls=':',label='12C/13C 40-100, log n>=3.8'),Line2D([],[],color='k',ls='--',label='12C/13C 20-200, any n')]
ax.legend(handles=h,fontsize=7.5,loc='upper center',frameon=False,bbox_to_anchor=(0.43,0.80))
plt.setp(ax.get_xticklabels(),visible=False)
st.axvspan(1.23,1.35,color='#CC79A7',alpha=.35); st.errorbar([1.503],[0],xerr=[0.017],fmt='ko',ms=5,capsize=2)
st.text(1.29,0.38,'3-2 best-fit equal-A\n1.23-1.35',fontsize=7.5,ha='center',va='bottom',color='#8a2d63'); st.text(1.47,0.0,'observed 3-2\n1.503$\\pm$0.017',fontsize=7.5,va='center',ha='right')
st.set_ylim(-0.6,1.2); st.set_yticks([]); st.set_xlabel('N/S contrast of HCN/HCO$^+$ (model 1-0 in the upper panel; 3-2 best-fit values only, no 3-2 profile computed)',fontsize=8)
# BF panel
pri=['unconstrained','log n>=3.8','log n>=4.2','log n 3-6 flat (alt prior)','log n>=3.8, L<=450pc, x<=1e-8','log n>=3.8, L<=100pc (N core), x<=1e-8','log n>=4.2, L<=100pc, x<=1e-8','log n>=3.8, L<=100pc, x<=3e-9 (tight)']
short=['any n','n>=6e3','n>=1.6e4','flat\nlog n 3-6','n>=6e3\nL<=450pc','n>=6e3\nL<=100pc','n>=1.6e4\nL<=100pc','n>=6e3, L<=100\nx<=3e-9']
ec={'errors halved':'#009E73','default errors':'#E69F00','errors doubled':'#CC79A7'}
for lab,xr,p,v in BF:
    i=pri.index(p); dx={'errors halved':-0.2,'default errors':0,'errors doubled':0.2}[lab]+(0.06 if xr=='20-200' else -0.06)
    bx.plot(i+dx,v,'o' if xr=='40-100' else '^',color=ec[lab],ms=6.5,mfc=ec[lab] if xr=='40-100' else 'white',mew=1.3)
bx.set_yscale('log'); bx.axhline(1,color='k',lw=.8); bx.axhline(2.58,color='0.4',ls='--',lw=.8); bx.text(7.6,2.75,'min 2.6',fontsize=8,ha='right',color='0.3')
bx.set_ylim(0.8,4000); bx.set_xticks(range(8)); bx.set_xticklabels(short,fontsize=6.5); bx.set_ylabel('Bayes factor  H$_{abun}$ / H$_{exc}$'); bx.minorticks_off()
bx.set_title('Same gas state, only [HCN]/[HCO$^+$] differs (H$_{abun}$) vs excitation only: 48 error-model x prior cells',fontsize=9)
h=[Line2D([],[],marker='o',ls='',color=c,label=l) for l,c in ec.items()]+[Line2D([],[],marker='o',ls='',color='k',label='12C/13C 40-100'),Line2D([],[],marker='^',ls='',color='k',mfc='white',label='12C/13C 20-200')]
bx.legend(handles=h,fontsize=7.5,ncol=2,loc='upper right',frameon=False)
fig.savefig('f4_radex_abundance.png',dpi=180,bbox_inches='tight'); fig.savefig('f4_radex_abundance.pdf',bbox_inches='tight')
with open('f4_profile_values.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['prior','hypothesis','model_1-0_contrast','delta_m2lnL_vs_Hfull_max','delta_vs_own_min']); w.writerows(rows)
with open('f4_bayes_factors.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['error_model','X_range_12C13C','physical_prior','BF_Habun_over_Hexc']); w.writerows(BF)
