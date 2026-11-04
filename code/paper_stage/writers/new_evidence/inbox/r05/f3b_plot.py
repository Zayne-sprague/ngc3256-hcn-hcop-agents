"""F3b: N/S line ratios, 1.1 mm continuum normalisation, Td band. Inputs: /workspace/out/q52/s320.json (line integrals K km/s +- stat, +-350 km/s),
s322.txt (1.1 mm continuum N/S at 1.4", 2.5"; Xb4d spw25_27_29_31 cont, single pixel at fixed N/S positions)."""
import json,csv,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=json.load(open('/workspace/out/q52/s320.json'))
def I(J,th,n,l):
    if l=='CN' and J=='10':
        a,b=R[f'{J}|{th}|{n}|CN_G1']['350'],R[f'{J}|{th}|{n}|CN_G2']['350']; return a[0]+b[0],np.hypot(a[1],b[1])
    if l=='CN': l='CN_G1'; a=R[f'{J}|{th}|{n}|CN_G1']['350']; b=R[f'{J}|{th}|{n}|CN_G2']['350']
    else: a=R[f'{J}|{th}|{n}|{l}']['350']; b=(0,0)
    return a[0],a[1]
# 3-2 CN: single 'CN' entry (whole 3-2 complex, one window) in s320.json
def I32(th,n,l):
    k=l; a=R[f'32|{th}|{n}|{k}']['350']; return a[0],a[1]
def NS(J,th,l):
    f=I if J=='10' else (lambda J_,t,n,l_: I32(t,n,l_))
    N,eN=f(J,th,'N',l); S,eS=f(J,th,'S',l); v=N/S; return v,v*np.hypot(eN/N,eS/S)
CONT={('32',1.4):1.25,('10',2.5):1.54}   # s322.txt: N/S of 1.1 mm continuum at 1.4" and 2.5" (N 64.4/S 51.3 mK; 33.6/21.8 mK)
sp=['HCO+','HCN','HNC','CN']
rows=[]
for J,th in [('10',2.5),('32',1.4)]:
    for l in sp:
        v,e=NS(J,th,l); c=CONT[(J,th)]
        rows.append(dict(J=J,beam=th,line=l,NS_line=v,err=e,cont_NS=c,NS_norm=v/c,err_norm=e/c))
# X/HCO+ ratio-of-ratios
rr={}
for J,th in [('10',2.5),('32',1.4)]:
    for l in sp[1:]:
        v,e=NS(J,th,l); h,eh=NS(J,th,'HCO+'); rr[(J,l)]=(v/h,v/h*np.hypot(e/v,eh/h))
cJ={'10':'#009E73','32':'#0072B2'}
fig,(ax,bx)=plt.subplots(1,2,figsize=(10,5),gridspec_kw={'width_ratios':[3,1.4],'wspace':0.18})
slot={('10','line'):-0.3,('10','norm'):-0.1,('32','line'):0.1,('32','norm'):0.3}
for r in rows:
    i=sp.index(r['line'])
    for kind,val,err in [('line',r['NS_line'],r['err']),('norm',r['NS_norm'],r['err_norm'])]:
        x=i+slot[(r['J'],kind)]
        if kind=='norm':
            ax.plot([x,x],[val*0.5,val*2],color=cJ[r['J']],lw=7,alpha=0.22,solid_capstyle='butt',zorder=1)
        ax.errorbar(x,val,err,fmt='s' if r['J']=='10' else 'o',color=cJ[r['J']],mfc=cJ[r['J']] if kind=='line' else 'white',ms=7,mew=1.5,capsize=2,zorder=3)
        ax.text(x,val*(1.0 if kind=='line' else 1.0),'',fontsize=6)
ax.axhline(1,color='k',lw=.8)
for (J,th),c in CONT.items():
    ax.axhline(c,color=cJ[J],ls='--',lw=1)
    ax.text(3.68,c,f'cont. {c:.2f}\n({th}")',color=cJ[J],fontsize=7.5,ha='left',va='center',clip_on=False)
ax.set_yscale('log'); ax.set_ylim(0.4,6.5); ax.set_xlim(-0.6,3.65)
ax.set_yticks([0.5,1,2,3,5]); ax.set_yticklabels(['0.5','1','2','3','5']); ax.minorticks_off()
ax.set_xticks(range(4)); ax.set_xticklabels(['HCO$^+$','HCN','HNC','CN'],fontsize=11)
ax.set_ylabel('N/S ratio')
ax.set_title('Line N/S and line/continuum N/S',fontsize=10)
h=[Line2D([],[],marker='s',ls='',color=cJ['10'],label='1-0, 2.5": line N/S'),Line2D([],[],marker='s',ls='',color=cJ['10'],mfc='white',label='1-0, 2.5": line/1.1 mm cont. N/S'),
   Line2D([],[],marker='o',ls='',color=cJ['32'],label='3-2, 1.4": line N/S'),Line2D([],[],marker='o',ls='',color=cJ['32'],mfc='white',label='3-2, 1.4": line/1.1 mm cont. N/S'),
   plt.Rectangle((0,0),1,1,color='0.6',alpha=.4,label='shade: Td$_N$/Td$_S$ = 0.5-2 (RJ: x0.5 to x2)')]
ax.legend(handles=h,fontsize=7.5,loc='upper left',frameon=False)
# panel B
for j,(J,th) in enumerate([('10',2.5),('32',1.4)]):
    for i,l in enumerate(sp[1:]):
        v,e=rr[(J,l)]; bx.errorbar(i+(-0.1 if J=='10' else 0.1),v,e,fmt='s' if J=='10' else 'o',color=cJ[J],ms=7,capsize=2)
bx.axhline(1,color='k',lw=.8); bx.set_yscale('log'); bx.set_ylim(0.8,3.2); bx.set_yticks([1,1.5,2,3]); bx.set_yticklabels(['1','1.5','2','3']); bx.minorticks_off()
bx.set_xticks(range(3)); bx.set_xticklabels(['HCN','HNC','CN']); bx.set_xlim(-0.5,2.5)
bx.set_ylabel('N/S of (line/HCO$^+$) = N/S of (line/cont.) relative to HCO$^+$'); bx.yaxis.label.set_size(7.5)
bx.set_title('Relative to HCO$^+$\n(continuum and Td cancel)',fontsize=9)
fig.savefig('f3b_ns_species.png',dpi=180,bbox_inches='tight'); fig.savefig('f3b_ns_species.pdf',bbox_inches='tight')
with open('f3b_values.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['J','beam_arcsec','species','N_over_S_line','stat_err','cont_1.1mm_N_over_S','line_over_cont_N_over_S','stat_err_norm','norm_Td_ratio_0.5','norm_Td_ratio_2','X_over_HCO+_NS','err_X_over_HCO+_NS'])
    for r in rows:
        x=rr.get((r['J'],r['line']),('',''))
        w.writerow([{'10':'1-0','32':'3-2'}[r['J']],r['beam'],r['line'],round(r['NS_line'],3),round(r['err'],3),r['cont_NS'],round(r['NS_norm'],3),round(r['err_norm'],3),round(r['NS_norm']*.5,3),round(r['NS_norm']*2,3)]+[round(t,3) if t!='' else '' for t in x])
for r in rows: print(r['J'],r['line'],round(r['NS_line'],2),round(r['err'],2),round(r['NS_norm'],2))
print(rr)
