import json,csv,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
J=json.load(open('f1_annuli.json'))
ANN=[(0,.15),(.15,.25),(.25,.5),(.5,1),(1,1.5),(1.5,2),(2,3)]
mid=np.array([(a+b)/2 for a,b in ANN]); lab=[f'{a:g}-{b:g}' for a,b in ANN]
def ratio(key,n,i=None):
    h=np.array(J[key.replace('LINE','HCN')][n]); o=np.array(J[key.replace('LINE','HCO+')][n]); return h/o
def prof(arr,t,var='med350'):
    return {n:ratio(f'{arr}|{var}|{t}|LINE',n) for n in 'NS'}
def stat(arr,t):
    r=prof(arr,t); e={}
    for n in 'NS':
        h=np.array(J[f'{arr}|med350|{t}|HCN'][n]);o=np.array(J[f'{arr}|med350|{t}|HCO+'][n])
        eh=np.array(J[f'{arr}|null|{t}|HCN'][n]);eo=np.array(J[f'{arr}|null|{t}|HCO+'][n])
        e[n]=r[n['0'] if False else n]*np.hypot(eh/h,eo/o)
    return e
def base(arr,t):
    c=prof(arr,t); e={n:np.zeros(len(ANN)) for n in 'NS'}; e['NS']=np.zeros(len(ANN))
    for v in ['med200','med500','lin350']:
        p=prof(arr,t,v)
        for n in 'NS': e[n]=np.maximum(e[n],abs(p[n]-c[n]))
        e['NS']=np.maximum(e['NS'],abs(p['N']/p['S']-c['N']/c['S']))
    return e
P=prof('Xb4b',0.3); Es=stat('Xb4b',0.3); Eb=base('Xb4b',0.3)
Eerr={n:np.hypot(Es[n],Eb[n]) for n in 'NS'}
P4=prof('Xb4d',1.05); Es4=stat('Xb4d',1.05); Eb4=base('Xb4d',1.05)
PB=prof('Xb4b',1.05)
ns=P['N']/P['S']; ns_stat=ns*np.hypot(Es['N']/P['N'],Es['S']/P['S'])
ns_base=Eb['NS']; ns_arr=abs(PB['N']/PB['S']/(P4['N']/P4['S'])-1)*ns
ns_res=np.where(mid<0.5,0.035*ns,0.)   # core-only (r<0.5") bias estimate, N/S biased <=3.5% low (s414); not measured beyond
# team-a q64 (0.3" circular, windows A,B,C median [min,max], HW=350 km/s); their annuli 0-.1,.1-.25,.25-.5,.5-1,1-2,2-3
R=[r for r in csv.DictReader(open('/snapshots/team-a-s33/out/q64r2/img_ratios.csv')) if r['ous']=='Xb4b' and r['HW']=='350' and r['beam']=='0.3' and r['win'] in 'ABC']
AA=[(0,.1),(.1,.25),(.25,.5),(.5,1),(1,2),(2,3)]; amid=np.array([(a+b)/2 for a,b in AA])
CN={(0,.1):'0_0.1',(.1,.25):'0.1_0.25',(.25,.5):'0.25_0.5',(.5,1):'0.5_1.0',(1,2):'1.0_2.0',(2,3):'2.0_3.0'}
TA={n:np.array([[float(r[f'{n}_a'+CN[(a,b)]]) for r in R] for a,b in AA]) for n in 'NS'}
TAm={n:np.median(TA[n],axis=1) for n in 'NS'}; TAlo={n:TA[n].min(1) for n in 'NS'}; TAhi={n:TA[n].max(1) for n in 'NS'}
cN,cS='#D55E00','#0072B2'
fig,(ax,ax2)=plt.subplots(2,1,figsize=(7.2,8),sharex=True,gridspec_kw={'height_ratios':[3,1.7],'hspace':0.06})
def band(a):
    a.axvspan(0.04,0.15,color='0.85',zorder=0); 
band(ax);band(ax2)
for n,c,m in [('N',cN,'o'),('S',cS,'s')]:
    ax.errorbar(mid,P[n],Eerr[n],fmt=m+'-',color=c,ms=6,lw=1.2,capsize=2,label=f'{n}  Xb4b 0.30" (stat$\\oplus$baseline/window)',zorder=3)
    ax.errorbar(mid*1.0,P4[n],np.hypot(Es4[n],Eb4[n]),fmt=m,mfc='none',color=c,ms=7,capsize=2,ls='none',label=f'{n}  Xb4d 1.05" beam (open)',zorder=2)
    d=0.93 if n=='N' else 1.07
    ax.errorbar(amid*d,TAm[n],[TAm[n]-TAlo[n],TAhi[n]-TAm[n]],fmt={'N':'^','S':'v'}[n],color=c,ms=3.5,alpha=.7,lw=.8,capsize=0,label=f'{n}  team-a q64 0.30" (range over windows A-C)',zorder=2)
ax.set_ylabel('HCN/HCO$^+$ 3-2 ratio (moment-0, same beam)'); ax.set_ylim(0.15,1.25)
ax.legend(fontsize=7,loc='upper right',frameon=False)
ax.text(0.06,0.165,'r<0.15" unsupported\n(0.25" native beam)',fontsize=7.5,color='0.3',va='bottom')
ax.annotate('',xy=(0.25/2,0.2),xytext=(0.0,0.2))
# lower
x=np.arange(len(ANN))
ax2.errorbar(mid,ns,ns_stat,fmt='ko',ms=5,capsize=2,label='N/S Xb4b 0.30" (stat)')
ax2.axhline(1,color='0.5',lw=.6)
ax2.set_ylabel('N/S ratio'); ax2.set_ylim(0,3.6); ax2.set_xscale('log'); ax2.set_xlim(0.055,3.6)
# budget bars: stacked-in-quadrature components drawn as horizontal ticks beside points
w=0.12
for i,(m_,s,b,a,rr) in enumerate(zip(mid,ns_stat,ns_base,ns_arr,ns_res)):
    xs=[m_*1.18,m_*1.18*1.0]
    for k,(v,c) in enumerate([(s,'#000000'),(b,'#009E73'),(a,'#CC79A7'),(rr,'#E69F00')]):
        ax2.plot([m_*(1.12+0.07*k)]*2,[ns[i],ns[i]+v] if v>0 else [ns[i]]*2,color=c,lw=3.2,solid_capstyle='butt')
from matplotlib.lines import Line2D
ax2.legend(handles=[Line2D([],[],marker='o',color='k',ls='none',ms=5,label='N/S, error bar = stat'),
  *[Line2D([],[],color=c,lw=3,label=l) for c,l in [('k','budget: stat'),('#009E73','baseline/window'),('#CC79A7','array (Xb4b/Xb4d at 1.05")'),('#E69F00','residual/restored (r<0.5" only)')]]],fontsize=7,ncol=2,loc='upper right',frameon=False)
ax2.set_xlabel('projected radius from nucleus  (arcsec)')
ax2.set_xticks([.1,.2,.3,.5,1,2,3]); ax2.set_xticklabels(['0.1','0.2','0.3','0.5','1','2','3'])
top=ax.secondary_xaxis('top',functions=(lambda r:r*194,lambda p:p/194)); top.set_xlabel('radius (pc), 1" = 194 pc')
top.minorticks_off(); top.set_xticks([20,50,100,200,400,600]); top.set_xticklabels(['20','50','100','200','400','600'])
# native beam marker: bar of width 0.25"(FWHM) in radius
for a in (ax,ax2):
    pass
ax.errorbar([0.0805],[0.62],xerr=[[0.0],[0.0]])  # placeholder removed below
ax.lines[-1].remove()
ax.plot([0.25,0.25],[0.2,0.24],color='k',lw=0)  # no-op
ax.annotate('',xy=(0.32,0.30),xytext=(0.07,0.30),arrowprops=dict(arrowstyle='<->',color='k',lw=1))
ax.text(0.15,0.31,'native beam FWHM 0.25"',fontsize=7,ha='center',va='bottom')
ax.set_xscale('log')
fig.savefig('f1_hcn_hcop_profiles.png',dpi=180,bbox_inches='tight'); fig.savefig('f1_hcn_hcop_profiles.pdf',bbox_inches='tight')
# CSV
with open('f1_values.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['series','annulus_arcsec','r_mid_arcsec','r_mid_pc','nucleus','ratio','err_stat','err_baseline','err_total_stat+baseline','n_s','n_s_stat','n_s_baseline','n_s_array','n_s_residual_core_only'])
    for n in 'NS':
        for i in range(len(ANN)):
            w.writerow([f'Xb4b_0.30',lab[i],mid[i],mid[i]*194,n,round(P[n][i],4),round(Es[n][i],4),round(Eb[n][i],4),round(Eerr[n][i],4),round(ns[i],3),round(ns_stat[i],3),round(ns_base[i],3),round(ns_arr[i],3),round(ns_res[i],3)])
    for n in 'NS':
        for i in range(len(ANN)):
            w.writerow(['Xb4d_1.05',lab[i],mid[i],mid[i]*194,n,round(P4[n][i],4),round(Es4[n][i],4),round(Eb4[n][i],4),round(np.hypot(Es4[n][i],Eb4[n][i]),4)]+['']*5)
    for n in 'NS':
        for i,(a,b) in enumerate(AA):
            w.writerow(['teamA_q64_0.30_median_ABC',f'{a:g}-{b:g}',amid[i],amid[i]*194,n,round(TAm[n][i],4),'', '', f'range {TAlo[n][i]:.4f}-{TAhi[n][i]:.4f}'])
for i in range(len(ANN)): print(lab[i],P['N'][i].round(3),P['S'][i].round(3),ns[i].round(2),ns_stat[i].round(2),ns_base[i].round(2),ns_arr[i].round(2))
