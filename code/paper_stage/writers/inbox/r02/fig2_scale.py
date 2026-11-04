"""R02 / Figure 2: HCN/HCO+ 3-2 T_b ratio and peak T_b of N, S (Bridge) vs beam size.
Estimator = that of out/q26/s171_summary.txt (Table 1/3): R_ap = aperture (r = target beam FWHM) mean of moment-0 at the HCO+ 3-2 peak
(window +-350 km/s about 2780 km/s), ratio of T_b (RJ). Peak T_b = 3-channel-smoothed peak channel T_b at the aperture centre pixel ('Tchan_peak' in s170 json).
Inputs: out/q26/s170_Xb4b.json, s170_Xb4d.json, s170_Xb4f.json.  usage: python3 -I fig2_scale.py"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='/workspace/out/q26/'; PC=194.
A=json.load(open(O+'s170_Xb4b.json')); D=json.load(open(O+'s170_Xb4d.json')); F=json.load(open(O+'s170_Xb4f.json'))
TH={'0.26"(50pc)':0.26,'55pc(0.2835")':0.2835,'120pc(0.6186")':0.6186,'250pc(1.2887")':1.2887}
C={'N':'#D55E00','S':'#0072B2','Bridge':'#009E73'}
plt.rcParams.update({'font.size':10})
fig,(a,b)=plt.subplots(1,2,figsize=(12.5,5.2))
xs=[TH[s] for s in TH]
for rn,c in C.items():
    y=np.array([A[s][rn]['R_ap'][0] for s in TH]); e=np.array([A[s][rn]['R_ap'][1] for s in TH]); x=np.array(xs)
    o=np.argsort(x)
    if rn!='Bridge': a.fill_between(x[o],0.86*y[o],1.14*y[o],color='0.6',alpha=0.22,lw=0,zorder=0,label='±14% flux scale' if rn=='N' else None)
    a.errorbar(x,y,e,color=c,marker='o',ms=6,lw=1.3,capsize=2,label=f'{rn}, Xb4b' if rn!='Bridge' else 'Bridge, Xb4b',zorder=3)
    for D_,mk,xx,lab in ((D['250pc(1.2887")'],'s',1.2887,'Xb4d'),(D['6.5"'],'^',6.5,'Xb4d 6.5"'),(F['6.5"'],'v',6.5,'Xb4f 6.5"')):
        a.errorbar([xx*(1.04 if mk=='s' else {'N':0.92,'S':1.0,'Bridge':1.08}[rn])],[D_[rn]['R_ap'][0]],[D_[rn]['R_ap'][1]],color=c,marker=mk,mfc='w',ms=7,lw=0,elinewidth=1.2,zorder=4)
    if rn=='Bridge':
        pass
a.set_xscale('log'); a.set_xlim(0.2,8); a.set_ylim(0,1.2)
a.set_xlabel('beam FWHM (arcsec)'); a.set_ylabel(r'HCN/HCO$^+$ J=3–2 $T_b$ ratio (same velocity window)')
sec=a.secondary_xaxis('top',functions=(lambda x:x*PC,lambda x:x/PC)); sec.set_xlabel('beam FWHM (pc), 1" = 194 pc')
from matplotlib.lines import Line2D
h,l=a.get_legend_handles_labels()
h+= [Line2D([],[],marker='s',mfc='w',mec='k',ls='',label='Xb4d (open square, 250 pc)'),Line2D([],[],marker='^',mfc='w',mec='k',ls='',label='Xb4d 6.5"'),Line2D([],[],marker='v',mfc='w',mec='k',ls='',label='Xb4f 6.5"')]
a.legend(handles=h,fontsize=8,loc='upper right',frameon=False,ncol=2)
a.text(0.215,0.04,'N/S contrast (aperture, stat. err. ~0.03):\n 2.25 (0.26")  2.32 (0.28")  2.11 (0.62")  1.90 (1.29")\n Xb4d: 1.82 (1.29")  1.10 (6.5")   Xb4f: 1.07 (6.5")\n adopted 50–120 pc: 2.1±0.2 (estimator spread 1.9–2.3)',fontsize=8,va='bottom',bbox=dict(fc='w',ec='0.7'))
a.set_title('(a) ratio vs beam: aperture r = beam FWHM',fontsize=10)
# right
for rn in ('N','S'):
    c=C[rn]
    for ln,ls,mk,mfc in (('HCO+','-','o',c),('HCN','--','D','w')):
        y=[A[s][rn][ln]['Tchan_peak'] for s in TH]
        b.plot(xs,y,ls,color=c,marker=mk,mfc=mfc,ms=6,lw=1.3,label=f'{ln} {rn}')
        for D_,xx in ((D['250pc(1.2887")'],1.2887*1.04),(D['6.5"'],6.5),(F['6.5"'],6.5)):
            pass
        # Xb4d/Xb4f as small open/filled markers
        for D_,xx,m2 in ((D['250pc(1.2887")'],1.2887*1.04,'s'),(D['6.5"'],6.5,'^'),(F['6.5"'],6.5,'v')):
            b.plot([xx],[D_[rn][ln]['Tchan_peak']],marker=m2,color=c,mfc='w' if ln=='HCN' else c,ms=6,ls='',alpha=0.8)
b.set_xscale('log'); b.set_yscale('log'); b.set_xlim(0.2,8); b.set_ylim(0.03,4)
b.set_xlabel('beam FWHM (arcsec)'); b.set_ylabel(r'peak $T_b$ (K), 3-channel smoothed, RJ')
sec2=b.secondary_xaxis('top',functions=(lambda x:x*PC,lambda x:x/PC)); sec2.set_xlabel('beam FWHM (pc)')
b.legend(fontsize=8,frameon=False,loc='lower left',ncol=2)
b.text(0.22,3.5,'peak N/S at 0.26": HCO$^+$ %.2f, HCN %.2f'%(A['0.26"(50pc)']['N']['HCO+']['Tchan_peak']/A['0.26"(50pc)']['S']['HCO+']['Tchan_peak'],A['0.26"(50pc)']['N']['HCN']['Tchan_peak']/A['0.26"(50pc)']['S']['HCN']['Tchan_peak']),fontsize=8.5,va='top')
b.set_title('(b) peak $T_b$ (solid: HCO$^+$, dashed/open: HCN)',fontsize=10)
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
for ax_ in (a,b):
    ax_.xaxis.set_major_locator(FixedLocator([0.3,0.5,1,2,3,5])); ax_.xaxis.set_major_formatter(FixedFormatter(['0.3','0.5','1','2','3','5'])); ax_.xaxis.set_minor_locator(NullLocator())
for s_ in (sec,sec2):
    s_.xaxis.set_major_locator(FixedLocator([100,200,500,1000])); s_.xaxis.set_major_formatter(FixedFormatter(['100','200','500','1000'])); s_.xaxis.set_minor_locator(NullLocator())
fig.tight_layout(); fig.savefig('/workspace/paper_answers/r02/fig2_ratio_vs_beam.png',dpi=180)
for rn in ('N','S'): print(rn,[(round(A[s][rn]['HCN']['Tchan_peak'],2),round(A[s][rn]['HCO+']['Tchan_peak'],2)) for s in TH],'Xb4d250',D['250pc(1.2887")'][rn]['HCN']['Tchan_peak'],D['250pc(1.2887")'][rn]['HCO+']['Tchan_peak'])
print({k:(D['6.5"'][k]['R_ap'][0],F['6.5"'][k]['R_ap'][0]) for k in ('N','S','Bridge')}, D['250pc(1.2887")']['N']['R_ap'][0]/D['250pc(1.2887")']['S']['R_ap'][0])
