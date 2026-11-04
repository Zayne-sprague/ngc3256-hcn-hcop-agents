"""F3a. Literature values: Nishimura+2024 (A&A 686, A48) Tables 1,2,4 (resolved HCN/HCO+ 3-2 = beam-sized aperture, L' ratio = T_b ratio), transcribed in
snapshots/team-b-s28/scripts/s300_q46_lit_context.py and re-checked here against the paper text in the library (lit46/nish_full.txt); errors added from Table 4.
CON flags: Nishimura Table 1 col (8) (from Falstad+2021 'Paper I'); HCN-vib detections: Falstad+2021 Table 3 and Nishimura Sect. 2/4.1 text.
Izumi+2016 (ApJ 818, 42) Table 1 high-resolution sample. NGC 3256 values: p0031/p0055 (3-2), p0055 and p0053 (4-3)."""
import csv,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
names=["IRAS 17208-0014","IRAS 09022-3615","IRAS 13120-5453","IRAS F14378-3651","IRAS F05189-2524","IRAS F10565+2448","IRAS 19542+1110","ESO 148-IG002","NGC 6240","IRAS F17138-1017","IRAS 17578-0400","ESO 173-G015","UGC 11763","NGC 3110","IC 4734","NGC 5135","ESO 221-IG10","IC 5179","UGC 2982","NGC 2369","ESO 286-G035","ESO 320-G030","NGC 5734"]
DL=[183,253,135,295,176,184,264,185,106,77.3,60.0,32.7,265,74.3,67.8,52.6,59,47.3,70.6,44.6,73.2,36,60.2]
aAGN=[.05,.19,.04,.05,.60,.04,.08,.19,.11,.07,.03,.03,np.nan,.10,.07,.24,.07,.10,.11,.10,.06,.03,.14]
beam=["0.35x0.32","0.33x0.28","0.36x0.32","0.29x0.26","0.42x0.27","0.33x0.22","0.34x0.29","0.40x0.31","0.32x0.26","0.38x0.29","0.36x0.29","0.33x0.30","0.35x0.31","0.40x0.27","0.38x0.29","0.33x0.27","0.34x0.31","0.28x0.23","0.36x0.30","0.28x0.25","0.29x0.26","0.34x0.29","0.31x0.30"]
res_pc=["298x272","382x324","225x200","388x347","343x220","281x187","423x360","344x267","161x130","139x106","102x82","51x47","397x351","142x96","126x96","83x68","95x87","63x52","117x98","59x53","101x90","62x53","89x86"]
R=[1.738,.793,2.671,1.323,1.630,1.320,1.789,.445,.691,.629,1.340,1.132,.728,1.490,1.076,1.028,.771,.564,np.nan,1.489,.666,2.344,np.nan]
eR=[.019,.009,.064,.022,.012,.035,.018,.015,.015,.123,.004,.010,.283,.330,.033,.147,.164,.127,np.nan,.110,.178,.018,np.nan]
G=[1.464,.555,1.706,1.140,1.534,.878,1.825,.387,.693,.587,1.283,.865,np.nan,1.044,1.063,1.100,.651,.842,np.nan,1.133,.587,2.254,np.nan]
CON={'IRAS 17208-0014','IRAS 17578-0400','ESO 320-G030'}                    # Nishimura Table 1 col (8)
VIB={'IRAS F14378-3651':'Falstad+21 Table 3: detected (not CON in Nishimura Table 1)','ESO 173-G015':'Falstad+21 Table 3 (1.2+-0.2); Nishimura Sect. 4.1: detected, not CON',
     'IRAS F05189-2524':'Nishimura Sect. 2 (faint HCN-vib, below CON criterion; Falstad+2019) -- second-hand'}
UNK={'NGC 6240'}   # Nishimura Sect. 2: broad HCO+ complicates HCN-vib evaluation
def vclass(n):
    if n in CON: return 'CON'
    if n in VIB: return 'vib'
    if n in UNK: return 'unk'
    return 'none'
def acl(a):
    return 'unknown' if np.isnan(a) else ('<0.10' if a<0.10 else ('0.10-0.25' if a<0.25 else '>=0.25'))
rows=[dict(n=n,D=D,a=a,R=r,e=e,G=g,beam=b,pc=p) for n,D,a,r,e,g,b,p in zip(names,DL,aAGN,R,eR,G,beam,res_pc)]
det=[r for r in rows if not np.isnan(r['R'])]
det.sort(key=lambda r:r['R'])
vals=np.array([r['R'] for r in det]); print(len(det))
col={'<0.10':'#0072B2','0.10-0.25':'#E69F00','>=0.25':'#D55E00','unknown':'#888888'}
mk={'CON':'*','vib':'D','none':'o','unk':'s'}; ms={'CON':15,'vib':7,'none':7,'unk':7}
ours={'N 0.28" (55 pc)':(0.99,'#CC0000','-'),'N 1.29" (250 pc)':(0.69,'#CC0000','--'),'S 0.28" (55 pc)':(0.44,'#000000','-'),'S 1.29" (250 pc)':(0.37,'#000000','--')}
def pct(x): return (vals<x).mean()*100, ((vals<x).sum())
pc={k:pct(v[0]) for k,v in ours.items()}
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,7.5),gridspec_kw={'width_ratios':[2.3,1],'wspace':0.3})
for i,r in enumerate(det):
    c=vclass(r['n']); a=acl(r['a'])
    ax.errorbar(r['R'],i,xerr=r['e'],fmt=mk[c],color=col[a],ms=ms[c],mec='k' if c=='CON' else col[a],mew=.6,elinewidth=1,capsize=1.5,zorder=3)
# NGC 5734 lower limit
j=len(det); ax.annotate('',xy=(3.2,j),xytext=(2.5,j),arrowprops=dict(arrowstyle='->',color=col['0.10-0.25'],lw=1.5)); 
for k,(v,c,ls) in ours.items():
    ax.axvline(v,color=c,ls=ls,lw=1.4,zorder=1)
ax.set_xscale('log'); ax.set_xlim(0.3,3.4); ax.set_ylim(-1,len(det)+1); ax.set_yticks(range(len(det)+1)); ax.set_yticklabels([r['n'] for r in det]+['NGC 5734 (HCN only)'],fontsize=7.5)
ax.set_xticks([0.4,0.6,1,1.5,2,3]); ax.set_xticklabels(['0.4','0.6','1','1.5','2','3']); ax.minorticks_off()
ax.set_xlabel("HCN/HCO$^+$ J=3-2 brightness-temperature ratio (resolved, beam-sized aperture)")
ax.set_title("Nishimura+2024: 21 U/LIRGs with both lines detected, ~0.3\" beams (50-420 pc), ranked",fontsize=9.5)
hl=[Line2D([],[],color=c,ls=ls,lw=1.4,label=f'NGC 3256 {k}: {v:.2f}, {pc[k][0]:.0f}th pct ({pc[k][1]}/21 below)') for k,(v,c,ls) in ours.items()]
leg1=ax.legend(handles=hl,fontsize=7.2,loc='upper left',frameon=True,framealpha=.9,title='NGC 3256 (T$_B$ ratio, this work)',title_fontsize=7.5); ax.add_artist(leg1)
h=[Line2D([],[],marker='o',ls='',color=col[a],label=f'$\\alpha_{{AGN}}$ {a}') for a in col]+[Line2D([],[],marker='*',ls='',color='w',mec='k',mfc='0.7',ms=13,label='CON (Table 1)'),Line2D([],[],marker='D',ls='',color='0.5',ms=6,label='HCN-vib detected, not CON'),Line2D([],[],marker='o',ls='',color='0.5',ms=6,label='HCN-vib not detected'),Line2D([],[],marker='s',ls='',color='0.5',ms=6,label='HCN-vib unclear (NGC 6240)')]
ax.legend(handles=h,fontsize=7,loc='lower right',frameon=True,framealpha=.9)
# ---- 4-3 panel
bx.axvspan(0.24,1.03,color='#56B4E9',alpha=.3); bx.axvspan(1.11,3.19,color='#D55E00',alpha=.2)
bx.text(0.50,5.35,'SB\n0.24-1.03',ha='center',fontsize=7.5,va='bottom'); bx.text(1.9,5.35,'AGN\n1.11-3.19',ha='center',fontsize=7.5,va='bottom')
pts=[('team-b p0055: 4.1"x2.7" (ACA, 800x520 pc)',0.52,0.29,0.0,0.0,4),('team-a p0053: 1.0" circ. (188 pc)',0.81,0.22,0.02,0.02,3),('team-a p0053: 1.9" circ. (357 pc)',0.67,0.18,0.0,0.0,2)]
for lab,n,s,en,es,y in pts:
    bx.errorbar(n,y,xerr=en,fmt='o',color='#CC0000',ms=7,capsize=2); bx.errorbar(s,y,xerr=es,fmt='o',color='k',ms=7,capsize=2)
    bx.text(0.17,y+0.33,lab,fontsize=6.6,va='bottom',ha='left')
bx.errorbar(0.95,1,xerr=[[0.14],[0.17]],fmt='o',mfc='white',color='#CC0000',ms=7,capsize=2); bx.errorbar(0.29,1,xerr=[[0.04],[0.06]],fmt='o',mfc='white',color='k',ms=7,capsize=2)
bx.text(0.17,1.33,'team-a p0053: 56 pc core, EXTRAPOLATED',fontsize=6.6,va='bottom',ha='left')
bx.set_xscale('log'); bx.set_xlim(0.15,3.6); bx.set_ylim(0.4,6.4); bx.set_yticks([]); bx.set_xticks([0.2,0.3,0.5,1,2,3]); bx.set_xticklabels(['0.2','0.3','0.5','1','2','3']); bx.minorticks_off()
bx.set_xlabel('HCN/HCO$^+$ J=4-3 T$_B$ ratio'); bx.set_title('Izumi+2016 high-resolution (<500 pc)\nbands vs NGC 3256 4-3 (different beams!)',fontsize=9)
bx.legend(handles=[Line2D([],[],marker='o',ls='',color='#CC0000',label='N'),Line2D([],[],marker='o',ls='',color='k',label='S')],fontsize=7.5,loc='lower right',frameon=False)
fig.savefig('f3a_literature_context.png',dpi=180,bbox_inches='tight'); fig.savefig('f3a_literature_context.pdf',bbox_inches='tight')
# ---- CSV
with open('f3a_literature_table.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['galaxy','D_L_Mpc','beam_arcsec','beam_pc','alpha_AGN','alphaAGN_class','CON_or_vib','HCN/HCO+_3-2_resolved','err','global','plotted','source_and_status'])
    for r in sorted(rows,key=lambda r:(np.isnan(r['R']),r['R'])):
        c=vclass(r['n']); lab={'CON':'CON','vib':'HCN-vib detected, non-CON','none':'HCN-vib not detected','unk':'unclear'}[c]
        src=f"Ratios: Nishimura+24 Table 4 (resolved/global); beam: Table 2; D_L, alpha_AGN, CON: Table 1 (alpha_AGN from Diaz-Santos+2017, CON from Falstad+2021: second-hand)"
        if c=='vib': src+='; vib: '+VIB[r['n']]
        if c=='none': src+='; vib non-detection: '+('Nishimura Sect. 2 text' if r['n'] in('IRAS F10565+2448','IRAS 19542+1110','ESO 148-IG002','UGC 11763') else 'Falstad+21 Table 3 upper limit')
        R_=r['R'];
        w.writerow([r['n'],r['D'],r['beam'],r['pc'],'' if np.isnan(r['a']) else r['a'],acl(r['a']),lab,'' if np.isnan(R_) else R_,'' if np.isnan(r['e']) else r['e'],'' if np.isnan(r['G']) else r['G'],'yes (ranked)' if not np.isnan(R_) else ('arrow: >2.5 (HCN only)' if r['n']=='NGC 5734' else 'no (neither line detected)'),src])
    for lab,rng in [('Izumi+16 SB high-res (M82, NGC253, NGC1097 ring, NGC1614, NGC7469 B,C,D)','0.24-1.03 (values .41 1.03 .82 .24 .75 .45 .48)'),('Izumi+16 AGN high-res (NGC1068 x5, NGC1097, NGC7469)','1.11-3.19 (values 1.53 2.84 3.19 3.14 2.58 2.01 1.11)')]:
        w.writerow([lab,'','','','','','','','','','band','Izumi+2016 Table 1 col R_HCN/HCO+ (4-3), read from library text; '+rng])
for k,(p,nb) in pc.items(): print(k,f'{p:.0f}th percentile, {nb}/21 below')
print('with NGC 5734 limit as 22nd:',{k:f'{(nb/22*100):.0f}' for k,(p,nb) in pc.items()})
