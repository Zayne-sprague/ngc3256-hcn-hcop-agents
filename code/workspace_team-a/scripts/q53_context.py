"""Q53: literature HCN/HCO+ table (values transcribed from papers, see src column) vs NGC 3256 N/S.
Run: python3 scripts/q53_context.py   -> out/q53/lit_table.csv, ratio_vs_class.png, summary.json"""
import csv, json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D_PC_PER_ARCSEC = 188.   # 0.25" = 47 pc in this project (Harada 1" ~ 0.19 kpc)
# (object/region, J, ratio, err, beam_pc, class, units, source)
L = []
def a(*x): L.append(x)
iz='Izumi+2016 ApJ818,42 Table1 (HCN/HCO+ 4-3, Tb scale)'
for n,r,e,b,c in [('NGC1068 AGN pos.',1.53,.34,35,'AGN'),('NGC1068 E-knot',2.84,.60,35,'AGN CND'),('NGC1068 W-knot',3.19,.71,35,'AGN CND'),
  ('NGC1068 CND-N',3.14,.72,35,'AGN CND'),('NGC1068 CND-S',2.58,.19,35,'AGN CND'),('NGC1097 AGN',2.01,.29,94,'AGN'),('NGC7469 AGN',1.11,.13,154,'AGN'),
  ('M82',.41,.12,353,'SB'),('NGC253',1.03,.22,262,'SB'),('NGC1097 SB ring',.82,.17,94,'SB'),('NGC1614',.24,.06,468,'SB'),
  ('NGC7469 SB ring B',.75,.14,154,'SB'),('NGC7469 SB ring C',.45,.08,154,'SB'),('NGC7469 SB ring D',.48,.06,154,'SB'),
  ('NGC1068 APEX 18"',1.85,.42,1257,'AGN'),('NGC1365 APEX',.59,.15,1475,'AGN (SB ring dominated)'),('NGC4418',1.64,.25,790,'buried AGN'),
  ('NGC4945 APEX',.78,.17,454,'AGN (SB ring dominated)'),('IRAS12127-1412',1.58,.58,1667,'buried AGN'),('NGC3256 APEX (whole galaxy)',.39,.11,3265,'SB'),
  ('NGC3628 APEX',.36,.11,672,'SB'),('NGC7552 APEX',.54,.15,1702,'SB'),('IRAS13242-5713',1.23,.29,3667,'SB')]:
    a(n,'4-3',r,e,b,c,'Tb',iz)
im='Imanishi+2023 ApJ950,75 Table4 (Jy km/s)'
for J,r,e in [('2-1',.56,.04),('3-2',.38,.07),('4-3',.31,.05)]: a('NGC1614 (1 kpc beam)',J,r,e,1000,'SB','Jy',im)
a('IRAS13120-5453 nucleus',  '4-3',2.8,None,250,'AGN+nuclear SB (Compton thick)','L\'', 'Privon+2017 ApJ835,213 Sec3.2 (325x180 pc element)')
a('IRAS13120-5453 integrated','4-3',1.77,.01,1500,'AGN+nuclear SB','L\'','Privon+2017 Sec3.2')
a('IRAS13120-5453 off-nucleus','4-3',1.0,None,500,'AGN+nuclear SB','L\'','Privon+2017 Sec3.2')
for k,v in [('AGN-dominated (GOALS, IRAM30m)',1.84),('composite',1.14),('SB-dominated',.88)]: a('Privon+2015 mean '+k,'1-0',v,None,None,k.split()[0],'L\'','Privon+2015 ApJ814,39 abstract (galaxy-integrated, ~kpc)')
a('NGC3256 N','1-0',.961,.006,320,'SB (this system)','K','Harada+2018 Table 13, 1.7"')
a('NGC3256 S','1-0',.645,.013,320,'AGN host (this system)','K','Harada+2018 Table 13, 1.7"')
a('NGC3256 N','1-0',.94,.13,400,'SB (this system)','K','Michiyama+2018 Table 4, 2.1"')
a('NGC3256 S','1-0',.60,.08,400,'AGN host (this system)','K','Michiyama+2018 Table 4, 2.1"')
# Imanishi+2019 ApJS241,19 Table 8: 3-2 ULIRG peak (beam<=500 pc) and integrated ratios (Jy km/s)
pk=[1.9,1.3,1.2,1.7,1.1,1.1,1.9,.6,1.5,1.9,1.5,1.5,1.3,.9,1.0,1.6,2.0,1.0,1.6,1.3,1.0,1.1,.8,1.0,1.4,1.8,.5,1.4,2.3]
ig=[1.7,.9,.8,1.0,1.0,1.0,1.1,1.1,1.6,1.1,1.1,1.5,.9,.9,1.0,1.6,1.2,1.0,1.1,.9,.7,1.1,.7,.9,1.4,2.0,.7,1.4,2.6]
a('ULIRG sample (29 nuclei, median; 16-84% in summary)','3-2',float(np.median(pk)),None,500,'ULIRG mixed','Jy','Imanishi+2019 Table 8 peak (beam<=0.2-0.9")')
# this work (K km/s Tb ratios): out/report.json Q1, Q1b, Q10, Q4c
T='this work (Q1/Q1b/Q10; Tb K km/s)'
for n,J,r,e,b,c in [('NGC3256 N core r<0.25"','3-2',.98,.02,56,'SB'),('NGC3256 S core r<0.25"','3-2',.42,.02,56,'AGN host'),
  ('NGC3256 N','3-2',.81,.03,200,'SB'),('NGC3256 S','3-2',.40,.02,200,'AGN host'),
  ('NGC3256 N','3-2',.52,.03,1220,'SB'),('NGC3256 S','3-2',.42,.03,1220,'AGN host'),
  ('NGC3256 N','2-1',.69,.05,1300,'SB'),('NGC3256 S','2-1',.49,.1,1300,'AGN host'),
  ('NGC3256 N deblended','4-3',.47,.01,620,'SB'),('NGC3256 S deblended','4-3',.25,.02,620,'AGN host'),
  ('NGC3256 N 10"','4-3',.33,.03,1880,'SB'),('NGC3256 S 10"','4-3',.30,.03,1880,'AGN host')]:
    a(n,J,r,e,b,c,'Tb',T)
with open('/workspace/out/q53/lit_table.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['object','J','ratio','err','beam_pc','class','units','source']); w.writerows(L)
# unit conversion Jy->Tb for same J
for J,(nh,no) in {'1-0':(88.6316,89.1885),'3-2':(265.886,267.558),'4-3':(354.505,356.734)}.items():
    print(J,'T-ratio / Jy-ratio =',round((no/nh)**2,4))
# 4-3 SB/AGN statistics from Izumi high-res (<500 pc)
hi=[x for x in L if x[7]==iz and x[4]<500]
sb=[x[2] for x in hi if x[5]=='SB']; ag=[x[2] for x in hi if x[5].startswith('AGN') and 'SB' not in x[5]]
# core-scale 4-3 prediction: scale 3-2 core by deblended 4-3/3-2 ratio-of-ratios (0.47/0.66, 0.25/0.36)
pN=.98*.47/.66; pS=.42*.25/.36
S=dict(izumi_highres_SB_range=[min(sb),max(sb)],izumi_highres_SB_median=float(np.median(sb)),izumi_highres_AGN_range=[min(ag),max(ag)],izumi_highres_AGN_median=float(np.median(ag)),
  imanishi19_peak_median=float(np.median(pk)),imanishi19_peak_16_84=list(np.percentile(pk,[16,84])),frac_peak_gt1=float(np.mean(np.array(pk)>1)),
  pred_core_4_3_N=pN,pred_core_4_3_S=pS)
json.dump(S,open('/workspace/out/q53/summary.json','w'),indent=1); print(json.dumps(S,indent=1))
# figure
fig,ax=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
col={'SB':'#2b6cb0','AGN':'#c0392b'}
for x in L:
    if x[1]!='4-3' or x[4] is None: continue
    c=col['SB'] if x[5].startswith('SB') else col['AGN'] if 'AGN' in x[5] else '#888'
    if 'NGC3256 ' in x[0] and x[7]==T: continue
    ax[0].errorbar(x[4],x[2],x[3] or 0,fmt='o' if x[7]==iz else 's',color=c,ms=5,alpha=.8,lw=1)
for x in L:
    if x[7]==T and x[1]=='4-3': ax[0].plot(x[4],x[2],'*',ms=14,color='k' if 'N' in x[0][7:9] else 'g',label=x[0].replace('NGC3256 ',''))
for P,m,cl in [(pN,'N','k'),(pS,'S','g')]: ax[0].plot(56,P,'*',ms=11,mfc='none',color=cl)
ax[0].set_xscale('log'); ax[0].set_xlabel('beam (pc)'); ax[0].set_ylabel('HCN/HCO+ (Tb ratio)'); ax[0].set_title('J=4-3: literature (blue SB, red AGN) vs NGC 3256')
ax[0].axhline(1,ls=':',c='gray'); ax[0].legend(fontsize=7)
for J,mk in [('3-2','o'),('2-1','^'),('1-0','s')]:
    for n,cl in [('NGC3256 N','k'),('NGC3256 S','g')]:
        pts=[(x[4],x[2],x[3]) for x in L if x[0].startswith(n) and x[1]==J and x[4]]
        for p in pts: ax[1].errorbar(p[0],p[1],p[2],fmt=mk,color=cl,mfc='none' if J!='3-2' else cl)
ax[1].plot([],[],'ko',label='N (black) 3-2 / open: 2-1 tri, 1-0 sq'); ax[1].plot([],[],'go',label='S (green)')
ax[1].set_xscale('log'); ax[1].set_xlabel('beam (pc)'); ax[1].set_title('NGC 3256 N, S: ratio vs beam, J=1-0..3-2'); ax[1].axhline(1,ls=':',c='gray'); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('/workspace/out/q53/ratio_vs_class.png',dpi=130)
