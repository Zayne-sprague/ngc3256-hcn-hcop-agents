# Q3 step 1: blue/red wing moment maps of HCN and HCO+ 3-2 (Xb4d, matched 1.05" beam), values at fixed positions incl. Michiyama+2018 outflow positions, peak lists with null-window peak counts.
import sys; sys.path.insert(0,'/workspace/scripts')
from q3_lib import *
import pandas as pd
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def off(ra_s,de_s):   # sexagesimal h m s / d m s -> offsets from (RA0,DE0) in arcsec
    h,m,s=ra_s; d,mm,ss=de_s; ra=15*(h+m/60+s/3600); de=-(d+mm/60+ss/3600)
    return ((ra-RA0)*np.cos(np.deg2rad(DE0))*3600,(de-DE0)*3600)
# Michiyama+2018 Table 1 (J2000)
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15),'Soutflow':(-0.21,-4.48),
 'M18 South-OF-Blue':off((10,27,51.17),(43,54,21.43)),'M18 OF-Red':off((10,27,51.20),(43,54,15.64)),'M18 North-OF-Blue':off((10,27,51.02),(43,54,13.12)),
 'M18 North-Nucleus':off((10,27,51.23),(43,54,14.00)),'M18 South-Nucleus':off((10,27,51.22),(43,54,19.20)),
 'mirror of Soutflow N of S':(-0.15+0.06,-2.55+1.93)}
W={'blue 200-350':(VS-350,VS-200),'blue 200-450':(VS-450,VS-200),'red 200-350':(VS+200,VS+350),'core +-150':(VS-150,VS+150)}
C_={sp:cube(sp,fw=1.05) for sp in ('HCN','HCO+')}
cT=Tb_factor(LINES['HCN 3-2@Xb4d'][2]/(1+VS/C),1.05)/Tb_factor(LINES['HCO+ 3-2@Xb4d'][2]/(1+VS/C),1.05)
rows=[];MAPS={};PK=[]
for wn,(lo,hi) in W.items():
    for sp,c in C_.items():
        m0,n=wsum(c,lo,hi); sg,nw=wnoise(c,lo,hi); MAPS[(wn,sp)]=(m0,sg)
        r=np.hypot(c['dx'],c['dy']); sn=m0/sg
        lm=(m0==ndi.maximum_filter(m0,size=int(1.5/c['pix'])))&(sn>5)&(r<10)
        for y,x in zip(*np.where(lm)): PK.append(dict(window=wn,line=sp,dRA=c['dx'][y,x],dDec=c['dy'][y,x],I=m0[y,x],snr=sn[y,x]))
        for pn,(x,y) in POS.items():
            rows.append(dict(window=wn,line=sp,pos=pn,x=x,y=y,I=samp(c,m0,x,y),e=samp(c,sg,x,y),n_null_windows=nw))
T=pd.DataFrame(rows); T['snr']=T.I/T.e
P=T.pivot_table(index=['window','pos'],columns='line',values=['I','e','snr']); P.columns=['_'.join(a) for a in P.columns]; P=P.reset_index()
R=[];L=[]
for _,r in P.iterrows():
    a,b,ea,eb=r['I_HCN'],r['I_HCO+'],r['e_HCN'],r['e_HCO+']
    if a/ea>=3 and b/eb>=3: R.append(a/b*cT);L.append('value')
    elif a/ea>=3: R.append(a/(3*eb)*cT);L.append('lower limit (HCO+<3sig)')
    elif b/eb>=3: R.append(3*ea/b*cT);L.append('upper limit (HCN<3sig)')
    else: R.append(np.nan);L.append('both <3 sigma')
P['HCN/HCO+ (Tb)']=R;P['kind']=L
P['eR']=np.where(P.kind=='value',P['HCN/HCO+ (Tb)']*np.hypot(P.e_HCN/P.I_HCN,P['e_HCO+']/P['I_HCO+']),np.nan)
pd.set_option('display.width',250); P.round(3).to_csv('/workspace/out/q3/q3_wing_positions.csv',index=False); print(P.round(3).to_string())
pk=pd.DataFrame(PK); pk.round(3).to_csv('/workspace/out/q3/q3_wing_peaks.csv',index=False); print(pk.round(2).to_string())
# null-window peak counts: same peak finder on same-width sums in line-free windows (HCN: only blue side 1380-2100 is free)
nrows=[]
for sp,c in C_.items():
    for wn in ('blue 200-350','red 200-350','blue 200-450'):
        lo,hi=W[wn]; n=int(((c['v']>=lo)&(c['v']<hi)).sum()); v=c['v']; r=np.hypot(c['dx'],c['dy'])
        for i in range(0,len(v)-n,max(1,n//2)):
            if any(v[i]>=a and v[i+n-1]<=b for a,b in FREE[sp]):
                m0=c['a'][i:i+n].sum(0)*c['dv']; sg=wnoise(c,lo,hi)[0]; sn=m0/sg
                lm=(m0==ndi.maximum_filter(m0,size=int(1.5/c['pix'])))&(sn>5)&(r<10)
                nrows.append(dict(line=sp,width=wn,v0=v[i],n_peaks_gt5=int(lm.sum()),max_snr=float(np.max(np.where(r<10,sn,-9))),n_gt8=int((lm&(sn>8)).sum())))
N=pd.DataFrame(nrows);N.to_csv('/workspace/out/q3/q3_null_windows.csv',index=False)
print(N.groupby(['line','width']).agg(nwin=('v0','size'),peaks5=('n_peaks_gt5','sum'),peaks8=('n_gt8','sum'),maxsn=('max_snr','max')).round(2))
fig,ax=plt.subplots(2,4,figsize=(17,8.5))
for j,wn in enumerate(['blue 200-350','core +-150','red 200-350','blue 200-450']):
    for i,sp in enumerate(('HCN','HCO+')):
        m0,sg=MAPS[(wn,sp)]; c=C_[sp]; a=ax[i,j]; ext=[c['dx'].min(),c['dx'].max(),c['dy'].min(),c['dy'].max()]
        a.imshow(m0/sg,origin='lower',extent=ext,vmin=-3,vmax=12 if 'core' not in wn else 60,cmap='RdBu_r' if 'core' not in wn else 'viridis')
        for pn,(x,y) in POS.items():
            if pn in('N','S','Soutflow','M18 OF-Red','M18 South-OF-Blue'): a.plot(x,y,'k+' if pn[0]!='M' else 'gx'); 
        a.set_xlim(8,-8);a.set_ylim(-9,7);a.set_title(f'{sp} {wn} S/N')
plt.tight_layout();plt.savefig('/workspace/out/q3/q3_wing_maps.png',dpi=90)
