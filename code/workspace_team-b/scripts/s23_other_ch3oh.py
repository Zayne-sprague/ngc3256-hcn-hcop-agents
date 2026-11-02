# Other CH3OH lines (Eu<100 K) in the Band 6 windows: predicted strength relative to 254.0154 for LTE T, and measured
# amplitude (matched filter with the same-pixel HCO+ 3-2 profile as template, velocity fixed at HCO+ v0 + 10 km/s) at S, Bridge, N in Xb4f.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import interp1d
from regions import REGIONS
C=299792.458; HCO=267.557626
ll=pd.read_csv('/workspace/out/linelist.csv')
ch=ll[ll.sp.eq('CH3OH')&(ll.upper_state_energy_K<100)&(((ll.freq_GHz>251.4)&(ll.freq_GHz<253.0))|((ll.freq_GHz>253.2)&(ll.freq_GHz<254.8))|((ll.freq_GHz>265.4)&(ll.freq_GHz<267.0))|((ll.freq_GHz>267.2)&(ll.freq_GHz<268.7)))].copy()
ch['f']=ch.freq_GHz.round(3)
ch=ch.drop_duplicates('f'); ch=ch[ch.freq_GHz<267.3]   # 267.46 is in the HCO+ 3-2 wing (267.558): excluded   # +/- parity duplicates of one line
ch['gu']=[int(r.resolved_QNs.strip()[0:2].replace('(','').replace(' ','') or 0) if False else 0 for r in ch.itertuples()]
print('CH3OH E/A lines with Eu<100 K in the Band 6 windows (rest GHz, Eu, log Aij):');print(ch[['resolved_QNs','freq_GHz','upper_state_energy_K','aij']].to_string(index=False))
def load(spw):
    h=fits.open(f'/data/member.uid___A001_X1273_Xb4f.ngc_3256_sci.spw{spw}.cube.I.pbcor.fits',memmap=True)
    H=h[0].header; d=np.squeeze(h[0].data); nu=(H['CRVAL3']+(np.arange(H['NAXIS3'])-H['CRPIX3']+1)*H['CDELT3'])/1e9; return d,H,nu
def sp(d,H,r):
    x,y=WCS(H).celestial.all_world2pix(*REGIONS[r][:2],0); return d[:,int(round(float(y))),int(round(float(x)))]*1e3
spws={'16':load('16'),'18':load('18'),'20':load('20'),'22':load('22')}
def which(f,v=2780):
    nuo=f/(1+ (v/C)/(1-v/C)) if False else f*(1-v/C)
    for k,(d,H,nu) in spws.items():
        if nu.min()+0.02<nuo<nu.max()-0.02: return k
rows=[]
dH,HH,nuH=spws['16']
for r in ('S','Bridge','N'):
    sH=sp(dH,HH,r); vH=C*(1-nuH/HCO); o=((vH>1500)&(vH<2100))|((vH>3500)&(vH<4100)); m=(vH>2200)&(vH<3400)
    prof=sH-np.polyval(np.polyfit(vH[o],sH[o],1),vH); prof/=prof[m].max()
    p0=np.polyfit(vH[m],prof[m],0)
    t=interp1d(vH[::-1],prof[::-1],bounds_error=False,fill_value=0.)
    # HCO+ v0 from the profile peak weighting
    v0=(prof[m].clip(0)*vH[m]).sum()/prof[m].clip(0).sum()
    for f,E,a,q in zip(ch.freq_GHz,ch.upper_state_energy_K,ch.aij,ch.resolved_QNs):
        k=which(f)
        if k is None: continue
        d,H,nu=spws[k]; s=sp(d,H,r); v=C*(1-nu/f); 
        sel=(v>1300)&(v<4300); line=np.abs(v-v0-10)<450
        pf=np.polyfit(v[sel&~line],s[sel&~line],1); res=s[sel&~line]-np.polyval(pf,v[sel&~line]); rms=res.std()
        s2=s-np.polyval(pf,v); tt=t(v-10)  # template shifted so its centre sits at v0+10
        tt=interp1d(vH-v0,prof,bounds_error=False,fill_value=0.)(v-(v0+10)); w=sel&(np.abs(v-v0-10)<500)
        A=(tt[w]*s2[w]).sum()/(tt[w]**2).sum(); eA=rms/np.sqrt((tt[w]**2).sum())
        rows.append(dict(region=r,f=f,Eu=E,logA=a,qn=q.strip(),spw=k,amp_mJy=A,err=eA,snr=A/eA))
df=pd.DataFrame(rows); pd.set_option('display.width',200)
df.to_csv('/workspace/out/s23_ch3oh_other_lines.csv',index=False)

A254=dict(S=df[(df.region=='S')&(df.f.round(3)==254.015)].amp_mJy.iloc[0],Bridge=df[(df.region=='Bridge')&(df.f.round(3)==254.015)].amp_mJy.iloc[0])
def Jq(q):
    import re; m=re.match(r'\s*(\d+)',q); return int(m.group(1)) if m else 0
for T in (8,10,15,25):
    pr=[]
    for r in ch.itertuples():
        if r.freq_GHz>254.1 or r.freq_GHz<254.0 or abs(r.freq_GHz-254.015)<1e-3:
            ratio=((2*Jq(r.resolved_QNs)+1)*10**r.aij*np.exp(-r.upper_state_energy_K/T)/r.freq_GHz**2)/(5*10**-4.72074*np.exp(-20.09/T)/254.0154**2)*(r.freq_GHz/254.0154)**2
            pr.append((round(r.freq_GHz,3),round(ratio*A254['S'],2)))
    print(f'LTE T={T} K: predicted peak (mJy/beam, S pixel; scaled from the 254 line, same width) per line:',pr)
print(df.pivot_table(index=['f','Eu','qn'],columns='region',values='snr').round(1).to_string())
# LTE predictions relative to 254.0154 (W ~ g A exp(-E/T)/nu^2, g=2J+1 for all)
import re
def J(q):
    m=re.match(r'\s*(\d+)',q); return int(m.group(1)) if m else 0
base=(5*10**-4.72074*np.exp(-20.09/10)/254.0154**2)
for T in (10,15,25):
    tot=0
    for r in ch.itertuples():
        if 251.6<r.freq_GHz<251.95 and r.aij>-4.5:
            tot+= (2*J(r.resolved_QNs)+1)*10**r.aij*np.exp(-r.upper_state_energy_K/T)/r.freq_GHz**2
    b=(5*10**-4.72074*np.exp(-20.09/T)/254.0154**2)
    print(f'T={T} K: sum of CH3OH A-type 3(3)-2(...) lines in 251.7-251.92 / 254.015 line (K km/s) = {tot/b:.3f}')
