# Q5 round 2: CH3OH 254.015 / HCO+ 3-2 / CS 7-6 in a matched 4.31x2.79" beam (Xb4d smoothed per channel to the Xb47 CS beam, s97a Xb4t; CS from Xb47 spw22 at the same sky positions: pixel of the CS beam
# map = peak-pixel value of the common-beam Xb47 moment map from s98 is NOT used; here CS W is the Xb47 aperture matched-filter W (b7_fits, 3.3" circle), a slightly different (aperture-averaged) estimator, flagged in the output).
# K km/s ratios; noise = max(formal, empirical sliding-template MAD).   usage: python3 -I s97d_q5_ch3oh_cs_matched.py
import sys; sys.path.insert(0,'/workspace/scripts')
import s91_q5_lib as L
from s91_q5_lib import *
L.REFLINE['Xb4t']=('spw25',267.55763); L.HCNLINE['Xb4t']=('spw27',265.88618)
O='Xb4t'; CH3=('CH3OH 2_0-1_-1 E',254.01538); SO2=('SO2 6(3,3)-6(2,4)',254.2805); SIS=('SiS 14-13',254.1031); NS=('NS 11/2-9/2 (f)',253.9685)
B=pd.read_csv('/workspace/out/q5/b7_fits.csv'); F=json.load(open('/workspace/out/q5/field_noise.json')) if False else None
import json; F=json.load(open('/workspace/out/q5/field_noise.json'))
TH=(4.312*2.788)
K=lambda S,nu: 1.222e6*S/(nu**2*TH)
rows=[]
for reg in ('N','S','Bridge','Sblue','W'):
    shp,vc,fw,snr=line_shape(O,reg,mode='profile',tname='ref')
    rr=joint_fit(O,'spw29',reg,[SO2,CH3,SIS,NS],shp,vc,window=1100,order=1,center=SO2[1]); W,sf=rr[CH3[0]]
    ms=slide_sigma(O,'spw29',reg,CH3[1],shp,vc,window=500); sig=max(sf,ms[0])
    h=joint_fit(O,'spw25',reg,[('HCO+ 3-2',267.55763)],shp,vc,window=1100,order=1)['HCO+ 3-2'][0]
    c=B[(B.region==reg)&(B.model=='CS1')&(B['shape']=='profile')].iloc[0]; cs=c.W; ecs=max(c.sig,F[f'Xb47|CS 7-6|{reg}']['mad']); cs_corr=cs-F[f'Xb47|CS 7-6|{reg}']['median']
    k=lambda x: K(x,251.7); kc=lambda x: K(x,339.77); kh=lambda x: K(x,265.1)
    rows.append(dict(region=reg,W_CH3OH_Jy=W,sig=sig,snr=W/sig,CH3OH_K=k(W),eCH3OH_K=k(sig),HCO_K=kh(h),CS_K=kc(cs_corr),eCS_K=kc(ecs),CH3OH_over_HCO=k(W)/kh(h),eCH3OH_over_HCO=k(sig)/kh(h),CH3OH_over_CS=k(W)/kc(cs_corr),eCH3OH_over_CS=np.hypot(k(sig)/k(W) if W!=0 else np.nan,ecs/cs_corr)*abs(k(W)/kc(cs_corr)),UL3_over_HCO=k(max(W,0)+3*sig)/kh(h),UL3_over_CS=k(max(W,0)+3*sig)/kc(cs_corr)))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/ch3oh_cs_matched_4p3x2p8.csv',index=False); pd.set_option('display.width',250); print(df.round(4).T.to_string())
