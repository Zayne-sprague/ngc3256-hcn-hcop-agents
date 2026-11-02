# Q5 Band 6 (Xb4f 5.7x3.5", Xb4d 0.94x0.74", Xb4b 0.25"): SO 6_5-5_4 (251.826; blended with CH3OH 3(3)-7(3) A ladder), 34SO 6_6-5_5, SO2 6(3,3)-6(2,4) (254.2805),
# CH3OH 2(0)-1(-1) E (254.015), SiS 14-13, NS (253.968), with HCO+ 3-2 and HCN 3-2 reference integrals; same fixed-position aperture for all lines of an OUS.
import sys; sys.path.insert(0,'/workspace/scripts')
from s91_q5_lib import *
REG=['N','S','Bridge','E1','E2','Sblue','W']
d=pd.read_csv('/workspace/out/linelist.csv'); m=d[(d.sp=='CH3OH')&(d.freq_GHz>251.70)&(d.freq_GHz<251.93)].copy(); m['k']=m.freq_GHz.round(3)
m=m.sort_values('sijmu2').groupby('k').tail(1)
def ladder(T): 
    w=(m.freq_GHz**3*m.sijmu2*np.exp(-m.upper_state_energy_K/T)).values; return list(zip(m.freq_GHz.values,w))
SO=('SO 6_5-5_4',251.82577); SO34=('34SO 6_6-5_5',253.2074); SO2=('SO2 6(3,3)-6(2,4)',254.2805); CH3=('CH3OH 2_0-1_-1 E',254.01538); SIS=('SiS 14-13',254.1031); NS=('NS 11/2-9/2 (f)',253.9685)
OUSL=['Xb4f','Xb4d','Xb4b']; MAP={'Xb4f':dict(a='spw22',b='spw20',h='spw16',n='spw18'),'Xb4d':dict(a='spw31',b='spw29',h='spw25',n='spw27'),'Xb4b':dict(a='spw31',b='spw29',h='spw25',n='spw27')}
rows=[]
def add(ous,reg,mode,model,line,W,sig,vc,fw,snr,rms=np.nan,infl=np.nan): rows.append(dict(ous=ous,region=reg,shape=mode,model=model,line=line,W=W,sig=sig,rms_mJy=rms,infl=infl,vcen=vc,fwhm=fw,ref_snr=snr))
for ous in OUSL:
  mp=MAP[ous]
  for reg in REG:
    for mode in ('profile','gauss'):
      try: shp,vc,fw,snr=line_shape(ous,reg,mode=mode,tname='ref')
      except Exception as e: continue
      if snr<6: continue
      try: 
        for nm,spw,f0,tn in (('HCO+ 3-2',mp['h'],267.55763,'ref'),('HCN 3-2',mp['n'],265.88618,'hcn')):
          s2,vc2,fw2,snr2=line_shape(ous,reg,mode=mode,tname=tn); r=joint_fit(ous,spw,reg,[(nm,f0)],s2,vc2,window=1100,order=1); add(ous,reg,mode,'ref',nm,*r[nm],vc2,fw2,snr2)
        # SO 6_5-5_4: (a) SO only; (b),(c) with the CH3OH A-ladder (LTE weights Tex 10, 25 K) as one extra component with free amplitude (total ladder integral)
        r=joint_fit(ous,mp['a'],reg,[SO],shp,vc,window=700,order=1); add(ous,reg,mode,'SO_only',SO[0],*r[SO[0]],vc,fw,snr,r['_rms'],r['_infl'])
        for T in (10,25):
            r=joint_fit(ous,mp['a'],reg,[SO],shp,vc,window=1000,order=1,extra=[(f'CH3OH A-ladder {T}K',ladder(T))],center=SO[1])
            for k in (SO[0],f'CH3OH A-ladder {T}K'): add(ous,reg,mode,f'SO+ladder{T}K',k,*r[k],vc,fw,snr,r['_rms'],r['_infl'])
        # limits (>=0, 3 and 5 sigma) on SO with the ladder free, and SO only
        for ns in (3,5):
            u=upper_limit(ous,mp['a'],reg,[SO],SO[0],shp,vc,nsig=ns,amax=1.0,window=700,order=1)[0]; add(ous,reg,mode,f'UL{ns}sig','SO_only',u,np.nan,vc,fw,snr)
        # other lines in spw20 / spw22 (joint with neighbours)
        r=joint_fit(ous,mp['b'],reg,[SO2,CH3,SIS,NS],shp,vc,window=1100,order=1,center=SO2[1])
        for k in (SO2[0],CH3[0],SIS[0],NS[0]): add(ous,reg,mode,'spw20joint',k,*r[k],vc,fw,snr,r['_rms'],r['_infl'])
        for ns in (3,5): add(ous,reg,mode,f'UL{ns}sig','SO2 6(3,3)-6(2,4)',upper_limit(ous,mp['b'],reg,[SO2,CH3,SIS,NS],SO2[0],shp,vc,nsig=ns,amax=1.0,window=1100,order=1,center=SO2[1])[0],np.nan,vc,fw,snr)
        r=joint_fit(ous,mp['b'],reg,[SO34],shp,vc,window=500,order=1); add(ous,reg,mode,'34SO_only',SO34[0],*r[SO34[0]],vc,fw,snr,r['_rms'],r['_infl'])
        for ns in (3,5): add(ous,reg,mode,f'UL{ns}sig','34SO 6_6-5_5',upper_limit(ous,mp['b'],reg,[SO34],SO34[0],shp,vc,nsig=ns,amax=1.0,window=500,order=1)[0],np.nan,vc,fw,snr)
      except Exception as e: print('fail',ous,reg,mode,repr(e)); 
  print(ous,'done',flush=True)
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/b6_fits.csv',index=False)
