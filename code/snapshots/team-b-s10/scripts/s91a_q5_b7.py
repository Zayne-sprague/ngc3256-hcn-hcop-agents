# Q5 Band 7 (Xb47, 4.1-4.3 x 2.7-2.8 arcsec): HCN 4-3, HCO+ 4-3, CS 7-6 (+H2CS 10_0-9_0, 29SiO 8-7 blend), SO 8_8-7_7 (+HC15N), 29SiO -> 28SiO.
import sys; sys.path.insert(0,'/workspace/scripts')
from s91_q5_lib import *
REG=['N','S','Bridge','E1','E2','Sblue','W']
CS=('CS 7-6',342.88285); H2CS=('H2CS 10_0,10-9_0,9',342.94644); SIO29=('29SiO 8-7',342.98066)
SO=('SO 8_8-7_7',344.31061); HC15N=('HC15N 4-3',344.20016)
rows=[]; ous='Xb47'
for reg in REG:
  for mode in ('profile','gauss'):
    shp,vc,fw,snr=line_shape(ous,reg,mode=mode,tname='ref')       # HCO+ 4-3 shape
    ref={}
    for nm,spw,f0,tn in (('HCO+ 4-3','spw18',356.73424,'ref'),('HCN 4-3','spw16',354.50547,'hcn')):
        s2,vc2,fw2,snr2=line_shape(ous,reg,mode=mode,tname=tn)
        r=joint_fit(ous,spw,reg,[(nm,f0)],s2,vc2,window=1100,order=1); ref[nm]=r[nm]; ref[nm+'_fw']=fw2
    fits={}
    fits['CS1']=joint_fit(ous,'spw22',reg,[CS],shp,vc,order=1)
    fits['CS3']=joint_fit(ous,'spw22',reg,[CS,H2CS,SIO29],shp,vc,order=1)
    fits['CS3q']=joint_fit(ous,'spw22',reg,[CS,H2CS,SIO29],shp,vc,order=2)
    fits['SO1']=joint_fit(ous,'spw20',reg,[SO],shp,vc,window=700,order=1)
    fits['SO2']=joint_fit(ous,'spw20',reg,[SO,HC15N],shp,vc,window=900,order=1,center=SO[1])
    # one-sided 3 and 5 sigma limits (amplitude >= 0, others refit) for the blend partners and SO
    lims={}
    for nm,spw,lines,ln,kw in (('H2CS','spw22',[CS,H2CS,SIO29],H2CS[0],{}),('29SiO','spw22',[CS,H2CS,SIO29],SIO29[0],{}),('SO','spw20',[SO,HC15N],SO[0],dict(window=900,center=SO[1])),('SO_only','spw20',[SO],SO[0],dict(window=700))):
        for ns in (3,5): lims[(nm,ns)]=upper_limit(ous,spw,reg,lines,ln,shp,vc,nsig=ns,amax=1.0,order=1,**kw)[0]
    # CS with H2CS and 29SiO both held at their 3-sigma limits (lowest CS compatible with positivity)
    rl=joint_fit(ous,'spw22',reg,[CS,H2CS,SIO29],shp,vc,order=1,fixed={H2CS[0]:lims[('H2CS',3)],SIO29[0]:lims[('29SiO',3)]})
    rows.append(dict(region=reg,shape=mode,model='CS_low(H2CS,29SiO at 3sig UL)',line='CS 7-6',W=rl['CS 7-6'][0],sig=rl['CS 7-6'][1],rms_mJy=rl['_rms'],infl=rl['_infl'],vcen=vc,fwhm=fw,ref_snr=snr))
    for (nm,ns),u in lims.items(): rows.append(dict(region=reg,shape=mode,model=f'UL{ns}sig',line=nm,W=u,sig=np.nan,rms_mJy=np.nan,infl=np.nan,vcen=vc,fwhm=fw,ref_snr=snr))
    for k,r in fits.items():
        for ln in [n for n in r if not n.startswith('_')]:
            W,s=r[ln]; rows.append(dict(region=reg,shape=mode,model=k,line=ln,W=W,sig=s,rms_mJy=r['_rms'],infl=r['_infl'],vcen=vc,fwhm=fw,ref_snr=snr))
    for nm in ('HCO+ 4-3','HCN 4-3'): rows.append(dict(region=reg,shape=mode,model='ref',line=nm,W=ref[nm][0],sig=ref[nm][1],rms_mJy=np.nan,infl=np.nan,vcen=vc,fwhm=ref[nm+'_fw'],ref_snr=snr))
    for lab,spw,f in (('CS 7-6','spw22',CS[1]),('SO 8_8-7_7','spw20',SO[1])):
        mad,sd,n,med=slide_sigma(ous,spw,reg,f,shp,vc)
        rows.append(dict(region=reg,shape=mode,model='slide',line=lab,W=med,sig=mad,rms_mJy=sd,infl=n,vcen=vc,fwhm=fw,ref_snr=snr))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/b7_fits.csv',index=False)
pd.set_option('display.width',250); print(df[df['shape']=='profile'].round(3).to_string())
