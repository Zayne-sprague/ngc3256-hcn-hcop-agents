# Line contamination of continuum images: (mfs - line-free cont) at N,S vs expected contribution if the spw-averaged line emission were included
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
C=299792.458
ct=pd.read_csv('/workspace/out/continuum_table.csv'); d=np.load('/workspace/out/region_spectra.npz')
SPWMAP={'Xb55':['spw16','spw18','spw20','spw22'],'Xb47':['spw16','spw18','spw20','spw22'],'Xb4f':['spw16','spw18','spw20','spw22'],'Xb4d':['spw25','spw27','spw29','spw31'],'Xb4b':['spw25','spw27','spw29','spw31']}
CONT={'Xb55':'Xb55.spw16_18_20_22.cont.I.pbcor','Xb47':'Xb47.spw16_18_20_22.cont.I.pbcor','Xb4f':'Xb4f.spw16_18_20_22.cont.I.pbcor','Xb4d':'Xb4d.spw25_27_29_31.cont.I.pbcor','Xb4b':'Xb4b.spw25_27_29_31.cont.I.pbcor'}
rows=[]
for ous,spws in SPWMAP.items():
    for src in ['N','S']:
        c=ct[(ct.file==CONT[ous])&(ct.src==src)].iloc[0]
        for spw in spws:
            m=ct[(ct.file==f'{ous}.{spw}.mfs.I.pbcor')&(ct.src==src)].iloc[0]
            fr=d[f'{ous}__{spw}__freq']; s=d[f'{ous}__{spw}__{src}__pk']
            dv=abs(np.median(np.diff(fr)))/fr.mean()*C; bw=len(fr)*dv
            lineavg=np.nansum(s)*dv/bw   # mean over whole spw of continuum-subtracted cube = line flux / bandwidth
            diff=m.peak_Jybeam-c.peak_Jybeam; e=np.hypot(m.rms,c.rms)
            rows.append(dict(ous=ous,src=src,spw=spw,cont_mJy=c.peak_Jybeam*1e3,mfs_mJy=m.peak_Jybeam*1e3,diff_mJy=diff*1e3,e_diff_mJy=e*1e3,diff_sigma=diff/e,
                             lineavg_if_included_mJy=lineavg*1e3,frac_of_lineavg=diff/lineavg if lineavg>1e-6 else np.nan))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/continuum_line_contamination.csv',index=False)
pd.set_option('display.width',200,'display.max_rows',200)
print(df.round(3).to_string())
