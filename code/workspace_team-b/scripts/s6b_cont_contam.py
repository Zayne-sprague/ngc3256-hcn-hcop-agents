# Is line emission leaking into the delivered mfs / cont images?  Expected contamination if the lines were included in an all-channel average
# (= W_line[peak pixel]/bandwidth_kms) versus the observed deviation of each per-spw mfs from a power law through all images of the same OUS.
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
c=pd.read_csv('/workspace/out/s6_continuum.csv'); f=pd.read_csv('/workspace/out/s7_line_fits.csv'); f=f[f.spec=='pix']
BW={'Xb47':1.79,'Xb4f':1.82,'Xb55':1.85}   # GHz mfs bandwidth (typical, from the README header table)
rows=[]
for ous in ('Xb47','Xb4f','Xb55','Xb4d','Xb4b'):
    s=c[c.file.str.startswith(ous+'.')&~c.file.str.contains('tapered')]
    for reg in ('N','S'):
        nu=s.nu_GHz.values; y=s[reg+'_peak_mJy_b'].values
        A=np.vstack([np.ones_like(nu),np.log(nu/260.)]).T; coef,*_=np.linalg.lstsq(A,np.log(np.clip(y,1e-3,None)),rcond=None)
        alpha=coef[1]; pred=np.exp(A@coef)
        for (_,r),p in zip(s.iterrows(),pred):
            spw=r.file.split('.')[1]
            if 'mfs' not in r.file: exp_c=0.0; tag='cont'
            else:
                w=f[(f.ous==ous)&(f.spw==spw)&(f.region==reg)].W.sum(); bwk=1.8/ (r.nu_GHz) *C; exp_c=w*1e3/bwk; tag='mfs'
            rows.append(dict(ous=ous,region=reg,file=r.file,kind=tag,nu_GHz=r.nu_GHz,peak_mJy_b=r[reg+'_peak_mJy_b'],powerlaw_fit=p,resid_mJy=r[reg+'_peak_mJy_b']-p,alpha_fit=alpha,expected_if_lines_included_mJy=exp_c,rms=r.rms_mJy_b))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/s6b_cont_contam.csv',index=False)
pd.set_option('display.width',250,'display.max_colwidth',45); print(df.round(3).to_string())
