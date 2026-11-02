import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from common import C
from lineset import VSYS
df=pd.read_csv('/workspace/out/line_table_raw.csv')
df['nu_sky']=df.f0*(1-VSYS/C)
df['K_per_Jybeam']=1.222e6/(df.nu_sky**2*df.bmaj*df.bmin)
df['I_Kkms']=df.I_sum*df.K_per_Jybeam; df['eI_Kkms']=df.eI_sum*df.K_per_Jybeam
df['I_lim3sig']=3*df.eI_sum; df['I_lim3sig_Kkms']=df.I_lim3sig*df.K_per_Jybeam
SPWR={('Xb55','spw22'):(163.244,165.110),('Xb55','spw20'):(164.536,166.402),('Xb55','spw16'):(174.607,176.469),('Xb55','spw18'):(176.366,178.227),
 ('Xb47','spw22'):(338.871,340.668),('Xb47','spw20'):(340.572,342.396),('Xb47','spw16'):(350.622,352.423),('Xb47','spw18'):(352.379,354.180)}
def edge_mhz(r):
    if (r.ous,r.spw) in SPWR: lo,hi=SPWR[(r.ous,r.spw)]
    elif r.ous in('Xb4f','Xb4d','Xb4b'):
        t={'spw16':(264.565,266.429),'spw25':(264.565,266.429),'spw18':(262.808,264.672),'spw27':(262.808,264.672),'spw20':(250.738,252.602),'spw29':(250.738,252.602),'spw22':(248.981,250.845),'spw31':(248.981,250.845)}[r.spw]; lo,hi=t
    return 1e3*min(r.nu_sky-lo,hi-r.nu_sky)
df['edge_MHz']=df.apply(edge_mhz,axis=1)
def status(r):
    if r.edge_MHz<100 and r.snr_I<20: return 'edge-flag'
    if r.snr_I>=5 and r.snr_peak>=3: return 'detected'
    if r.snr_I>=3 and r.snr_peak>=2.5: return 'marginal'
    return 'limit'
df['status']=df.apply(status,axis=1)
df['peak_mJy']=df.peak*1e3; df['rms_mJy']=df.rms_chan*1e3
# aperture flux density (Jy km/s): mean spectrum integral (Jy/beam km/s) * Npix/beam_pix ; Npix = pi (1.5/pix)^2
pixs={'Xb55':1.2,'Xb47':0.53,'Xb4f':0.70,'Xb4d':0.15,'Xb4b':0.047}
df['Npix_ap']=np.where(df.kind=='aperture',np.pi*(df.aperture_arcsec/df.ous.map(pixs))**2,np.nan)
df['flux_Jykms']=np.where(df.kind=='aperture',df.I_sum*df.Npix_ap/df.beam_pix,np.nan)
# shared-pixel flag
key=df[df.kind=='peakpix'].groupby(['ous','spw','line','x','y']).region.apply(lambda s:'+'.join(s))
df['shared_with']=df.apply(lambda r: key.get((r.ous,r.spw,r.line,r.x,r.y),''),axis=1)
df.loc[~df.shared_with.str.contains(r'\+'),'shared_with']=''
cols=['ous','band','spw','region','kind','line','f0','nu_sky','status','I_sum','eI_sum','snr_I','I_Kkms','eI_Kkms','I_lim3sig','I_lim3sig_Kkms','edge_MHz','peak_mJy','rms_mJy','snr_peak','v_mom1','W50','g_amp','g_v0','e_g_v0','g_fwhm','e_g_fwhm','flux_Jykms','bmaj','bmin','dv','ra','dec','shared_with','blends','null_mean','eI_stat','eI_emp','n_null']
df[cols].to_csv('/workspace/out/line_table.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',900,'display.max_colwidth',40)
pk=df[(df.kind=='peakpix')]
# detection matrix: status per line x (ous,region)
mat=pk.pivot_table(index=['band','line'],columns=['ous','region'],values='snr_I',aggfunc='first')
mat.round(1).to_csv('/workspace/out/detection_snr_matrix.csv')
print(df.groupby(['band','line','status']).size().unstack(fill_value=0).to_string())
