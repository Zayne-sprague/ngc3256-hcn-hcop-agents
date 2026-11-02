# Standard line table v2 = old table (out/line_table.csv) minus the wrong 'unid 163.54' rows and the superseded CH3OH rows, plus the CH3OH fits (out/ch3oh_table.csv).
import pandas as pd, numpy as np
old=pd.read_csv('/workspace/out/line_table.csv')
old=old[~old.line.str.startswith('unid 163.54')&~old.line.str.startswith('CH3OH')]
new=pd.read_csv('/workspace/out/ch3oh_table.csv')
new=new[~new.line.str.contains('9\\(-0,9\\)')&~new.line.str.contains('7\\(3\\)-7\\(2\\)')]   # 9(0)-8(1) sits on the HCO+ wing (Xb4d S: 16 sigma = wing, not CH3OH); 7(3)-7(2) negative/SO-blended
rows=pd.DataFrame(dict(ous=new.ous,band=new.ous.map({'Xb55':'B5','Xb4f':'B6','Xb4d':'B6','Xb4b':'B6'}),spw=new.spw,region=new.region,kind='peakpix',line=new.line,
  nu_sky=new.nu_sky,status=new.status,I_sum=new.I_Jykms,eI_sum=new.eI_Jykms,snr_I=new.snr,I_Kkms=new.I_Kkms,eI_Kkms=new.eI_Kkms,I_lim3sig=3*new.eI_Jykms,I_lim3sig_Kkms=3*new.eI_Kkms,
  rms_mJy=new.rms_mJy,g_v0=new.v0,g_fwhm=new.fwhm,blends='CH3OH LTE-weighted blend template; v and FWHM fixed to the region HCO+ values (same OUS)'))
out=pd.concat([old,rows],ignore_index=True); out.to_csv('/workspace/out/line_table_v2.csv',index=False); print(len(old),len(rows),len(out))
