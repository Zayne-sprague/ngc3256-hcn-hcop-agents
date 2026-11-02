# Q5: CS 7-6/HCN 4-3 and SO limits from the single peak-pixel spectrum (native resolution, one beam) at the fixed region positions (Xb47). Same estimator as s91a with kind='pix'.
import sys; sys.path.insert(0,'/workspace/scripts')
from s91_q5_lib import *
rows=[]
for reg in ('N','S','Bridge','E1'):
    shp,vc,fw,snr=line_shape('Xb47',reg,kind='pix',mode='profile',tname='ref')
    hc,vh,fh,_=line_shape('Xb47',reg,kind='pix',mode='profile',tname='hcn')
    rh=joint_fit('Xb47','spw16',reg,[('HCN 4-3',354.50547)],hc,vh,window=1100,kind='pix')
    rc=joint_fit('Xb47','spw22',reg,[('CS 7-6',342.88285)],shp,vc,window=1100,kind='pix')
    rs=joint_fit('Xb47','spw20',reg,[('SO 8_8-7_7',344.31061)],shp,vc,window=700,kind='pix')
    nu,y=spec('Xb47','spw22',reg,'pix'); v=vel(nu,342.88285)-vc; pk=np.convolve(y,np.ones(3)/3,'same')[abs(v)<200].max()
    R=rc['CS 7-6'][0]/rh['HCN 4-3'][0]
    rows.append(dict(region=reg,W_HCN=rh['HCN 4-3'][0],W_CS=rc['CS 7-6'][0],sig_CS=rc['CS 7-6'][1],W_SO=rs['SO 8_8-7_7'][0],sig_SO=rs['SO 8_8-7_7'][1],R_Jy=R,R_K=R*(354.50547/342.88285)**2,rms_chan_mJy_CS=rc['_rms'],rms_chan_mJy_SO=rs['_rms'],peak_CS_3ch_mJy=pk,fwhm_ref=fw))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q5/peakpix_b7.csv',index=False); print(df.round(3).to_string())
