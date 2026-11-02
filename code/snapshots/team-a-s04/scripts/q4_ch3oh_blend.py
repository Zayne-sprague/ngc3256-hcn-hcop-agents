# Q4f: does blending with SO2 / NH2CHO lines in the Band-5 Q-branch window explain the 1.5-1.6x "excess" of CH3OH Q J=1-3 over the LTE/RADEX expectation?
# Spectrum: Xb55 spw22 (+spw20 for J>=7 not needed), position (RA0,DE0), native 9.36x6.28". Window 165.00-165.45 GHz rest at v=2777 (J=1-5 Q lines).
# Model: baseline(1) + A_Q * sum_J w_J(T) G_J (CH3OH Q J=1..5, LTE weights nu*Smu2*exp(-Eu/T)) + A_SO2 * sum w(T_so2) G (SO2 9(4,6)-10(3,7) 165.1236, 5(2,4)-5(1,5) 165.1447, 7(1,7)-6(0,6) 165.2254; LTE weights, Eu from the curated catalogue)
#        + A_NH2CHO G(165.1766; 8(2,6)-8(1,7) hyperfine sum). Common FWHM, v0 (free within 2777+-25, FWHM 150+-30 = HCO+ 2-1 value); amplitudes >= 0 (NNLS). Profile in A_SO2: chi2 and the Q J=1-3 integral vs A_SO2.
# Output out/q4/ch3oh_blend.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_lib import *
from scipy.optimize import nnls
RA0,DE0=156.9634583333,-43.90461111111
cur=pd.read_csv('/workspace/out/inband_lines_curated.csv')
SO2=[(165.123634,80.63684,2.8055),(165.144651,23.58866,6.55135),(165.225451,27.08443,4.429)]   # f0, Eu, Smu2 (CDMS rows of the curated list)
NH2=[(165.17655729,48.80722,1.14947),(165.1766232,48.80723,1.46083),(165.1767557,48.80723,3.97007)]
Q=[(165.050175,23.3673,5.38627),(165.061130,28.01088,8.94855),(165.099240,34.97606,12.46907),(165.190475,44.2628,15.92557),(165.369341,55.87117,19.29095)]
def comp(v,lines,T,v0,fw,fref=165.050175):
    w=np.array([f*s*np.exp(-E/T) for f,E,s in lines]); t=np.zeros_like(v)
    for (f,E,s),wi in zip(lines,w): t+=wi*np.exp(-0.5*((v-(1-f/fref)*C-v0)/(fw/2.355))**2)
    return t/w.sum()    # unit-area weights -> amplitude is the peak of the summed profile / sum of weights ; integral = amp*fw*1.0645 (per unit total weight)
fr,s,_,nb=spectrum('Xb55','spw22','N',None,pos=(RA0,DE0)); v=vel(fr,165.050175); m=np.abs(v-2777)<650; vv=v[m]; y=s[m]*1e3
dvv=abs(np.median(np.diff(v))); kq=kb(165.05,9.36,6.28)*1e-3
rows=[]
for TQ in (6.,9.,15.):
  for TS in (10.,30.,80.):
    for v0 in (2777.,):
      best={}
      for fw in (130.,150.,170.):
        for dv0 in (-20,0,20):
          B=np.vstack([np.ones_like(vv),(vv-2777)/100,comp(vv,Q,TQ,v0+dv0,fw),comp(vv,SO2,TS,v0+dv0,fw),comp(vv,NH2,20.,v0+dv0,fw)]).T
          for fixso2 in (False,True):
            Bc=B[:,2:] if not fixso2 else B[:,[2,4]]
            # baseline free sign: solve with baseline via augment (nonneg only for line columns): use lstsq for baseline + nnls on residual iterations
            Aall=np.hstack([B[:,:2],-B[:,:2],Bc]); x,rn=nnls(Aall,y)
            chi=rn**2; amp=x[4:]
            key=('noSO2' if fixso2 else 'SO2free')
            if key not in best or chi<best[key][0]: best[key]=(chi,amp,fw,dv0)
      for key,(chi,amp,fw,dv0) in best.items():
        AQ=amp[0]; ASO=amp[1] if key=='SO2free' else 0.; ANH=amp[-1]
        # integrals K km/s: per unit weight-normalised profile integral = amp * fw*1.0645 * (peak of sum of weights is not 1) -> convert via numeric integral of the component
        gridv=np.arange(2300,3300,1.); 
        IQ=np.sum(comp(gridv,Q[:3],TQ,v0+dv0,fw))*1.*AQ*kq*(sum(f*s_*np.exp(-E/TQ) for f,E,s_ in Q[:3])/sum(f*s_*np.exp(-E/TQ) for f,E,s_ in Q))
        ISO=np.sum(comp(gridv,SO2[1:2],TS,v0+dv0,fw))*ASO*kq*((SO2[1][0]*SO2[1][2]*np.exp(-SO2[1][1]/TS))/sum(f*s_*np.exp(-E/TS) for f,E,s_ in SO2)) if ASO>0 else 0.
        rows.append(dict(T_CH3OH=TQ,T_SO2=TS,model=key,chi2=chi,ndata=int(m.sum()),fwhm=fw,dv0=dv0,I_Q123_Kkms=IQ,I_SO2_5_2_4_Kkms=ISO,A_NH2CHO=ANH))
R=pd.DataFrame(rows); R.round(4).to_csv(OUT+'ch3oh_blend.csv',index=False); pd.set_option('display.width',250)
g=R.groupby(['T_CH3OH','model']).agg(chi2=('chi2','min'),IQ_min=('I_Q123_Kkms','min'),IQ_max=('I_Q123_Kkms','max'),ISO_max=('I_SO2_5_2_4_Kkms','max')).reset_index(); print(g.round(3).to_string())
print('rms (mJy) of residual ~ sqrt(chi2/ndata):',np.sqrt(R.chi2.min()/R.ndata.iloc[0]))
