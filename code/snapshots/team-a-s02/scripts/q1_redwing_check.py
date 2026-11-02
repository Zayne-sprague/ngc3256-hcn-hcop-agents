# Q1: inspect HCN vs HCO+ 3-2 profiles at N and S in the red wing (+100..+300 km/s): is the high ratio a real line-shape difference, a neighbouring line (HCN side), or absorption?
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
exec(open('/workspace/scripts/q1_velratio.py').read().split("SP=spectra")[0])
import pandas as pd
rows=[]
for nm in ('N','S','Bridge'):
    S=spectra('Xb4d','spw27','spw25',{nm:POS[nm]})[nm]; (vH,sH),(vO,sO)=S['HCN'],S['HCO+']
    eH,eO=noise(vH,sH),noise(vO,sO)
    print(nm,'noise HCN/HCO+ (mJy/beam/chan):',eH*1e3,eO*1e3)
    for lo in range(-100,420,25):
        mH=(vH-VSYS>=lo)&(vH-VSYS<lo+25); mO=(vO-VSYS>=lo)&(vO-VSYS<lo+25)
        a=sH[mH].mean()*1e3; b=sO[mO].mean()*1e3; rows.append(dict(region=nm,v=lo+12.5,HCN_mJy=a,HCO_mJy=b,snrH=a/(eH*1e3/np.sqrt(mH.sum())),snrO=b/(eO*1e3/np.sqrt(mO.sum()))))
R=pd.DataFrame(rows); print(R.round(2).to_string()); R.to_csv('/workspace/out/q1/q1_redwing_profiles.csv',index=False)
# which known lines fall at +100..+250 km/s of HCN (sky GHz) ?
for dv in (100,150,200,250):
    print(dv,'km/s above HCN systemic -> sky freq %.4f GHz'%(FH*(1-(VSYS+dv)/C)))
