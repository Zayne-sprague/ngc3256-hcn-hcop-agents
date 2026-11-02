# Rest-frame coverage (at VSYS=2780 km/s radio) of every cube spw and which of the proposal's / literature's lines fall inside (full profile +-300 km/s)
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd
T={'HCN 2-1':177.26111,'HCO+ 2-1':178.37506,'HCN 3-2':265.88618,'HCO+ 3-2':267.55763,'HCN 4-3':354.50547,'HCO+ 4-3':356.73424,'H13CN 4-3':345.33977,'H13CO+ 4-3':346.99835,
   'CS 7-6':342.88285,'SO 6_5-5_4':219.94944,'SO 8_8-7_7':344.31061,'SO 5_6-4_5':251.82577,'SO 6_6-5_5':258.25583,'SiO 6-5':260.51802,'SiO 8-7':347.33058,'CO 3-2':345.79599,'13CO 3-2':330.58797,
   'HNC 4-3':362.63030,'CN 3-2 (J=5/2-3/2)':340.2478,'CH3OH 2-1 E ladder (165.05)':165.06113,'H26a':353.62277,'H29a':256.30205,'H13CO+ 3-2':260.25534,'HN13C 3-2':261.26331,'CCH 3-2':262.00426,'C34S 7-6':337.39646,'HC3N 28-27':254.6997,'NS 11/2-9/2 f':253.96793}
rows=[]; seen=set()
for fn in files('*.cube.I.pbcor.fits')+files('*tapered*cube*'):
    h=fits.getheader(fn); n=h['NAXIS3']; nu=h['CRVAL3']+(np.arange(n)+1-h['CRPIX3'])*h['CDELT3']
    lo,hi=rest_of(nu.min())/1e9,rest_of(nu.max())/1e9
    key=(short(fn).split('.')[0],short(fn).split('.')[1]) if 'tapered' not in fn else None
    inside=[]
    for k,f in T.items():
        half=300/C*f
        if f-half>=lo and f+half<=hi: inside.append(k)
        elif f>=lo-half and f<=hi+half: inside.append(k+' (partial/edge)')
    rows.append(dict(file=short(fn),obs_GHz=f'{nu.min()/1e9:.3f}-{nu.max()/1e9:.3f}',rest_GHz_at_VSYS=f'{lo:.3f}-{hi:.3f}',dchan_kms=abs(C*np.median(np.diff(nu))/nu.mean()),lines_covered='; '.join(inside)))
df=pd.DataFrame(rows).drop_duplicates('file'); df.to_csv('/workspace/out/s10_coverage.csv',index=False)
pd.set_option('display.width',250,'display.max_colwidth',120,'display.max_rows',100); print(df.to_string())
