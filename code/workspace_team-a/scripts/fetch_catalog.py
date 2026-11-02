# Query Splatalogue (CDMS+JPL) for every spw sky-frequency range, converted to rest at z=0.0105 margin (v up to ~3300 km/s)
import warnings; warnings.filterwarnings('ignore')
import numpy as np, astropy.units as u
from astroquery.splatalogue import Splatalogue
ranges={'B5':[(163.2,166.4),(174.6,178.3)],'B7':[(338.8,342.5),(350.5,354.3)],'B6':[(248.9,252.7),(262.7,266.5)]}
zlo,zhi=0.0085,0.0102
import pandas as pd
rows=[]
for b,rs in ranges.items():
    for lo,hi in rs:
        t=Splatalogue.query_lines(lo*(1+zlo)*u.GHz,hi*(1+zhi)*u.GHz,line_lists=['CDMS','JPL'],show_upper_degeneracy=False)
        df=t.to_pandas()
        df['band']=b; rows.append(df)
        print(b,lo,hi,len(df))
df=pd.concat(rows)
df=df[['name','chemical_name','resolved_QNs','linelist','orderedfreq','upper_state_energy_K','aij','sij','sijmu2','transition_in_space','band']]
df.to_csv('/workspace/out/splatalogue_raw.csv',index=False)
