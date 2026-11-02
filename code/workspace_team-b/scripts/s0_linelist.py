# Line list for the rest-frame ranges covered: Splatalogue (CDMS+JPL), curated species relevant to extragalactic mm spectra,
# plus hydrogen recombination lines computed from the Rydberg formula.
import astropy.units as u, pandas as pd, numpy as np, sys, re
from astroquery.splatalogue import Splatalogue
rng=[(164.4,179.8),(251.2,268.1),(341.9,357.3)]
KEEP=r'^(CO|13CO|C18O|C17O|13C18O|HCN|H13CN|HC15N|HCO\+|H13CO\+|HC18O\+|HC17O\+|HNC|HN13C|H15NC|CS|13CS|C34S|C33S|SO|34SO|33SO|S18O|SO\+|SiO|29SiO|30SiO|SiS|CN|13CN|C15N|CH3OH|13CH3OH|CH3CN|CH313CN|13CH3CN|HC3N|H13CCCN|HC13CCN|HCC13CN|CCH|C2H|13CCH|C13CH|H2CO|H213CO|H2CS|H2S|H234S|c-C3H2|c-HCCCH|c-C3H|l-C3H|CH3CCH|HNCO|SO2|34SO2|OCS|O13CS|NS|NO|HOC\+|N2H\+|HCS\+|HCO|CH3CHO|CH3OCH3|HCOOCH3|NH2CHO|H2CCO|HDO|C4H|CCS|HC5N|SiC2|HOCO\+|H2COH\+|NH2CN|NH2D|HNO|OCS)$'
out=[]
for a,b in rng:
    for lo in np.arange(a,b,1.0):
        hi=min(lo+1.0,b)
        t=Splatalogue.query_lines(lo*u.GHz,hi*u.GHz,line_lists=['CDMS','JPL'])
        df=t.to_pandas(); df['name']=df.name.str.replace(r'<[^>]+>','',regex=True).str.strip()
        df['sp']=df.name.str.split().str[0]; df['vib']=df.name.str.contains('v|v2|v1',regex=True)
        df=df[df.sp.str.match(KEEP)]; out.append(df); print(lo,len(t),len(df),file=sys.stderr)
df=pd.concat(out).drop_duplicates('lineid')
df['freq_GHz']=df.orderedfreq/1e3
df=df[['name','sp','resolved_QNs','linelist','freq_GHz','upper_state_energy_K','aij','sij','sijmu2']]
# H recombination lines (n+1 -> n) Rydberg formula with reduced electron mass
RH=3.2898419603e15*(1/(1+5.4461702e-4))   # Hz
rows=[]
for n in range(20,60):
    for dn in (1,2):
        f=RH*(1/n**2-1/(n+dn)**2)/1e9
        if any(a<=f<=b for a,b in rng): rows.append(dict(name=f'H{n}{"α" if dn==1 else "β"} recomb.',sp='H-RRL',resolved_QNs='',linelist='calc',freq_GHz=f,upper_state_energy_K=0,aij=np.nan,sij=np.nan,sijmu2=np.nan))
df=pd.concat([df,pd.DataFrame(rows)]).sort_values('freq_GHz')
df.to_csv('/workspace/out/linelist.csv',index=False); print(len(df)); print(pd.DataFrame(rows)[['name','freq_GHz']])
