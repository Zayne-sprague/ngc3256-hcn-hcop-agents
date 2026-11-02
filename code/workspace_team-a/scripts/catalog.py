import pandas as pd, numpy as np, re
C=299792.458
def load(emax=400):
    df=pd.read_csv('/workspace/out/splatalogue_raw.csv')
    df['sp']=df.name.str.replace(r'<[^>]+>','',regex=True).str.replace(r'\s+',' ',regex=True).str.strip()
    df['f0']=df.orderedfreq/1e3  # GHz
    df=df.rename(columns={'upper_state_energy_K':'Eu','aij':'logA'})
    df['qn']=df.resolved_QNs.astype(str).str.strip()
    df=df.drop_duplicates(['sp','qn','f0'])
    return df[df.Eu<emax]
def cands(df,fsky,vlo=2500,vhi=3050,minlogA=-5.5):
    # radio convention: f=f0(1-v/c)
    f0lo=fsky/(1-vlo/C); f0hi=fsky/(1-vhi/C)
    s=df[(df.f0>=min(f0lo,f0hi))&(df.f0<=max(f0lo,f0hi))&(df.logA>minlogA)].copy()
    s['v']=(1-fsky/s.f0)*C
    return s.sort_values('logA',ascending=False)
if __name__=='__main__':
    import sys
    df=load()
    for fsky in [265.980,352.394,177.311,176.958,176.932,163.539,164.053,251.641,249.505,342.1,338.9]:
        s=cands(df,fsky,2450,3100)
        s=s[s.Eu<300]
        print('\n== fsky',fsky,'GHz'); print(s[['sp','qn','f0','Eu','logA','v','linelist']].head(12).to_string(index=False))
