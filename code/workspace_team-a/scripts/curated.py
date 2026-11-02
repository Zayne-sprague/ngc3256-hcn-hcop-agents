# Curated catalogue: species plausibly detectable in LIRGs (Martin+2006,2011,2021; Harada+2018; Aladro+2015), Eu<250 K, from Splatalogue (CDMS/JPL) raw dump
import pandas as pd, numpy as np, re, sys
sys.path.insert(0,'/workspace/scripts')
from catalog import load
SPEC=['CO','13CO','C18O','C17O','HCN','H13CN','HC15N','HCO+','H13CO+','HC18O+','HOC+','HNC','HN13C','CS','C34S','13CS','C33S','SO','34SO','SO2','SiO','CH3OH','13CH3OH','HC3N','CH3CCH','CH3CN','H2CS','H2CO','H2S','OCS','NS','NO','CN','13CN','CCH','N2H+','HNCO','HCS+','CH3CHO','c-C3H2','H2CCO','NH2CHO','HC3N, v7=1','HC3N v7=1','HC3N v6=1','HC3N v5=1/v7=3','HC3N v7=2','HCN v2=1','HCN, v2=1','HCO+ v2=1','HNC v2=1','CH3SH','C2S','CCS','H2C18O','H2C17O','HDO','H2O','CH2NH','NH3','HCCCN','CCD','SiS','PN','PO','HCP','HNCS','CH3NH2','HC5N','l-C3H','HOCO+','H2COH+','SO+','H3O+','H13CO+','DCN','DCO+','NH2D','CH3OCH3','HCOOCH3','CH3CH2OH','C2H5CN','HCOOH','CH2CO','N2D+','NH2CN','OH','OCS','CH3CHO','CH3CN v8=1']
def curated(emax=250):
    df=load(emax)
    df=df[~((df.Eu==0)&(df.logA==0))]
    keep=df.sp.apply(lambda s: any(s.startswith(x+' ') or s==x or s.startswith(x+',') for x in SPEC))
    df=df[keep]
    return df
if __name__=='__main__':
    df=curated(); print(len(df)); print(df.sp.value_counts().head(60))
