# Bug-free line identification: match observed sky frequency (at v_sys, tolerance +-dvtol km/s) against the full Splatalogue CDMS+JPL dump
# (all species, Eu<Emax, NO Aij cut), rank by relative LTE strength  S*mu2 * nu * exp(-Eu/T) (T = 30 K; optically thin, partition function ignored)
# JPL/CDMS duplicates are merged on (species, rounded frequency).  Unit test: python3 -I ident.py --test
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from catalog import load
from curated import SPEC
SPECSET={x.split()[0].strip(',') for x in SPEC}   # first token: names carry electronic-state/vib labels (e.g. 'SO 3Sigma v=0')
C=299792.458
_df=None
def cat(emax=250):
    global _df
    if _df is None:
        d=load(emax); d=d[~((d.Eu==0)&(d.logA==0))].copy()
        d['sp']=d.sp.str.replace(r'\s+v\s*=.*$','',regex=True).str.replace(r'\s+vt\s*=.*$','',regex=True)
        d=d[d.sijmu2>0]
        d['key']=d.sp+'_'+d.f0.round(4).astype(str)
        d=d.sort_values('linelist').drop_duplicates('key')   # CDMS first
        _df=d
    return _df
def rel_strength(d,T=30.): return d.sijmu2*d.f0*np.exp(-d.Eu/T)
def ident(fsky,v0=2780.,dvtol=100.,emax=200.,T=30.,top=6):
    d=cat(); f0lo,f0hi=sorted([fsky/(1-(v0-dvtol)/C),fsky/(1-(v0+dvtol)/C)])
    s=d[(d.f0>=f0lo)&(d.f0<=f0hi)&(d.Eu<emax)].copy()
    s['v']=(1-fsky/s.f0)*C; s['rel']=rel_strength(s,T)
    # normalise per species to the strongest in-band line of that species at the fitted T? -> keep absolute rel, then sort
    s['cur']=s.sp.str.split().str[0].str.strip(',').isin(SPECSET)   # species known/expected in LIRGs (Martin+06,11,21; Harada+18; Aladro+15); others only as fallback
    s=s.sort_values(['cur','rel'],ascending=False)
    return s.head(top)
def fmt(s): return '; '.join(f'{"" if r.cur else "[exotic] "}{r.sp} {str(r.qn).strip()[:24]} ({r.f0:.4f}, Eu={r.Eu:.0f}K, logA={r.logA:.2f}, v={r.v:+.0f}, rel={r.rel:.1f})' for r in s.itertuples())
if __name__=='__main__':
    if '--test' in sys.argv:
        ok=True
        def has(fsky,sp,v0=2780.,dv=60):
            s=ident(fsky,v0,dv,emax=250,top=40); return s.sp.str.startswith(sp).any()
        for name,f0 in [('HCO+',267.557626),('HCN',265.886434),('CS',342.883000),('SO',178.605403),('SO',251.825770),('NS',253.968393),('HOC+',268.451094),('HC3N',254.6995),('HCO+',356.734223),('HCN',354.505478),('HCN',177.261111),('HCO+',178.375056)]:
            r=has(f0*(1-2780/C),name); ok&=r; print('recover',name,f0,r)
        for J,f0 in [(1,165.050175),(2,165.061130),(3,165.099240),(4,165.190475),(5,165.369341),(2.5,254.015377)]:
            r=has(f0*(1-2780/C),'CH3OH',dv=5); ok&=r; print('recover CH3OH',f0,r)
        print('ALL OK' if ok else 'FAIL'); sys.exit(0 if ok else 1)
    pd.set_option('display.width',250,'display.max_colwidth',200)
    for f in map(float,sys.argv[1:]): print(f, fmt(ident(f)))
