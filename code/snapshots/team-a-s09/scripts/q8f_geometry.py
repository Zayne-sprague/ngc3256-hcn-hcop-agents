# Q8f: per-line source sizes, beam filling and compact share at N,S,Bridge from the joint uv-forward-model fits (q8c_uvmodel.py), unfiltered (intrinsic) geometry.
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
rows=[]; names=None
BEAMS={'Xb4d 0.95x0.75':(0.95,0.75),'Xb47 4.2x2.75':(4.2,2.75),'Xb4f 5.9x3.6':(5.9,3.6),'Xb55 9x6':(9.0,6.0)}
for V in ('taper_phys','taper_raw_phys','taper','taper_raw'):
    Z=np.load(f'/workspace/out/q8c/chain_{V}.npz',allow_pickle=True); fl=Z['flat']; nm=list(Z['names']); ix={n:i for i,n in enumerate(nm)}
    th=fl[np.random.default_rng(3).choice(len(fl),6000,replace=False)]
    for reg in ('N','S','Bridge'):
        for l in ('HCN','HCO+'):
            s,se,Ac,Ae=[np.exp(th[:,ix[f'{reg}_{l}_{k}']]) for k in ('lns','lnse','lnAc','lnAe')]
            q=lambda x:np.percentile(x,[16,50,84])
            rows.append(dict(family=V,region=reg,line=l,beam='size FWHM compact (")',p16=q(s)[0],p50=q(s)[1],p84=q(s)[2]))
            rows.append(dict(family=V,region=reg,line=l,beam='size FWHM extended (")',p16=q(se)[0],p50=q(se)[1],p84=q(se)[2]))
            rows.append(dict(family=V,region=reg,line=l,beam='intrinsic peak Tb compact (K)',p16=q(Ac)[0],p50=q(Ac)[1],p84=q(Ac)[2]))
            rows.append(dict(family=V,region=reg,line=l,beam='intrinsic peak Tb extended (K)',p16=q(Ae)[0],p50=q(Ae)[1],p84=q(Ae)[2]))
            for bn,b in BEAMS.items():
                fc=s**2/np.sqrt((s**2+b[0]**2)*(s**2+b[1]**2)); fe=se**2/np.sqrt((se**2+b[0]**2)*(se**2+b[1]**2)); sh=Ac*fc/(Ac*fc+Ae*fe)
                rows.append(dict(family=V,region=reg,line=l,beam=f'filling compact, {bn}',p16=q(fc)[0],p50=q(fc)[1],p84=q(fc)[2]))
                rows.append(dict(family=V,region=reg,line=l,beam=f'filling extended (u), {bn}',p16=q(fe)[0],p50=q(fe)[1],p84=q(fe)[2]))
                rows.append(dict(family=V,region=reg,line=l,beam=f'compact share of unfiltered peak, {bn}',p16=q(sh)[0],p50=q(sh)[1],p84=q(sh)[2]))
G=pd.DataFrame(rows); G.to_csv('/workspace/out/q8c/geometry_uvmodel.csv',index=False); pd.set_option('display.width',200,'display.max_rows',900)
print(G[G.beam.str.contains('size|Xb4f|Xb47')&G.family.isin(['taper_phys','taper'])&(G.line=='HCN')].round(3).to_string())
