import sys; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
import numpy as np
# independent: native Xb4d cube HCO+ 3-2, beam-area-rescaled flux in r=1" aperture, field-subtract w/ ring offsets
for mol in ('HCO+','HCN'):
    R=load_any('std_Xb4d',mol,1.05,12.0); v=R['v']; dv=R['dv']
    m=apmask(R['n'],R['ps'],R['pxy'](*LOBE),1.0); sp=aspec(R['c'],R['ba'],m)
    # off-source apertures in field
    fm=field_pixels(R); ys,xs=np.where(fm); Y,X=np.mgrid[:R['n'],:R['n']]
    rng=np.random.default_rng(5); fs=[]
    while len(fs)<80:
        i=rng.integers(len(ys)); mm=np.hypot(X-xs[i],Y-ys[i])*R['ps']<=1.0
        if fm[mm].all(): fs.append(aspec(R['c'],R['ba'],mm))
    fs=np.array(fs); sp2=sp-fs.mean(0)
    w=(v>=-405)&(v<=-150)
    print(mol,'raw',sp[w].sum()*dv,'fieldsub',sp2[w].sum()*dv,'field sigma',np.std([f[w].sum()*dv for f in fs]), 'dv',dv)
    # continuum-like baseline from far channels
    far=(v<-700)|(v>500); print(' far-chan mean flux/chan',sp[far].mean(), 'rms',sp[far].std())
