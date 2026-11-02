# Q40e: systematics of the matched-uv (u>8, 10 kl; common 8.64x6.0 beam) decomposition + the final single table (reviewer fix 1,5).
# Variants: nuclear size 3", clump size 4"/8", shifts of E1/W by +-1", +Bridge component, RMAX 11/17"; noise from filtered null windows (baseline variant);
# 10% flux calibration per band is added only to ratios across different bands (double ratios), not to same-J contrasts.
import sys; sys.path.insert(0,'/workspace/scripts')
from q40d_uvfilter import *
import q40d_uvfilter as D
import json
NAMES0=('N','S','E1','E2','W')
VARS={'base':{}, 'nuc3':dict(sizes={'N':3,'S':3}), 'clump4':dict(sizes={'E1':4,'E2':4,'W':4}), 'clump8':dict(sizes={'E1':8,'E2':8,'W':8}),
      'bridge':dict(names=('N','S','Bridge','E1','E2','W')), 'rmax11':dict(rmax=11.), 'rmax17':dict(rmax=17.)}
for c in ('E1','W'):
    for sx,sy in ((1,0),(-1,0),(0,1),(0,-1)): VARS[f'sh_{c}_{sx}{sy}']=dict(shift={c:(sx,sy)})
def templ2(names,umin,shift,sizes):
    out=[]
    for n in names:
        ex,ey=B.EXY[n]; sx,sy=shift.get(n,(0,0)); Si=np.linalg.inv(B.cov2(*TB,s=sizes.get(n,0.)))
        X=np.stack([XG-ex-sx,YG-ey-sy],-1); g=np.exp(-0.5*np.einsum('...i,ij,...j->...',X,Si,X))
        F=np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(g))); T=0.5*(1+erf((RK-umin)/1.5)) if umin>0 else 1.0
        out.append(np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F*T)))))
    return out
def fit2(y,cols,rmax):
    s=np.hypot(XG,YG)<rmax; A=np.array([c[s] for c in cols]).T; x,_=nnls(A,y[s]); return x
OUT={}
DS={}
for J,nm in SETS:
    for ln,lab in ((0,'HCN'),(1,'HCO+')): DS[(nm,lab)]=B.getmaps(nm,ln)
for umin in (6,8,10):
    FM={}
    for key,(d,m,ns) in DS.items():
        FM[key]=(filt(d,m,umin,TB)[0],[filt(d,np.nan_to_num(n),umin,TB)[0] for n in ns])
    for vn,v in VARS.items():
        names=v.get('names',NAMES0); cols=templ2(names,umin,v.get('shift',{}),v.get('sizes',{})); rmax=v.get('rmax',14.)
        for J,nm in SETS:
            for lab in ('HCN','HCO+'):
                y,nl=FM[(nm,lab)]; a=fit2(y,cols,rmax); e=[np.nan]*len(names)
                if vn=='base':
                    dd=[fit2(y+yn,cols,rmax)-a for yn in nl]; e=np.std(dd,axis=0)
                OUT[(umin,vn,nm,lab)]=dict(names=list(names),amp=a.tolist(),err=np.array(e).tolist())
import pickle; pickle.dump(OUT,open('/workspace/out/q40b/q40e_sys.pkl','wb'))
print('done',len(OUT))
