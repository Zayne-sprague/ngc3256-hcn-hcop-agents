"""q01: dilution test. Xb4b (0.25") HCN/HCO+ 3-2 moment-0 (med350 baseline, +-350 km/s about 2780, as team-b s410/Table C1)
convolved to circular beams; N/S = (HCN/HCO+)_N/(HCN/HCO+)_S for
 (a) Table C1 estimator: annulus sums about the 0.24" continuum peaks, r<0.25" and r<1"
 (b) uvtaper (s192) aperture: mean in r = beam FWHM aperture centred on HCO+ moment peak within 0.4" of REGIONS N,S.
Null scatter: 5 off-line pseudo-moments (+-900,1000,1100 km/s) -> stat error. Reuses s410 functions verbatim. Output q01_values.json"""
import numpy as np, json, sys
sys.path.insert(0,'/workspace/scripts')
from regions import REGIONS
src=open('/workspace/scripts/s410_q64_imageplane_calib.py').read().split('out={}')[0]
exec(src)
BEAMS=[None,0.3,0.6,1.05,1.5,2.2]
def mom_null(h,v,c,hw=350):
    d=h.data[0]; w=abs(v-2780-c)<=hw
    o=(abs(v-2780)>600)&(abs(v-2780)<1500)&~(abs(v-2780-c)<=hw+60)
    b=np.nanmedian(d[o],axis=0); return np.nansum(d[w]-b,axis=0)*abs(np.diff(v).mean())
def apsum(m,H,hco,rad):
    w=WCS(H).celestial; y,x=np.mgrid[:m.shape[0],:m.shape[1]]; pix=abs(H['CDELT1'])*3600; r={}
    for k in 'NS':
        px,py=[float(t) for t in w.world_to_pixel_values(*REGIONS[k][:2])]
        mm=np.where(np.hypot(x-px,y-py)*pix<0.4,hco,-np.inf); j=np.argmax(mm); cx,cy=x.flat[j],y.flat[j]
        r[k]=float(np.mean(m[np.hypot(x-cx,y-cy)*pix<=rad]))
    return r
cubes={s:load('Xb4b',s) for s in REST}
names={'spw27':'HCN','spw25':'HCO+'}
def measure(M,t,peakmap=None):
    """M: dict spw->moment; returns dict est->{line:{N,S}}"""
    C={}
    for s in REST:
        m,b=conv(M[s],cubes[s][1],t); C[s]=m
    H=cubes['spw27'][1]; bm=t if t else 0.246
    out={'ann025':{},'ann1':{},'ap':{}}
    for s in REST:
        a=ann_sums(C[s],H); n=names[s]
        out['ann025'][n]={k:sum(a[k][:2]) for k in 'NS'}; out['ann1'][n]={k:sum(a[k][:4]) for k in 'NS'}
    pk=peakmap if peakmap is not None else C['spw25']
    for s in REST: out['ap'][names[s]]=apsum(C[s],H,pk,bm)
    return out
def ns(o): return (o['HCN']['N']/o['HCO+']['N'])/(o['HCN']['S']/o['HCO+']['S'])
res={}
M={s:mom(cubes[s][0],cubes[s][2],'med',350) for s in REST}
NM={c:{s:mom_null(cubes[s][0],cubes[s][2],c) for s in REST} for c in (-1000,-900,900,1000,1100)}
for t in BEAMS:
    o=measure(M,t); e={}
    nl=[]
    for c in NM:
        # null: pseudo-line maps; aperture centre taken from the real HCO+ map (peak position fixed)
        pk=conv(M['spw25'],cubes['spw25'][1],t)[0]
        nl.append(measure(NM[c],t,peakmap=pk))
    for est in o:
        # propagate null scatter of each of the 4 integrals onto N/S
        sig={(l,k):np.std([x[est][l][k] for x in nl],ddof=1) for l in ('HCN','HCO+') for k in 'NS'}
        r=ns(o[est]); e[est]=r*np.sqrt(sum((sig[(l,k)]/o[est][l][k])**2 for l in ('HCN','HCO+') for k in 'NS'))
    res[str(t)]={est:dict(NS=ns(o[est]),err=e[est],HCN_HCOp_N=o[est]['HCN']['N']/o[est]['HCO+']['N'],HCN_HCOp_S=o[est]['HCN']['S']/o[est]['HCO+']['S']) for est in o}
    print(t,{k:(round(v['NS'],3),round(v['err'],3)) for k,v in res[str(t)].items()},flush=True)
json.dump(res,open('/workspace/paper_answers/q01/q01_values.json','w'),indent=1)
