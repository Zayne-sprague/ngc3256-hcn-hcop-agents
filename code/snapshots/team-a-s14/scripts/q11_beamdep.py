# Q11c: beam-shape dependence (core spill into the 8-12" annulus) of the annulus HCN/HCO+ ratios; far-field zero level (30-45", Xb4f/Xb55) and none.
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
ANN=[(0,4),(4,8),(8,12)]
def sm(d,tb):
    k=kernel(d['hd'],tb,d['beam']); 
    return None if k is None else fftconvolve(np.nan_to_num(d['m0']),k,mode='same')*OM(tb)/OM(d['beam'])
def ann(d,s,fl):
    r=np.hypot(d['dx']+0.075,d['dy']); omp=OM(tb_)/d['pix']**2
    z=np.nanmean(s[(r>=fl[0])&(r<fl[1])]) if fl else 0
    return np.array([(s[(r>=a)&(r<b)]-z).sum()/omp for a,b in ANN])
sets={'2-1 Xb55':('HCN 2-1','HCO+ 2-1'),'3-2 Xb4f':('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f'),'3-2 Xb4d':('HCN 3-2@Xb4d','HCO+ 3-2@Xb4d')}
beams={'7.0 circ':(7.,7.,0.),'8.5 circ':(8.5,8.5,0.),'10 circ':(10.,10.,0.),'12 circ':(12.,12.,0.),'9.36x6.28 (Xb55)':(9.36,6.28,-88.),'9.36x6.28 rot 90':(9.36,6.28,2.)}
for s,(kh,ko) in sets.items():
    H=load(kh,nullspec=False); O=load(ko,nullspec=False)
    for bn,tb in beams.items():
        tb_=tb; a=sm(H,tb); b=sm(O,tb)
        if a is None or b is None: continue
        fl=(30,45) if s!='3-2 Xb4d' else (20,32)
        for lab,f_ in (('none',None),('ff',fl)):
            h,o=ann(H,a,f_),ann(O,b,f_)
            print(f'{s:9s} beam {bn:18s} zero {lab:4s} HCN {np.round(h,2)} HCO+ {np.round(o,2)} R {np.round(h/o,3)}')
