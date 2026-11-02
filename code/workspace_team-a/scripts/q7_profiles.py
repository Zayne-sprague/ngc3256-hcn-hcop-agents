# Q7: HCN/HCO+ 3-2 radial profiles (annuli) about the N and S peaks, matched beams, several velocity windows, null-window noise.
# python3 -I q7_profiles.py   (needs out/q7/cube_*.npz from q7_cache.py)
import sys; sys.path.insert(0,'/workspace/scripts')
from q7_lib import *
import pandas as pd
rng=np.random.default_rng(7)
cb={(o,l):Cube(o,l) for o in ('Xb4b','Xb4d') for l in ('HCN','HCO+')}
WINS={'pm350':(-350,350),'pm165':(-165,165),'core100':(-100,100),'blue':(-350,-100),'red':(100,350)}
CASES={'Xb4b@0.30':('Xb4b',(0.30,0.30,0.)),'Xb4b@1.05':('Xb4b',(1.05,1.05,0.)),'Xb4d@1.05':('Xb4d',(1.05,1.05,0.))}
EDGES={'Xb4b@0.30':[0,0.15,0.25,0.5,1,1.5,2,3],'Xb4b@1.05':[0,0.5,1,1.5,2,3],'Xb4d@1.05':[0,0.5,1,1.5,2,3]}
NULLW=(-760,-540)   # null window (200 km/s); baseline for the null refit excludes it
def null_m0(c,w,wn=NULLW):
    m=np.zeros(len(c.v),bool)
    for a,b in BASE: m|=(c.v>=a)&(c.v<=b)
    m&=~((c.v>=wn[0]-10)&(c.v<=wn[1]+10)); raw=c.cube  # cube is already baseline-subtracted with the full fit; add back is not needed: refit residual baseline
    X=np.vstack([np.ones(m.sum()),c.v[m]/500.]).T; P=np.linalg.pinv(X); cc=np.tensordot(P,raw[m],axes=(1,0))
    k=(c.v>=wn[0])&(c.v<=wn[1]); Xk=np.vstack([np.ones(k.sum()),c.v[k]/500.]).T
    res=raw[k]-np.tensordot(Xk,cc,axes=(1,0)); return res.sum(0)*c.dv*np.sqrt((w[1]-w[0])/(wn[1]-wn[0]))
def centres(c,tb,x0,y0,which):
    img=smooth(c.m0(-350,350),c,tb) if tb else c.m0(-350,350)
    return peak_near(c,img,x0,y0)
REG={'N':(0.0,2.55),'S':(-0.15,-2.55)}
rows=[]
for case,(ous,tb) in CASES.items():
    H,O=cb[(ous,'HCN')],cb[(ous,'HCO+')]
    # centres: peak of the HCO+ and HCN +-350 maps at this beam; adopted = HCO+ peak (variants tested at 0.30)
    cen={}
    for reg,(x0,y0) in REG.items():
        cen[reg]={'HCO+':peak_near(O,smooth(O.m0(-350,350),O,tb),x0,y0),'HCN':peak_near(H,smooth(H.m0(-350,350),H,tb),x0,y0)}
    for wn,w in WINS.items():
        mH,mO=smooth(H.m0(*w),H,tb),smooth(O.m0(*w),O,tb)
        nH=smooth(null_m0(H,w),H,tb); nO=smooth(null_m0(O,w),O,tb)
        # noise per annulus: annulus sums at random centres in the null maps
        edges=EDGES[case]; sig={}
        pts=rng.uniform(-4.5,4.5,(300,2))
        for a,b in zip(edges[:-1],edges[1:]):
            sig[(a,b)]=tuple(1.4826*np.median(np.abs(np.array([annulus_flux(H,n,x,y,a,b,tb)[0] for x,y in pts])-np.median([annulus_flux(H,n,x,y,a,b,tb)[0] for x,y in pts]))) for n in (nH,nO))
        for reg in REG:
            for cl in ('HCO+','HCN'):
                if cl=='HCN' and (wn!='pm350' or case!='Xb4b@0.30'): continue
                x0,y0=cen[reg][cl]
                for a,b in zip(edges[:-1],edges[1:]):
                    fH,n=annulus_flux(H,mH,x0,y0,a,b,tb); fO,_=annulus_flux(O,mO,x0,y0,a,b,tb); eH,eO=sig[(a,b)]
                    # cumulative aperture too
                    cH,_=annulus_flux(H,mH,x0,y0,0,b,tb); cO,_=annulus_flux(O,mO,x0,y0,0,b,tb)
                    rows.append(dict(case=case,win=wn,reg=reg,centre=cl,x0=round(x0,3),y0=round(y0,3),r1=a,r2=b,npix=n,HCN=fH,eHCN=eH,HCOp=fO,eHCOp=eO,R=fH/fO,eR=abs(fH/fO)*np.hypot(eH/max(abs(fH),1e-9),eO/fO),Rcum=cH/cO))
df=pd.DataFrame(rows); df.to_csv('/workspace/out/q7/q7_profiles.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500)
print(df[df.centre=='HCO+'].round(3).to_string())
print(df[df.centre=='HCN'].round(3).to_string())
