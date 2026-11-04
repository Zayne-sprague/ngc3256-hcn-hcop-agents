"""F1 helper: re-run team-b s410/s412 annulus sums with the requested 7 annuli (adds 2-3") for Xb4b 0.30" and Xb4b/Xb4d 1.05",
4 baseline/window variants + null (off-line pseudo-moment) scatter. Reuses s410/s412 functions verbatim. Output f1_annuli.json"""
import numpy as np, json
src=open('/workspace/scripts/s410_q64_imageplane_calib.py').read().split('out={}')[0]
exec(src)
ANN=[(0,.15),(.15,.25),(.25,.5),(.5,1),(1,1.5),(1.5,2),(2,3)]
def mom_null(h,v,c,hw=350):
    d=h.data[0]; w=abs(v-2780-c)<=hw
    o=(abs(v-2780)>600)&(abs(v-2780)<1500)&~(abs(v-2780-c)<=hw+60)
    b=np.nanmedian(d[o],axis=0); return np.nansum(d[w]-b,axis=0)*abs(np.diff(v).mean())
out={}
for arr,tg in [('Xb4b',[0.3,1.05]),('Xb4d',[1.05])]:
    cubes={s:load(arr,s) for s in REST}
    for s,name in [('spw27','HCN'),('spw25','HCO+')]:
        H=cubes[s][1]
        for t in tg:
            npb=Beam(t*u.arcsec).sr.to(u.arcsec**2).value/(abs(H['CDELT1'])*3600)**2
            for base,hw in VAR:
                m,b=conv(mom(cubes[s][0],cubes[s][2],base,hw),H,t); a=ann_sums(m,H)
                out[f'{arr}|{base}{hw}|{t}|{name}']={k:(np.array(a[k])/npb).tolist() for k in 'NS'}
            L=[]
            for c in (-1000,-900,900,1000,1100):
                m,b=conv(mom_null(cubes[s][0],cubes[s][2],c),H,t); a=ann_sums(m,H); L.append({k:np.array(a[k])/npb for k in 'NS'})
            out[f'{arr}|null|{t}|{name}']={k:np.std([x[k] for x in L],axis=0,ddof=1).tolist() for k in 'NS'}
    print(arr,flush=True)
json.dump(out,open('/workspace/paper_answers/r01/f1_annuli.json','w'))
