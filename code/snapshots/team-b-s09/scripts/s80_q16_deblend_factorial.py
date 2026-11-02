# Q16r2 (5): which analysis choice moves the deblended S core ladder from ~0.3 (team-a p0012: 0.268) to ~0.44?  Full factorial on the OBSERVED maps with the s79 NNLS template deblend:
#  template: clipped (>=0, noise included) vs masked (>4 sigma);  baseline: none vs local plane (ring 12-22");  map error: formal (0.25) vs empirical ring rms (0.6/1.0);
#  fit radius 10/14";  core radius 1.0/1.5/2.0";  Outer component on/off.   usage: python3 -I s80_q16_deblend_factorial.py -> out/q16b/s80.txt/json
import sys,json,itertools; sys.path.insert(0,'/workspace/scripts')
from s79_q16_injection import *
OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
res={}
for mol in ('HCN','HCO+'):
    P(f'# {mol} (template = Xb4d {mol} 3-2). team-a: '+('N 0.395 S 0.268' if mol=='HCN' else 'N 0.482 S 0.442'))
    for tmpl,ring in itertools.product(('clip','mask'),(None,(12,22))):
        D=setup(mol,mol,thr=None if tmpl=='clip' else 4.0,ring=ring)
        for errm,Rfit,rc,outer in itertools.product(('formal','emp'),(10,14),(1.0,1.5,2.0),(True,False)):
            nz=0.25 if errm=='formal' else (0.6 if mol=='HCN' else 1.0)
            r=deblend(D,D['d3'],D['d4'],rc=rc,Rfit=Rfit,noise=nz,outer=outer)
            res[f'{mol}|{tmpl}|{ring}|{errm}|{Rfit}|{rc}|{outer}']={k:float(v) for k,v in r.items()}
        # summary lines per (tmpl,ring): median over the other factors
        sub=[v for k,v in res.items() if k.startswith(f'{mol}|{tmpl}|{ring}|')]
        for cmp_ in ('N','S','Bridge'):
            pass
        for errm in ('formal','emp'):
            sb=[v for k,v in res.items() if k.startswith(f'{mol}|{tmpl}|{ring}|{errm}|')]
            q=lambda c: np.nanpercentile([x.get(c,np.nan) for x in sb],[0,50,100])
            P(f'  template {tmpl:4s} baseline {str(ring):9s} errors {errm:6s}: N {q("N")[1]:.3f} [{q("N")[0]:.3f},{q("N")[2]:.3f}]  S {q("S")[1]:.3f} [{q("S")[0]:.3f},{q("S")[2]:.3f}]  Bridge {q("Bridge")[1]:.3f} [{q("Bridge")[0]:.3f},{q("Bridge")[2]:.3f}]   (median [min,max] over Rfit, rc, outer)')
json.dump(res,open(OUT+'s80.json','w'),indent=1); open(OUT+'s80.txt','w').write('\n'.join(log))
