# Q16 (2): T21/T32/T43 at N, S, Bridge in MATCHED ELLIPTICAL beams (no circular smoothing: keeps the native N-S resolution). 
# (A) common beam = Xb4f native (J=3 and J=4 only; Xb47 4-3 is smaller than the Xb4f beam)  (B) common beam = Xb55 native 8.6x6.0" (J=2,3,4).
# Position read = bilinear value of the moment-0 (|v-vsys|<350 km/s, line-free baseline subtracted) map, Tb ratio = (S_J/S_3) * (nu_3/nu_J)^2 (RJ, same beam).
# Contamination weight of the other nucleus = Gaussian response of the matched beam at the N-S offset (5.1", position angle from the nuclear positions).
# usage: python3 -I s72_q16_native_beams.py -> out/q16/s72.txt/json
import sys,json,numpy as np; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
from s30_matchedJ_ratios import *
from regions import REGIONS
from scipy.ndimage import map_coordinates
OUT='/workspace/out/q16/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
JOBS={'J2':('Xb55',2),'J3':('Xb4f',3),'J4':('Xb47',4)}
raw={}
for k,(ous,J) in JOBS.items():
    for mol,spw in zip(('HCN','HCO+'),S.SPW[(ous,J)]):
        raw[(k,mol)]=S.load(ous,spw,S.F0[mol][J])
def pa_of(b): return b[2]
P('# native beams (maj,min,PA deg): '+'; '.join(f'{k} {raw[(k,"HCN")][6][0]:.2f}x{raw[(k,"HCN")][6][1]:.2f} PA {raw[(k,"HCN")][6][2]:.0f}' for k in JOBS))
def moment(k,mol,tgt):
    cube,v,ps,wc,x0,y0,b0=raw[(k,mol)]
    cc=S.conv(cube,ps,b0,tgt) if tgt is not None else cube
    c,core,dv,ref=S.products(cc,v,350); m=abs(v)<350
    m0=c[m].sum(0)*dv; sig=c[ref].std(0)*np.sqrt(m.sum())*dv
    return m0,sig,ps,wc,x0,y0
def read(m0,sig,ps,wc,x0,y0,ra,dec):
    px,py=wc.world_to_pixel_values(ra,dec); px-=x0; py-=y0
    return float(map_coordinates(m0,[[py],[px]],order=1)[0]),float(map_coordinates(sig,[[py],[px]],order=1)[0])
nuS=lambda mol,J:S.F0[mol][J]/1e9
res={}
def run(label,tgtfun,keys):
    P(f'\n# common beam: {label}')
    for mol in ('HCN','HCO+'):
        out={}
        for p in ('N','S','Bridge'):
            ra,dec=REGIONS[p][0],REGIONS[p][1]; vals={}
            for k in keys:
                tgt=tgtfun(k); m0,sig,ps,wc,x0,y0=moment(k,mol,tgt); vals[k]=read(m0,sig,ps,wc,x0,y0,ra,dec)
            out[p]=vals
        for p,v in out.items():
            Jn={'J2':2,'J3':3,'J4':4}; tb={k:v[k][0]*(1/nuS(mol,Jn[k])**2) for k in v}; te={k:v[k][1]/nuS(mol,Jn[k])**2 for k in v}
            s=f'  {mol:5s} {p:6s} '+' '.join(f'S_{k}={v[k][0]:.3f}+-{v[k][1]:.3f}' for k in v)
            if 'J4' in v and 'J3' in v: s+=f'  T43/T32={tb["J4"]/tb["J3"]:.3f}+-{abs(tb["J4"]/tb["J3"])*np.hypot(te["J4"]/tb["J4"],te["J3"]/tb["J3"]):.3f}'
            if 'J2' in v: s+=f'  T32/T21={tb["J3"]/tb["J2"]:.3f}+-{abs(tb["J3"]/tb["J2"])*np.hypot(te["J2"]/tb["J2"],te["J3"]/tb["J3"]):.3f}'
            P(s); res[f'{label}|{mol}|{p}']=dict(tb=tb,te=te)
b4f=raw[('J3','HCN')][6]; b55=raw[('J2','HCN')][6]
run('Xb4f native %.2fx%.2f PA%.0f (J3,J4)'%b4f,lambda k: None if k=='J3' else tuple(b4f),('J3','J4'))
run('Xb55 native %.2fx%.2f PA%.0f (J2,J3,J4)'%b55,lambda k: None if k=='J2' else tuple(b55),('J2','J3','J4'))
# beam response at the other nucleus
cd=np.cos(np.radians(43.9)); dx=(REGIONS['N'][0]-REGIONS['S'][0])*cd*3600; dy=(REGIONS['N'][1]-REGIONS['S'][1])*3600
def resp(b):
    th=np.radians(b[2]); u=np.array([dx,dy]); e_maj=np.array([np.sin(th),np.cos(th)]); e_min=np.array([np.cos(th),-np.sin(th)])
    return float(np.exp(-4*np.log(2)*((u@e_maj/b[0])**2+(u@e_min/b[1])**2)))
P(f'\n# N-S offset ({dx:.2f}",{dy:.2f}"); beam response at the other nucleus: Xb4f native {resp(b4f):.3f}; Xb55 native {resp(b55):.3f}; circular 6.5" {np.exp(-4*np.log(2)*(np.hypot(dx,dy)/6.5)**2):.3f}')
json.dump(res,open(OUT+'s72.json','w'),indent=1); open(OUT+'s72.txt','w').write('\n'.join(log))
