# Q3 round 2, fix 2: primary-beam gain per band at each position, S/N per line, and the effect of PB correction on the 10" fluxes.
# The cubes are pipeline pbcor products (no .pb files delivered) -> the PB response is inferred from the per-pixel noise of the pbcor maps (sigma_pbcor ~ sigma_raw/PB):
#   ln sigma(x,y) = a + 4 ln2 [(x-x0)^2+(y-y0)^2]/FWHM_eff^2 fitted to the 10" error maps (out/s30_maps_*_10as.npz, HCN and HCO+).
# Then PB(position) = exp(-4 ln2 r^2 / FWHM_eff^2) relative to the fitted centre. Test of PB handling: smooth(pbcor) vs smooth(pbcor*PB)/PB at the position (native maps -> 10").
# usage: python3 -I s55a_pb_snr.py -> out/q3b/s55a_pb_snr.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
from scipy.optimize import least_squares
from scipy.ndimage import gaussian_filter,map_coordinates
from regions import REGIONS
OUT='/workspace/out/q3b/'; log=[]
def P(s=''): print(s); log.append(s)
CEN=[0.5*(156.96346+156.96340),0.5*(-43.90388-43.90533)]
def xy_of(name,d):   # pixel position in the cropped map (RA increases to the left)
    ra,dec=REGIONS[name][:2] if name in REGIONS else CEN
    dx=(ra-CEN[0])*np.cos(np.deg2rad(CEN[1]))*3600/float(d['ps']); dy=(dec-CEN[1])*3600/float(d['ps'])
    return float(d['cx'])-dx,float(d['cy'])+dy
BANDS={'Xb55 (B5, J=2)':'Xb55_J2','Xb4f (B6, J=3)':'Xb4f_J3','Xb47 (B7, J=4)':'Xb47_J4'}
pb={}; res={'pb':{},'snr':{},'pbtest':{}}
for lab,tg in BANDS.items():
    d=np.load(f'/workspace/out/s30_maps_{tg}_10as.npz'); ps=float(d['ps']); n=d['eHCN'].shape[0]; Y,X=np.mgrid[:n,:n]
    out=[]
    for key in ('eHCN','eHCOp'):
        e=d[key]; ok=np.isfinite(e)&(e>0)&(d['rr']<22)
        def f(p): a,x0,y0,b=p; return (a+b*((X-x0)**2+(Y-y0)**2)*ps**2)[ok]-np.log(e[ok])
        s=least_squares(f,[np.log(np.median(e[ok])),n/2,n/2,1e-3]); a,x0,y0,b=s.x
        fw=np.sqrt(4*np.log(2)/b) if b>0 else np.nan; out.append((a,x0,y0,b,fw,np.sqrt(np.mean(s.fun**2))))
    pb[lab]=out
    P(f'{lab}: noise-map PB fit HCN: FWHM_eff {out[0][4]:.1f}" (centre offset from N-S midpoint {np.hypot(out[0][1]-d["cx"],out[0][2]-d["cy"])*ps:.1f}", rms resid {out[0][5]:.3f}); HCO+: FWHM_eff {out[1][4]:.1f}" (offset {np.hypot(out[1][1]-d["cx"],out[1][2]-d["cy"])*ps:.1f}")')
P('\n# PB gain (relative to the fitted pointing centre) at each position; S/N of the 10" |v|<350 km/s flux per line (stat only)')
R30={r['tag']:r for r in json.load(open('/workspace/out/s30_ratios.json'))}
hdr='pos     r_mid["]  '+'  '.join(f'PB[{b.split()[0]}]' for b in BANDS)+'   | S/N HCN 2-1,3-2,4-3 | S/N HCO+ 2-1,3-2,4-3'
P(hdr)
for pos in ('N','S','Bridge','E1','E2','W','Sblue'):
    g=[];sn_h=[];sn_o=[]
    for lab,tg in BANDS.items():
        d=np.load(f'/workspace/out/s30_maps_{tg}_10as.npz'); x,y=xy_of(pos,d); a,x0,y0,b,fw,_=pb[lab][0]
        g.append(float(np.exp(-b*((x-x0)**2+(y-y0)**2)*float(d['ps'])**2)))
        r=R30[tg+'_10as']['rows'][pos+'_pos']; sn_h.append(r['HCN']/r['eHCN']); sn_o.append(r['HCOp']/r['eHCOp'])
    rr=np.hypot(*(((np.array(REGIONS[pos][:2])-np.array(CEN))*[np.cos(np.deg2rad(CEN[1])),1])*3600))
    P(f'{pos:7s} {rr:6.1f}   '+'   '.join(f'{v:6.2f}' for v in g)+'     | '+' '.join(f'{v:6.1f}' for v in sn_h)+' | '+' '.join(f'{v:6.1f}' for v in sn_o))
    res['pb'][pos]=g; res['snr'][pos]=dict(HCN=sn_h,HCO=sn_o)
P('\n# PB handling test: native pbcor map smoothed to 10" (A) vs native pbcor*PB smoothed to 10" then /PB at the position (B), HCN and HCO+; ratio A/B')
for lab,tg in BANDS.items():
    nat=np.load(f'/workspace/out/s30_maps_{tg}_native.npz'); d10=np.load(f'/workspace/out/s30_maps_{tg}_10as.npz'); ps=float(nat['ps']); bm=nat['beam']; ksig=np.sqrt(max(100-bm[0]*bm[1],0))/2.3548/ps
    n=nat['HCN'].shape[0]; Y,X=np.mgrid[:n,:n]; a,x0,y0,b,fw,_=pb[lab][0]; PBm=np.exp(-b*((X-x0)**2+(Y-y0)**2)*ps**2)
    for pos in ('N','E1','E2','W'):
        x,y=xy_of(pos,nat); row=[]
        for m in ('HCN','HCOp'):
            A=gaussian_filter(np.nan_to_num(nat[m]),ksig); B=gaussian_filter(np.nan_to_num(nat[m])*PBm,ksig)/map_coordinates(PBm,[[y],[x]],order=1)[0]
            row.append(map_coordinates(A,[[y],[x]],order=1)[0]/map_coordinates(B,[[y],[x]],order=1)[0])
        P(f'  {lab:15s} {pos:3s}  A/B  HCN {row[0]:.3f}  HCO+ {row[1]:.3f}')
        res['pbtest'][f'{lab}|{pos}']=row
json.dump(res,open(OUT+'s55a_pb_snr.json','w'),indent=1,default=float); open(OUT+'s55a_pb_snr.txt','w').write('\n'.join(log))
