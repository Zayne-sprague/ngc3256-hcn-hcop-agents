# Q20 step b: matched-beam (4.5"x3.0" PA 87) window maps for the 3-2 (Xb4d) and 4-3/7-6 (Xb47) lines, plus the SO(6_5-5_4) and H13CN windows.
# Maps are baseline-subtracted Jy/beam*km/s (per-pixel baseline = median of line-free channels, listed below), windows: core |v|<150, red +180..320, blue -320..-180 (channel edges honoured),
# flux rescaled to the common beam; B7 maps resampled (bilinear) onto the B7 pixel grid (0.53").   Saves npz for s110c.
# usage: python3 -I s110b_q20_maps.py
import sys; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.ndimage import map_coordinates
OUT='/workspace/out/q20'; CB=(4.5,3.0,87.0)
LINES={ # name: (ous,spw,f0,baseline channel ranges (km/s))
 'HCN32':('Xb4d','spw27',265.88618e9,[(-900,-450)]),'HCO+32':('Xb4d','spw25',267.55763e9,[(-900,-450)]),
 'HCN43':('Xb47','spw16',354.50547e9,[(-900,-450)]),'HCO+43':('Xb47','spw18',356.73424e9,[(-630,-450),(450,880)]),
 'CS76':('Xb47','spw22',342.88285e9,[(-840,-450),(450,720)]),'SO65':('Xb4d','spw31',251.82577e9,[(-900,-450),(450,600)]),
 'H13CN43':('Xb47','spw20',345.33977e9,[(450,1000)])}
txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
def win(v,dv,a,b): return np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv
maps={}
# reference B7 grid (HCN43)
c0,v0,ps7,wc7,x07,y07,b07=load_cube('Xb47','spw16',354.50547e9,-100,100,half=24.0); n7=c0.shape[1]
Y,X=np.mgrid[:n7,:n7]; wra,wdec=wc7.pixel_to_world_values(X+x07,Y+y07)
for nm,(o,sp,f0,bl) in LINES.items():
    c,v,ps,wc,x0,y0,b0=load_cube(o,sp,f0,-1500,1500,half=24.0); dv=abs(np.median(np.diff(v)))
    m=np.zeros(len(v),bool)
    for a,b in bl: m|=(v>a)&(v<b)
    if nm=='H13CN43':   # CO 3-2 red-wing tail: linear baseline from 450..1000 (see s110d for systematics)
        mm=(v>450)&(v<1000); A=np.vstack([np.ones(mm.sum()),v[mm]]).T; co=np.linalg.lstsq(A,c[mm].reshape(mm.sum(),-1),rcond=None)[0]
        c=c-(np.outer(np.ones(len(v)),co[0])+np.outer(v,co[1])).reshape(c.shape)
    else: c=c-np.median(c[m],axis=0)
    out={}
    for wn,(a,b) in {'core':(-150,150),'red':(180,320),'blue':(-320,-180),'red2':(200,300),'blue2':(-300,-200),'redf':(400,540),'bluef':(-540,-400),**{f'r{a}':(a,a+40) for a in range(120,400,40)},**{f'b{a}':(-a-40,-a) for a in range(120,400,40)}}.items():
        k=win(v,dv,a,b); mp=np.tensordot(k,c,axes=(0,0))*dv
        if o=='Xb4d': mp=conv(mp[None],ps,b0,CB)[0]
        else: mp=conv(mp[None],ps,b0,CB)[0]
        if o=='Xb4d':
            px=pxy_fn(wc,x0,y0)(wra,wdec); mp=map_coordinates(mp,[px[1],px[0]],order=1,cval=np.nan)
        out[wn]=mp.astype(np.float32)
    # off-source sigma of the window map (field 9..14" from MID)
    mid=offset_pos(S0,0.1,2.55); mx,my=pxy_fn(wc7,x07,y07)(*mid); rr=np.hypot(X-mx,Y-my)*ps7
    # channel rms in native cube
    sg={wn:float(np.nanstd(out[wn][(rr>9)&(rr<14)])) for wn in out}
    P(f'{nm}: native beam {b0[0]:.2f}x{b0[1]:.2f}, dv {dv:.2f}, baseline chans {int(m.sum()) if nm!="H13CN43" else "lin450-1000"}; off-source sd of window maps at 4.5x3.0 (Jy/beam km/s): '+', '.join(f'{k} {s:.4f}' for k,s in sg.items()))
    maps[nm]=out; maps[nm+'_sg']=sg
np.savez(f'{OUT}/s110b_maps.npz',**{f'{nm}|{w}':maps[nm][w] for nm in LINES for w in maps[nm]},ps=ps7,x0=x07,y0=y07,**{f'sg|{nm}|{w}':maps[nm+'_sg'][w] for nm in LINES for w in maps[nm+'_sg']},wra=wra,wdec=wdec)
open(OUT+'/s110b.txt','w').write('\n'.join(txt))
