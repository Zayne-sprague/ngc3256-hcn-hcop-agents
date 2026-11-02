# Q2 round 2, fix 5: the red feature (+180..320 km/s, ~3.4" N of S, 1.8" S of N) treated like the blue lobe.
# (i) a-posteriori trials: r=1.0" aperture grid (0.5" step) within 8" of the N-S midpoint (excluding 1.0" around N and S); the feature flux vs the distribution in the same region of (a) the same window, (b) line-free windows of equal width (local residual+sidelobe noise), (c) the opposite-sign windows.
# (ii) flux profile along the N-S axis (PA 185) in red sub-windows and spectrum of the N-S line; azimuth ring around N at the feature radius (1.8") in the red window and in line-free windows.
# (iii) leakage: same ring/axis tests with the S blue window (-405..-150) (no N/Bridge emission there) as the leakage reference; (iv) N wing limits for 150-350 km/s (blue and red), r=1.0" at N, field noise, as fraction of core.
# Data: Xb4d@1.05" pbcor, per-pixel baseline median(-900..-450).  Output out/q2r2/s77f.txt/json, s77f_maps.png
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); rng=np.random.default_rng(9)
def wflux(v,dv,sp,lo,hi):
    ov=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return float((sp*ov).sum()*dv)
MID=offset_pos(S0,0.08,2.6); res={}
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0); n=R['n']; ps=R['ps']
    def W(lo,hi): k=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv     # Jy/beam km/s map
    wins={'red':(180,320),'blue':(-405,-150),'free_p':(380,520),'free_m':(-800,-660),'free_p2':(540,680)}
    M={k:W(*w) for k,w in wins.items()}
    cx,cy=R['pxy'](*MID); sx,sy=R['pxy'](*S0); nx,ny=R['pxy'](*N0); Y,X=np.mgrid[:n,:n]
    # aperture sums via convolution with a r=1" top-hat (sum/beam area -> Jy km/s)
    from scipy.signal import fftconvolve
    hat=(np.hypot(*np.mgrid[-int(1.2/ps):int(1.2/ps)+1,-int(1.2/ps):int(1.2/ps)+1])*ps<=1.0).astype(float)/R['ba']
    A={k:fftconvolve(m,hat,mode='same') for k,m in M.items()}
    fm=field_pixels(R,rmin=3.5,rmax=22.0); sig={k:float(np.std(A[k][fm])) for k in A}
    P(f'### {mol}: field sigma of r=1" aperture flux (Jy km/s) {json.dumps({k:round(s,3) for k,s in sig.items()})}')
    # trial region
    reg=(np.hypot(X-cx,Y-cy)*ps<8.0)&(np.hypot(X-sx,Y-sy)*ps>1.0)&(np.hypot(X-nx,Y-ny)*ps>1.0)
    ap_area_pix=np.pi*(1.0/ps)**2; nind=reg.sum()*ps*ps/np.pi   # independent r=1" apertures (area/pi arcsec2)
    fpos=offset_pos(S0,0.15,3.3); fx,fy=R['pxy'](*fpos); fv=A['red'][int(round(fy)),int(round(fx))]
    P(f'  feature aperture (r=1" at 3.3" N of S): red {fv:.3f} Jy km/s = {fv/sig["red"]:.1f} field-sigma; region r<8" of the midpoint holds ~{nind:.0f} independent apertures')
    for k in ('red','free_p','free_m','free_p2','blue'):
        a=A[k][reg]; P(f'  region apertures {k:7s}: max {a.max():+.3f} min {a.min():+.3f}; rms {a.std():.3f}; N(>= {fv:.3f})={int((a>=fv).sum())} grid points (~{(a>=fv).sum()*0.0+int((a>=fv).sum())*ps*ps/np.pi:.1f} independent); N(<= -{fv:.3f})={int((a<=-fv).sum())}')
    # local rms from line-free windows in the region (all three) => local significance
    loc=np.sqrt(np.mean([A[k][reg].var()+A[k][reg].mean()**2 for k in ('free_p','free_m','free_p2')]))
    P(f'  local rms (line-free windows, same region) {loc:.3f} -> feature {fv/loc:.1f} sigma_local; one-sided Gaussian p={0.5*__import__("math").erfc(fv/loc/2**0.5):.1e}; x{nind:.0f} trials = {nind*0.5*__import__("math").erfc(fv/loc/2**0.5):.1e}')
    res[f'{mol}|trials']=dict(fv=float(fv),sig_field=sig['red'],loc=float(loc),nind=float(nind),nge=int((A['red'][reg]>=fv).sum()))
    # position of the maximum of the red-window aperture map inside the region (excluding 1" at nuclei)
    a=np.where(reg,A['red'],-1e9); y,x=np.unravel_index(np.argmax(a),a.shape); P(f'  maximum red-window aperture in region at ({-(x-sx)*ps:+.2f}" E, {(y-sy)*ps:+.2f}" N of S): {A["red"][y,x]:.3f}')
    # (ii) flux along the N-S axis
    th=np.deg2rad(185.0); ex=np.array([np.sin(th),np.cos(th)])
    offs=np.arange(-3,9.01,0.5); prof={}
    for nm,(lo,hi) in {'red180-250':(180,250),'red250-320':(250,320),'red180-320':(180,320),'blue-405..-150':(-405,-150),'free+380..520':(380,520)}.items():
        mm=W(lo,hi); row=[]
        for o in offs:
            pos=offset_pos(S0,-o*ex[0],-o*ex[1]); xx,yy=R['pxy'](*pos); m=apmask(n,ps,(xx,yy),1.0); row.append(float(mm[m].sum()/R['ba']))
        prof[nm]=row
    res[f'{mol}|axis']=dict(offs=offs.tolist(),prof=prof)
    P(f'  along N-S axis (offset N of S, r=1" apertures, Jy km/s):')
    for nm,row in prof.items(): P(f'    {nm:15s} '+' '.join(f'{o:+.1f}:{x:+.2f}' for o,x in zip(offs,row)))
    # azimuth ring around N at 1.8" and S-lobe ring: red window, line-free, blue
    for cen,cn,rad in ((N0,'N',1.8),(N0,'N',3.4),(S0,'S',1.9)):
        px,py=R['pxy'](*cen); row={}
        for nm in ('red','blue','free_p','free_m'):
            r_=[]
            for az in np.arange(0,360,30):
                x=px+rad/ps*np.sin(np.deg2rad(az)); y=py+rad/ps*np.cos(np.deg2rad(az)); m=apmask(n,ps,(x,y),1.0); r_.append(float(M[nm][m].sum()/R['ba']))
            row[nm]=r_
        res[f'{mol}|ring|{cn}|{rad}']=row
        P(f'  ring r={rad}" around {cn} (az from N, E-ward, 30 deg steps) red: '+' '.join(f'{x:+.2f}' for x in row['red'])+' | blue: '+' '.join(f'{x:+.2f}' for x in row['blue'])+' | free(+380..520): '+' '.join(f'{x:+.2f}' for x in row['free_p']))
    # (iv) N wing limits 150-350 and core
    for pn,pos in (('N',N0),('S',S0)):
        m=apmask(n,ps,R['pxy'](*pos),1.0); sp=aspec(c,R['ba'],m); core=wflux(v,dv,sp,-150,150)
        out={}
        for nm,(lo,hi) in {'blue150-350':(-350,-150),'red150-350':(150,350),'red180-320':(180,320),'blue-405..-150':(-405,-150)}.items():
            f=wflux(v,dv,sp,lo,hi); s_=float(np.std(A['red'][fm])); out[nm]=(f,f/core)
        P(f'  {pn} r=1": core(-150..150) {core:.2f} Jy km/s; '+'; '.join(f'{nm} {f:+.3f} ({100*fr:+.2f}% of core; field sigma {np.std(A["red"][fm]):.3f} => 3sigma upper {3*np.std(A["red"][fm]):.3f} = {100*3*np.std(A["red"][fm])/core:.2f}%)' for nm,(f,fr) in out.items()))
        res[f'{mol}|wing|{pn}']=dict(core=core,**{k:v for k,v in out.items()},sigma=float(np.std(A['red'][fm])))
    if mol=='HCN':
        fig,ax=plt.subplots(1,3,figsize=(18,6)); ext=[-(0-sx)*ps+0,0,0,0]
        for a_,(k,t) in zip(ax,(('red','HCN +180..320'),('blue','HCN -405..-150'),('free_p','HCN +380..520 (line-free)'))):
            sub=M[k]; xs0,xs1=int(sx-9/ps),int(sx+9/ps); ys0,ys1=int(sy-4/ps),int(sy+9/ps)
            a_.imshow(sub[ys0:ys1,xs0:xs1],origin='lower',extent=[-(xs0-sx)*ps,-(xs1-sx)*ps,(ys0-sy)*ps,(ys1-sy)*ps],cmap='viridis',vmin=-0.1,vmax=0.4); a_.set_title(t+' Jy/b km/s'); a_.plot([0,-(nx-sx)*ps],[0,(ny-sy)*ps],'r+',ms=12)
            a_.contour(-(np.arange(xs0,xs1)-sx)*ps,(np.arange(ys0,ys1)-sy)*ps,M['red' if k!='blue' else 'blue'][ys0:ys1,xs0:xs1],levels=[3*sig[k]/1.0*0+0.1,0.2,0.3],colors='w',linewidths=.6)
            a_.set_xlabel('E offset from S (")'); a_.set_ylabel('N offset (")')
        plt.tight_layout(); plt.savefig(OUT+'/s77f_maps.png',dpi=70)
json.dump(res,open(OUT+'/s77f.json','w'),indent=1,default=float); open(OUT+'/s77f.txt','w').write('\n'.join(txt))
