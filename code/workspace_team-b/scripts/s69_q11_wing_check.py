# Q11: do team-a p0005-r02/r03, p0006, p0009, p0010 claims hold in our data?  (i) blue-half HCN/HCO+ 3-2 in fixed 2" apertures (Soutflow, N, Bridge, S) in Xb4d and Xb4b@1.05" with
# local-centroid velocities; (ii) wing ratio vs beam at N (Xb4b native 0.25", 0.30", 1.05"; velocity slope, radial dependence); (iii) blend/identity test of the HCN excess over k*HCO+ at HCN-frame
# 2375-2575 km/s; (iv) N wing limits.  Reuses s30 loaders (cropped pbcor cubes, analytic Fourier convolution with flux rescaled by beam-area ratio).
# Baseline: env Q11_BASE A (default; per-pixel median of -900<v-VSYS<-450), B (linear fit incl. 480..700), C (median v<-800) (cube covers -1000..+700); ratios are integrated Jy km/s ratios; R_Tb = R_Jy * 1.0126.
# Output: out/q11/*.json, out/q11/s69.txt   (run: python3 -I scripts/s69_q11_wing_check.py)
import sys, json, os; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
from s30_matchedJ_ratios import *
from scipy.optimize import curve_fit
from scipy.ndimage import binary_dilation
S.HALF=10.0
S.SPW[('Xb4b',3)]=('spw27','spw25')
OUT='/workspace/out/q11'; os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(11)
KT=1.0126
txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
SB_ofs=(-0.2,-1.85)
def pos(name):
    if name=='Soutflow':
        r0=REGIONS['S'][:2]; return (r0[0]+SB_ofs[0]/3600/np.cos(np.deg2rad(r0[1])), r0[1]+SB_ofs[1]/3600)
    return REGIONS[name][:2]
POS=['Soutflow','N','Bridge','S']

BASE=os.environ.get('Q11_BASE','A')   # A: per-pixel median of -900<v<-450 (local, blue side); B: per-pixel linear fit to -900<v<-450 and 480<v<700 ; C: median v<-800 (s30-like)
def baseline(cube,v):
    if BASE=='A': return np.median(cube[(v>-900)&(v<-450)],axis=0)
    if BASE=='C': return np.median(cube[v<-800],axis=0)
    ref=((v>-900)&(v<-450))|((v>480)&(v<700)); A=np.vstack([np.ones(ref.sum()),v[ref]/1000.]).T
    co=np.linalg.pinv(A)@cube[ref].reshape(ref.sum(),-1); full=np.vstack([np.ones(len(v)),v/1000.]).T@co
    return full.reshape(cube.shape).astype(np.float32)
def prep(ous,fwhm):
    """cubes (baseline subtracted, v<-800 ref) of HCN and HCO+ at circular beam fwhm (None=native)"""
    o={}
    for mol,spw in zip(('HCN','HCO+'),S.SPW[(ous,3)]):
        cube,v,ps,wc,x0,y0,b0=S.load(ous,spw,S.F0[mol][3])
        if fwhm: cube=S.conv(cube,ps,b0,(fwhm,fwhm,0.0)); tgt=(fwhm,fwhm,0.0)
        else: tgt=b0
        cube=cube-baseline(cube,v)
        n=cube.shape[1]
        o[mol]=dict(c=cube,v=v,dv=abs(np.median(np.diff(v))),ps=ps,ba=1.1331*tgt[0]*tgt[1]/ps**2,beam=tgt,
                    pxy=(lambda ra,dec,wc=wc,x0=x0,y0=y0:(wc.world_to_pixel_values(ra,dec)[0]-x0,wc.world_to_pixel_values(ra,dec)[1]-y0)))
    return o
def apmask(R,xy,rad,rin=0.0):
    n=R['c'].shape[1]; Y,X=np.mgrid[:n,:n]; d=np.hypot(X-xy[0],Y-xy[1])*R['ps']; return (d<=rad)&(d>=rin)
def spec(R,m): return R['c'][:,m].sum(1)/R['ba']        # Jy/beam summed over aperture / beam area -> Jy (per channel)
def win(R,sp,a,b,vloc=0.0):
    v=R['v']-vloc; k=(v>=a)&(v<=b); return sp[k].sum()*R['dv'],k.sum()
def gfit(v,s,c0=0.0):
    f=lambda x,a,m,w:a*np.exp(-0.5*((x-m)/w)**2)
    try:
        p,_=curve_fit(f,v,s,p0=[s.max(),c0,50.0],bounds=([0,-200,15],[np.inf,200,250])); return p
    except Exception: return np.array([np.nan,np.nan,np.nan])
def randpos(o,n,rad,rmin,rmax,excl=()):
    H=o['HCN']; Y,X=np.mgrid[:H['c'].shape[1],:H['c'].shape[1]]
    m0=H['c'][np.abs(H['v'])<350].sum(0)*H['dv']; sg=np.median(np.abs(m0-np.median(m0)))*1.4826
    bright=binary_dilation(m0>4*sg,iterations=int(1.5/H['ps']))
    cx,cy=H['pxy'](*S.CEN); out=[]
    while len(out)<n:
        x,y=rng.uniform(0,H['c'].shape[1],2); d=np.hypot(x-cx,y-cy)*H['ps']
        if rmin<d<rmax and not bright[int(y),int(x)] and all(np.hypot(x-e[0],y-e[1])*H['ps']>2*rad for e in excl): out.append((x,y))
    return out
def noise_vec(o,RP,rad,wins):
    """paired window sums (HCN,HCO+) at random positions, windows = list of (a,b) relative to VSYS"""
    res=np.zeros((len(RP),len(wins),2))
    for i,xy in enumerate(RP):
        for j,mol in enumerate(('HCN','HCO+')):
            R=o[mol]; m=apmask(R,xy,rad); sp=spec(R,m)
            for k,(a,b) in enumerate(wins): res[i,k,j]=win(R,sp,a,b)[0]
    return res
def ratio_mc(h,c,nz,N=20000):
    """MC ratio from paired empirical noise samples nz (n,2); returns median,16,84 + P"""
    i=rng.integers(0,len(nz),N); hh=h+nz[i,0]; cc=c+nz[i,1]
    ok=(hh>0)&(cc>0); r=hh[ok]/cc[ok]; q=np.percentile(r,[16,50,84]) if ok.sum()>100 else [np.nan]*3
    return q,ok.mean()
RES={}
# =========================================================== Part 1: fixed 2" aperture blue-half ratios
P('# Part 1: fixed aperture r=1.0" (2" diameter), local-centroid velocity (Gaussian centroid of HCO+ spectrum, +-250 km/s), blue half v_loc-400..-100, red half +100..+400, core |v_loc|<100. R_Tb=R_Jy*1.0126.')
P('# Noise: 60 random line-free apertures 5-9.5" from the N-S midpoint (bright emission masked), same windows, paired HCN/HCO+ draws (includes baseline/correlated noise).')
P1={}
for ous,fw in (('Xb4d',1.05),('Xb4d',None),('Xb4b',1.05)):
    o=prep(ous,fw); tag=f'{ous}_'+(f'{fw}as' if fw else 'native'); H,O=o['HCN'],o['HCO+']
    RP=randpos(o,60,1.0,5.0,9.5,excl=[H['pxy'](*pos(p)) for p in POS])
    P(f'## {tag}: beam {H["beam"][0]:.2f}x{H["beam"][1]:.2f} (HCO+ {O["beam"][0]:.2f}x{O["beam"][1]:.2f}), noise positions {len(RP)}')
    P1[tag]={}
    for pn in POS:
        xy=H['pxy'](*pos(pn)); mH=apmask(H,xy,1.0); mO=apmask(O,xy,1.0); sH=spec(H,mH); sO=spec(O,mO)
        sel=np.abs(O['v'])<250; p=gfit(O['v'][sel],sO[sel]); vloc=float(p[1]) if np.isfinite(p[1]) else 0.0
        pH=gfit(H['v'][np.abs(H['v'])<250],sH[np.abs(H['v'])<250]);
        wins=[(vloc-400,vloc-100),(vloc+100,vloc+400),(vloc-100,vloc+100)]
        # noise windows relative to VSYS are fine (noise stationary); use windows relative to v=0 of same widths
        nz=noise_vec(o,RP,1.0,[(-400,-100),(100,400),(-100,100)])
        row=dict(vloc=vloc,vloc_HCN=float(pH[1]),sigma_HCO=float(p[2]))
        for k,nm in enumerate(('blue','red','core')):
            h,_=win(H,sH,wins[k][0],wins[k][1]); c,_=win(O,sO,wins[k][0],wins[k][1]); eh,ec=nz[:,k,0].std(),nz[:,k,1].std()
            q,fok=ratio_mc(h,c,nz[:,k,:])
            row[nm]=dict(HCN=h,eHCN=eh,HCO=c,eHCO=ec,R_Tb=h/c*KT if c>0 else np.nan,R_Tb_q16_50_84=[float(x*KT) for x in q],frac_pos=float(fok))
        # additional: 'VSYS-referenced' -400..-150 blue window (our p0005 recipe) and 2375-2575 km/s excess window (-405..-205)
        for nm,(a,b) in (('blue_p0005_-400_-150',(-400,-150)),):
            h,_=win(H,sH,a,b); c,_=win(O,sO,a,b); nzz=noise_vec(o,RP,1.0,[(a,b)])[:,0,:]
            q,fok=ratio_mc(h,c,nzz); row[nm]=dict(HCN=h,eHCN=nzz[:,0].std(),HCO=c,eHCO=nzz[:,1].std(),R_Tb=h/c*KT if c>0 else np.nan,R_Tb_q16_50_84=[float(x*KT) for x in q])
        P1[tag][pn]=row
        b,r_,c_=row['blue'],row['red'],row['core']
        P(f'  {pn:8s} vloc={vloc:+6.1f} (HCN {pH[1]:+.0f}, sigma_HCO+ {p[2]:.0f})  blue: HCN {b["HCN"]:.3f}+-{b["eHCN"]:.3f} ({b["HCN"]/b["eHCN"]:.1f}s) HCO+ {b["HCO"]:.3f}+-{b["eHCO"]:.3f} ({b["HCO"]/b["eHCO"]:.1f}s) R_Tb={b["R_Tb"]:.2f} [{b["R_Tb_q16_50_84"][0]:.2f},{b["R_Tb_q16_50_84"][1]:.2f},{b["R_Tb_q16_50_84"][2]:.2f}]  | red R_Tb={r_["R_Tb"]:.2f} ({r_["HCN"]/r_["eHCN"]:.1f}s,{r_["HCO"]/r_["eHCO"]:.1f}s) | core R_Tb={c_["R_Tb"]:.3f}+-{c_["R_Tb"]*np.hypot(c_["eHCN"]/c_["HCN"],c_["eHCO"]/c_["HCO"]):.3f}   p0005-recipe(-400..-150 VSYS) R_Tb={row["blue_p0005_-400_-150"]["R_Tb"]:.2f}')
    # difference tests (Soutflow vs N, Bridge) on blue half via paired MC with independent noise
    for ref in ('N','Bridge','S'):
        a=P1[tag]['Soutflow']['blue']; bb=P1[tag][ref]['blue']
        def draw(d,N=20000):
            h=rng.normal(d['HCN'],d['eHCN'],N); c=rng.normal(d['HCO'],d['HCO']*0+d['eHCO'],N); ok=(h>0)&(c>0); return np.where(ok,h/np.where(c>0,c,np.nan),np.nan)*KT
        ra,rb=draw(a),draw(bb); dd=ra-rb; dd=dd[np.isfinite(dd)]
        P(f'    Soutflow - {ref} blue-half R_Tb: median {np.median(dd):+.2f} [16-84 {np.percentile(dd,16):+.2f},{np.percentile(dd,84):+.2f}], P(Soutflow>{ref})={np.mean(dd>0):.2f}; with 10% systematic floor on each ratio the z-score is {np.median(dd)/np.hypot(0.5*(np.percentile(dd,84)-np.percentile(dd,16)),0.1*np.hypot(a["R_Tb"],bb["R_Tb"])):+.2f}')
    if tag=='Xb4d_1.05as': KEEP=o
    del o
RES['part1']=P1
# =========================================================== Part 2: beam test at N with Xb4b
P('\n# Part 2: N, Xb4b. Ratio vs local-centroid velocity in 40 km/s bins (+-140), linear slope per 100 km/s; same data at native 0.25", 0.30", 1.05".')
P2={}
BINS=[(a,a+40) for a in range(-140,140,40)]
for fw,rad in ((None,0.25),(0.30,0.25),(1.05,0.5),(1.05,0.25)):
    o=prep('Xb4b',fw); H,O=o['HCN'],o['HCO+']; tag=('native' if fw is None else f'{fw}as')+f'_r{rad}'
    xy=H['pxy'](*pos('N')); m=apmask(H,xy,rad); sH=spec(H,m); sO=spec(O,apmask(O,xy,rad))
    p=gfit(O['v'][np.abs(O['v'])<250],sO[np.abs(O['v'])<250]); vloc=p[1]
    RP=randpos(o,80,max(rad,0.25),4.0,9.5,excl=[H['pxy'](*pos(q)) for q in POS])
    wins=[(a,b) for a,b in BINS]+[(100,260),(-260,-100),(-60,60)]
    nz=noise_vec(o,RP,rad,wins)
    R=[];E=[];VC=[]
    for k,(a,b) in enumerate(BINS):
        h,_=win(H,sH,a,b,vloc); c,_=win(O,sO,a,b,vloc); r=h/c*KT
        e=r*np.hypot(nz[:,k,0].std()/h,nz[:,k,1].std()/c); R.append(r);E.append(e);VC.append(0.5*(a+b))
    R=np.array(R);E=np.array(E);VC=np.array(VC)
    w=1/E**2; A=np.vstack([np.ones_like(VC),VC/100.]).T; cov=np.linalg.inv(A.T@(A*w[:,None])); sl=cov@(A.T@(w*R))
    stk={}
    for nm,(a,b),k in (('red100_260',(100,260),len(BINS)),('blue-260_-100',(-260,-100),len(BINS)+1),('core',(-60,60),len(BINS)+2)):
        h,_=win(H,sH,a,b,vloc); c,_=win(O,sO,a,b,vloc); r=h/c*KT
        stk[nm]=dict(HCN=h,eHCN=nz[:,k,0].std(),HCO=c,eHCO=nz[:,k,1].std(),R=r,eR=abs(r)*np.hypot(nz[:,k,0].std()/h,nz[:,k,1].std()/c))
    hb=win(H,sH,100,260,vloc)[0]+win(H,sH,-260,-100,vloc)[0]; cb=win(O,sO,100,260,vloc)[0]+win(O,sO,-260,-100,vloc)[0]
    ehb=(nz[:,len(BINS),0]+nz[:,len(BINS)+1,0]).std(); ecb=(nz[:,len(BINS),1]+nz[:,len(BINS)+1,1]).std()
    rb=hb/cb*KT; erb=abs(rb)*np.hypot(ehb/hb,ecb/cb); rc=stk['core']['R']; erc=stk['core']['eR']
    stk['both100_260']=dict(HCN=hb,eHCN=ehb,HCO=cb,eHCO=ecb,R=rb,eR=erb,wing_over_core=rb/rc,e_wing_over_core=rb/rc*np.hypot(erb/rb,erc/rc))
    P(f'      both-sides 100-260 stack: R_wing={rb:.2f}+-{erb:.2f} (HCN {hb/ehb:.1f}s, HCO+ {cb/ecb:.1f}s)  R_core={rc:.2f}+-{erc:.2f}  wing/core={rb/rc:.2f}+-{rb/rc*np.hypot(erb/rb,erc/rc):.2f}')
    P2[tag]=dict(vloc=float(vloc),R_bins=R.tolist(),eR=E.tolist(),slope_per100=float(sl[1]),eslope=float(np.sqrt(cov[1,1])),intercept=float(sl[0]),stack=stk)
    P(f'## {tag}: beam {H["beam"][0]:.2f}, vloc(HCO+)={vloc:+.0f}; bins R_Tb '+' '.join(f'{r:.2f}+-{e:.2f}' for r,e in zip(R,E))+f' | slope {sl[1]:+.3f}+-{np.sqrt(cov[1,1]):.3f} per 100 km/s')
    P('      stacks: '+'; '.join(f'{k} R={v["R"]:.2f}+-{v["eR"]:.2f} (HCN {v["HCN"]/v["eHCN"]:.1f}s, HCO+ {v["HCO"]/v["eHCO"]:.1f}s)' for k,v in stk.items()))
    if fw is None:
        # radial dependence of the 100-260 (both sides) wing ratio at native resolution
        P('      radial wing/core ratio at native 0.25": annuli (arcsec) -> R_core(|v|<60), R_wing(60<|v|<260 stacked both sides, v_loc), S/N HCN,HCO+ wing')
        for rin,rout in ((0,0.25),(0.25,0.5),(0.5,1.0),(1.0,2.0)):
            mH=apmask(H,xy,rout,rin); mO=apmask(O,xy,rout,rin); sh=spec(H,mH); so=spec(O,mO)
            core=(win(H,sh,-60,60,vloc)[0],win(O,so,-60,60,vloc)[0])
            wing=(win(H,sh,100,260,vloc)[0]+win(H,sh,-260,-100,vloc)[0],win(O,so,100,260,vloc)[0]+win(O,so,-260,-100,vloc)[0])
            RPa=randpos(o,60,rout,4.0,9.5,excl=[H['pxy'](*pos(q)) for q in POS]) if rout<=1.0 else []
            if RPa:
                nzr=np.zeros((len(RPa),2))
                for i,xx in enumerate(RPa):
                    for j,R_ in enumerate((H,O)):
                        s_=spec(R_,apmask(R_,xx,rout,rin)); nzr[i,j]=win(R_,s_,100,260)[0]+win(R_,s_,-260,-100)[0]
                eh,ec=nzr[:,0].std(),nzr[:,1].std()
            else: eh=ec=np.nan
            rr=wing[0]/wing[1]*KT if wing[1]>0 else np.nan
            P(f'        {rin:.2f}-{rout:.2f}": core R={core[0]/core[1]*KT:.2f}  wing HCN {wing[0]:.3f}+-{eh:.3f} HCO+ {wing[1]:.3f}+-{ec:.3f} -> R_wing={rr:.2f}+-{abs(rr)*np.hypot(eh/wing[0],ec/wing[1]):.2f}')
            P2[tag].setdefault('radial',[]).append(dict(rin=rin,rout=rout,Rcore=core[0]/core[1]*KT,HCNw=wing[0],eHCNw=eh,HCOw=wing[1],eHCOw=ec,Rwing=rr))
    del o
RES['part2']=P2
# ---- mixing prediction: convolve native-beam radial profile model.  The convolved result IS the data; here a model-independent decomposition:
# the 1.05" r<0.5" aperture wing ratio vs the ratio obtained by summing the native cube over the same area (identical up to beam-area weighting)
json.dump(RES,open(OUT+f'/s69_part12_base{BASE}.json','w'),indent=1,default=float)
open(OUT+f'/s69_part12_base{BASE}.txt','w').write('\n'.join(txt))
