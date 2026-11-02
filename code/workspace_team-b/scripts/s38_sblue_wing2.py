# Q1 (round 2): calibrated test of the Sblue blue-wing HCN/HCO+ 3-2 ratio (R~1.4) in the S-nucleus outflow direction. Xb4d native cubes (beams 0.94x0.74/0.93x0.74).
# Aperture centre = S + (-0.2",-1.85") (s32), r=0.75" (primary), per-pixel baseline from line-free channels (s30 products). dv=4.4 km/s. Fluxes in Jy km/s per beam (aperture sum / beam area in pixels).
# T1 null: same pipeline at 400 random line-free positions (noise) and at bright positions N,S,Bridge,E1,E2 (windows blue -400..-150 relative to VSYS): distribution of wing R, false-positive rates.
# T2 posterior of R_wing (independent truncated Gaussians h,c>0, flat prior, sigma = empirical aperture scatter at random positions), 68/95% intervals, P(R_wing>R_core), lower limit on R_wing-R_core.
# T3 robustness: window variants, aperture centre +-0.2", radius; spectra (aperture + peak pixel) -> out/s38_spectra.png/.txt
# T4 replication with Xb4b convolved to the Xb4d beam (s11 maps) and Xb4f (blended) ; T5 neighbouring lines in the windows (out/linelist.csv); T6 'peak-pixel-per-line' recipe bias by injection into noise patches
# T7 RADEX: abundance ratio needed for R_K=1.4 vs core at the same excitation.   Output out/s38_sblue2.json/.txt
import sys, json, subprocess, os; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
from s30_matchedJ_ratios import *
from scipy.ndimage import binary_dilation
rng=np.random.default_rng(3)
o=S.run('Xb4d',3,None,windows=(350,)); H,O=o['HCN'],o['HCO+']
ps=H['ps']; ba=H['beamarea']; n=H['c'].shape[1]; Y,X=np.mgrid[:n,:n]; dv=H['dv']
S0=REGIONS['S'][:2]; cosd=np.cos(np.deg2rad(S0[1]))
SB=(S0[0]-0.2/3600/cosd, S0[1]-1.85/3600)
txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
_wc={}
def wmap(R,a,b):
    k=(id(R),a,b)
    if k not in _wc: m=(R['v']>=a)&(R['v']<=b); _wc[k]=R['c'][m].sum(0)*dv
    return _wc[k]
def apsum(m,x,y,rad): return m[np.hypot(X-x,Y-y)*ps<=rad].sum()/ba
def pxy(R,ra,dec): return R['pxy'](ra,dec)
bright=binary_dilation(H['m0_350']>5*H['err_350'],iterations=int(2.5*0.94/ps))
cx,cy=H['pxy'](*CEN); RND=[]
while len(RND)<400:
    x,y=rng.uniform(0,n,2)
    if 6<np.hypot(x-cx,y-cy)*ps<18 and not bright[int(y),int(x)]: RND.append((x,y))
def noise(a,b,rad):
    return [ (apsum(wmap(H,a,b),x,y,rad),apsum(wmap(O,a,b),x,y,rad)) for x,y in RND]
W0=(-400,-150); RAD=0.75
res={}
# ---- T1 null
nz=np.array(noise(*W0,RAD)); sH,sO=nz[:,0].std(),nz[:,1].std()
P(f'## T1 noise at 400 random line-free positions, window {W0}, r={RAD}": sigma HCN {sH:.3f}, HCO+ {sO:.3f} Jy km/s/beam; mean {nz.mean(0)}; corr(HCN,HCO+) {np.corrcoef(nz.T)[0,1]:.2f}')
# noise-only R: fraction of random positions with both >3 sigma, and R>1.2
both=(nz[:,0]>3*sH)&(nz[:,1]>3*sO); P(f'   random positions with both wings >3 sigma: {both.sum()}/400; with HCN>3s&HCO+>2s: {((nz[:,0]>3*sH)&(nz[:,1]>2*sO)).sum()}; R>1.2 among both>2sigma positive: {((nz[:,0]>2*sH)&(nz[:,1]>2*sO)&(nz[:,0]/np.maximum(nz[:,1],1e-9)>1.2)).sum()}/400')
# bright positions: wing leakage / real wings
BP={k:pxy(H,*REGIONS[k][:2]) for k in ('N','S','Bridge','E1','E2')}; BP['Sblue']=pxy(H,*SB)
P('## T1b bright positions, same windows/aperture: wing flux and R (nuclear core leakage / real wing)')
for k,(x,y) in BP.items():
    h=apsum(wmap(H,*W0),x,y,RAD); c=apsum(wmap(O,*W0),x,y,RAD); hc=apsum(wmap(H,-150,150),x,y,RAD); cc=apsum(wmap(O,-150,150),x,y,RAD)
    R_=h/c if (h>0 and c>0) else np.nan
    P(f'  {k:7s} wing HCN {h:+.3f} ({h/sH:+.1f}s) HCO+ {c:+.3f} ({c/sO:+.1f}s) R_wing={R_:.2f}   core R={hc/cc:.2f}   wing/core HCN {h/hc:.3f} HCO+ {c/cc:.3f}')
    res[f'bright_{k}']=dict(h=h,c=c,hc=hc,cc=cc)
# detections-selected null: take all pixels positions (grid) where both wings detect >4 sigma (aperture r=0.75) in the field -> count of independent such positions
mh=wmap(H,*W0); mc=wmap(O,*W0)
from scipy.ndimage import uniform_filter
k=int(2*RAD/ps)|1; ah=uniform_filter(mh,k)*np.pi*(RAD/ps)**2/ba; ac=uniform_filter(mc,k)*np.pi*(RAD/ps)**2/ba
sel=(ah>4*sH)&(ac>4*sO)&(np.hypot(X-cx,Y-cy)*ps<15)
lab_n=0
from scipy.ndimage import label
lb,nl=label(sel); P(f'## T1c field-wide: contiguous regions (r<15" of midpoint) with both blue-wing aperture fluxes >4 sigma: {nl}')
for i in range(1,nl+1):
    m=lb==i; j=np.argmax(np.where(m,ah/sH+ac/sO,-1e9)); yy,xx=np.unravel_index(j,ah.shape)
    P(f'     region {i}: {m.sum()} px  centre offset from N (E,N) arcsec ({-(xx-BP["N"][0])*ps:+.1f},{(yy-BP["N"][1])*ps:+.1f}) HCN {ah[yy,xx]:.3f} HCO+ {ac[yy,xx]:.3f} R={ah[yy,xx]/ac[yy,xx]:.2f}')
# ---- T2 posterior
def post(h,c,sh,sc,N=400000,rcore=None,ecore=None):
    a=rng.normal(h,sh,N*3); a=a[a>0][:N]; b=rng.normal(c,sc,N*3); b=b[b>0][:N]; m=min(len(a),len(b)); r=a[:m]/b[:m]
    q=np.percentile(r,[2.5,16,50,84,97.5]); out=dict(q=q.tolist(),P_gt1=float((r>1).mean()),P_gt1_2=float((r>1.2).mean()))
    if rcore is not None:
        rc=rng.normal(rcore,ecore,m); d=r-rc; out['diff_q']=np.percentile(d,[5,16,50,84]).tolist(); out['P_wing_gt_core']=float((d>0).mean())
    return out
x,y=BP['Sblue']; h=apsum(wmap(H,*W0),x,y,RAD); c=apsum(wmap(O,*W0),x,y,RAD)
hc=apsum(wmap(H,-150,150),x,y,RAD); cc=apsum(wmap(O,-150,150),x,y,RAD)
# core noise: scatter of core-window aperture sums at the random positions
nzc=np.array([(apsum(wmap(H,-150,150),x_,y_,RAD),apsum(wmap(O,-150,150),x_,y_,RAD)) for x_,y_ in RND]); sHc,sOc=nzc[:,0].std(),nzc[:,1].std()
Rcore=hc/cc; eRcore=Rcore*np.hypot(sHc/hc,sOc/cc)
pp=post(h,c,sH,sO,rcore=Rcore,ecore=eRcore)
P(f'## T2 Sblue r={RAD}": wing HCN {h:.3f}+-{sH:.3f} ({h/sH:.1f}s) HCO+ {c:.3f}+-{sO:.3f} ({c/sO:.1f}s); R_wing point {h/c:.2f}; posterior (positive fluxes, flat prior) quantiles 2.5/16/50/84/97.5: {np.round(pp["q"],2)}; P(R>1)={pp["P_gt1"]:.3f}; P(R>1.2)={pp["P_gt1_2"]:.3f}')
P(f'   core: HCN {hc:.3f} HCO+ {cc:.3f} R_core {Rcore:.2f}+-{eRcore:.2f}; R_wing-R_core quantiles 5/16/50/84: {np.round(pp["diff_q"],2)}; P(R_wing>R_core)={pp["P_wing_gt_core"]:.3f}')
res['posterior']=dict(h=h,c=c,sH=sH,sO=sO,Rcore=Rcore,eRcore=eRcore,**pp)
# ---- T3 variants
P('## T3 robustness (R_wing point, +-stat; sigma from random positions in that window/aperture)')
var=[]
for (a,b) in ((-400,-150),(-350,-100),(-300,-200),(-450,-150),(-500,-150)):
    for rad in (0.5,0.75,1.0):
        zz=np.array(noise(a,b,rad)); sh_,so_=zz[:,0].std(),zz[:,1].std()
        hh=apsum(wmap(H,a,b),x,y,rad); co=apsum(wmap(O,a,b),x,y,rad); r=hh/co if hh>0 and co>0 else np.nan
        var.append(dict(win=(a,b),rad=rad,h=hh,c=co,sh=sh_,sc=so_,R=r)); P(f'  win {a},{b} r={rad}: HCN {hh:.3f}({hh/sh_:.1f}s) HCO+ {co:.3f}({co/so_:.1f}s) R={r:.2f}+-{r*np.hypot(sh_/hh,so_/co) if r==r else np.nan:.2f}')
for dx,dy in ((0.2,0),(-0.2,0),(0,0.2),(0,-0.2),(0.4,0),(-0.4,0),(0,0.4),(0,-0.4)):
    xx,yy=x+dx/ps,y+dy/ps; hh=apsum(wmap(H,*W0),xx,yy,RAD); co=apsum(wmap(O,*W0),xx,yy,RAD)
    P(f'  centre shift ({dx:+.1f},{dy:+.1f})": HCN {hh:.3f}({hh/sH:.1f}s) HCO+ {co:.3f}({co/sO:.1f}s) R={hh/co if co>0 and hh>0 else np.nan:.2f}')
    var.append(dict(shift=(dx,dy),h=hh,c=co,R=hh/co if co>0 and hh>0 else np.nan))
res['variants']=var
# spectra
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(2,2,figsize=(11,6),sharex=True); lines=[]
for j,(nm,R) in enumerate((('HCN 3-2',H),('HCO+ 3-2',O))):
    m=np.hypot(X-x,Y-y)*ps<=RAD; asp=R['c'][:,m].sum(1)/ba
    d=np.hypot(X-x,Y-y)*ps<=1.0; sub=np.where(d,wmap(R,*W0),-1e9); iy,ix=np.unravel_index(np.argmax(sub),sub.shape); pk=R['c'][:,iy,ix]
    ms=np.hypot(X-BP['S'][0],Y-BP['S'][1])*ps<=RAD; sp=R['c'][:,ms].sum(1)/ba
    vv=R['v']; ax[0,j].step(vv,asp*1e3,where='mid',label='Sblue ap r=0.75"'); ax[0,j].step(vv,sp*0.1*1e3,where='mid',alpha=.5,label='S nucleus ap x0.1')
    ax[0,j].axvspan(-400,-150,color='y',alpha=.15); ax[0,j].set_title(nm+' aperture spectrum [mJy km/s... Jy/beam x1e3]'); ax[0,j].set_xlim(-600,400); ax[0,j].legend(fontsize=7)
    ax[1,j].step(vv,pk*1e3,where='mid',color='C2'); ax[1,j].axvspan(-400,-150,color='y',alpha=.15); ax[1,j].set_title(f'{nm} peak-pixel spectrum (blue-wing peak within 1")  [mJy/beam]'); ax[1,j].set_xlabel('v-VSYS [km/s]')
    ax[1,j].axhline(0,color='k',lw=.5)
        # binned (3 ch) wing profile
    k3=3; lines.append((nm,[(float(vv[i:i+k3].mean()),float(asp[i:i+k3].mean()*1e3),float(pk[i:i+k3].mean()*1e3)) for i in range(0,len(vv)-k3,k3) if -600<vv[i]<400]))
plt.tight_layout(); plt.savefig('/workspace/out/s38_spectra.png',dpi=90)
with open('/workspace/out/s38_spectra.txt','w') as f:
    for nm,L in lines:
        f.write(f'# {nm}: v-VSYS[km/s], aperture(r=0.75") spectrum [1e-3 Jy km/s per beam per... mean of 3 ch, Jy/beam*1e3 summed/beamarea], peak pixel [mJy/beam]\n')
        for t in L: f.write('%7.1f %8.2f %8.2f\n'%t)
# wing kinematic structure: wing flux in 50 km/s bins
P('## T3b blue-wing profile in 50 km/s bins at Sblue aperture (r=0.75"), Jy km/s per beam: HCN | HCO+')
for a in range(-450,-100,50):
    P(f'  {a:5d}..{a+50:5d}  {apsum(wmap(H,a,a+50),x,y,RAD):+.3f} | {apsum(wmap(O,a,a+50),x,y,RAD):+.3f}   (1sigma ~ {np.array(noise(a,a+50,RAD))[:,0].std():.3f} | {np.array(noise(a,a+50,RAD))[:,1].std():.3f})')
# ---- T4 replication: Xb4b (0.25") convolved to the Xb4d beam (s11 maps, window -450..-150) -> compare with Xb4d same window
fnb=files('*Xb4b.ngc_3256_sci.spw25.cube.I.pbcor.fits')[0]; wcb=WCS(fits.getheader(fnb)).celestial; psb=abs(fits.getheader(fnb)['CDELT1'])*3600
T4={}
for mol,Rr in (('HCN',H),('HCO+',O)):
    z=np.load(f'/workspace/out/s11_Xb4b_{mol}_smoothed_maps.npz'); M=z['blue']; nb=M.shape[0]; Yb,Xb=np.mgrid[:nb,:nb]
    px,py=wcb.world_to_pixel_values(*SB); px-=int(z['x0']); py-=int(z['y0']); bab=1.1331*0.94*0.74/psb**2
    f=lambda x_,y_,rad=RAD:M[np.hypot(Xb-x_,Yb-y_)*psb<=rad].sum()/bab
    # random positions in the smoothed map, away from bright blue-wing/core emission: use central region r 4-11" from midpoint and |M| below 3x robust sigma at centre
    cxb,cyb=wcb.world_to_pixel_values(*CEN); cxb-=int(z['x0']); cyb-=int(z['y0'])
    zr=[]; tries=0
    ex=[wcb.world_to_pixel_values(*pp_) for pp_ in (REGIONS['N'][:2],REGIONS['S'][:2],SB,REGIONS['Bridge'][:2])]; ex=[(e0-int(z['x0']),e1-int(z['y0'])) for e0,e1 in ex]
    while len(zr)<300 and tries<50000:
        tries+=1; x_,y_=rng.uniform(0,nb,2); d=np.hypot(x_-cxb,y_-cyb)*psb
        if 2<d<11 and all(np.hypot(x_-e0,y_-e1)*psb>2.5 for e0,e1 in ex): zr.append(f(x_,y_))
    zr=np.array(zr); T4[mol]=(f(px,py),1.4826*np.median(abs(zr-np.median(zr))),len(zr),float(np.median(zr)))
hh,sh_=T4['HCN'][:2]; co,so_=T4['HCO+'][:2]
P(f'## T4 Xb4b convolved to Xb4d beam, window -450..-150, Sblue ap r=0.75": HCN {hh:.3f}+-{sh_:.3f} HCO+ {co:.3f}+-{so_:.3f} R={hh/co:.2f}+-{hh/co*np.hypot(sh_/hh,so_/co):.2f} (Xb4d same window: HCN 0.485 HCO+ 0.332 R=1.46); noise: MAD-sigma of {T4["HCN"][2]} aperture sums at random positions 2-11" from midpoint (>2.5" from N,S,Sblue,Bridge); medians {T4["HCN"][3]:.3f},{T4["HCO+"][3]:.3f}')
res['T4']=T4
# Xb4f native 5.7x3.5" integrated wing (blended with S core)
o4f=S.run('Xb4f',3,None,windows=(350,)); Hf,Of=o4f['HCN'],o4f['HCO+']; nf=Hf['c'].shape[1]; Yf,Xf=np.mgrid[:nf,:nf]
def apf(R,a,b,pos,rad=2.5):
    m=(R['v']>=a)&(R['v']<=b); mp=R['c'][m].sum(0)*R['dv']; px,py=R['pxy'](*pos); return mp[np.hypot(Xf-px,Yf-py)*R['ps']<=rad].sum()/R['beamarea'], mp
P('## T4b Xb4f native (5.7x3.5") wing -400..-150, aperture r=2.5" (blend with S core): HCN | HCO+ | R ; noise from the same-window scatter at 200 random positions r 10-18"')
for nm,pos in (('Sblue',SB),('S',REGIONS['S'][:2]),('N',REGIONS['N'][:2])):
    out_=[]
    for R in (Hf,Of):
        val,mp=apf(R,-400,-150,pos); cxf,cyf=R['pxy'](*CEN); zz=[]
        while len(zz)<200:
            x_,y_=rng.uniform(0,nf,2)
            if 10<np.hypot(x_-cxf,y_-cyf)*R['ps']<18: zz.append(mp[np.hypot(Xf-x_,Yf-y_)*R['ps']<=2.5].sum()/R['beamarea'])
        out_.append((val,np.std(zz)))
    (h_,eh_),(c_,ec_)=out_; P(f'  {nm:6s} HCN {h_:.3f}+-{eh_:.3f} HCO+ {c_:.3f}+-{ec_:.3f}  R={h_/c_ if c_>0 and h_>0 else np.nan:.2f}')
    res[f'Xb4f_{nm}']=dict(h=h_,eh=eh_,c=c_,ec=ec_)
# ---- T5 neighbouring lines (rest freq) in the wing windows
import csv
P('## T5 catalogue lines (out/linelist.csv, rest GHz) falling in the blue-wing window -400..-150 km/s (and 0..+/-150 core for the same line, radio def about VSYS 2780) of HCN 3-2 and HCO+ 3-2, with E_up<500 K')
L=list(csv.DictReader(open('/workspace/out/linelist.csv')))
for nm,f0 in (('HCN 3-2',265.88618),('HCO+ 3-2',267.55763)):
    # observed freq of line at velocity u (radio): nu=f0*(1-(VSYS+u)/c) ; a catalogue line at rest f is seen at nu_obs=f/(1+z), z from optical def.
    zopt=VSYS/C/(1-VSYS/C)  # radio->optical (small diff)
    lo=f0*(1-(VSYS-150)/C); hi=f0*(1-(VSYS-400)/C)  # observed freq range of the wing window (blue = higher freq)
    lo,hi=sorted((f0*(1-(VSYS-150)/C),f0*(1-(VSYS-400)/C)))
    hits=[l for l in L if l['freq_GHz'] and lo<=float(l['freq_GHz'])/(1+zopt)<=hi and float(l['upper_state_energy_K'] or 1e9)<500]
    P(f'  {nm}: obs freq window {lo:.4f}-{hi:.4f} GHz (rest-frame-in-galaxy shift z={zopt:.5f}); {len(hits)} catalogue lines E_up<500K')
    for l in hits[:12]: P(f'     {l["name"]} {l["resolved_QNs"][:30]} f_rest(gal frame)={float(l["freq_GHz"])/(1+zopt):.4f} Eup {float(l["upper_state_energy_K"]):.0f}K logA {l["aij"]}')
# ---- T6 peak-pixel-per-line bias (p0001 recipe) by injection into noise patches
P('## T6 peak-pixel-per-line recipe bias: inject beam-sized Gaussian blob of true R into real noise patches (random positions), take max pixel within 1" of centre for each line separately')
bm=lambda amp,x0_,y0_: amp*np.exp(-4*np.log(2)*(((X-x0_)*ps/0.84)**2+((Y-y0_)*ps/0.84)**2))
wH,wO=wmap(H,*W0),wmap(O,*W0)
for h0,c0 in ((0.449,0.320),(0.25,0.18),(0.15,0.15)):
    # amplitudes are aperture sums r=.75; convert to peak: peak ~ aperture sum*beam/(aperture integral)
    k_=bm(1,0,0); 
    Rr=[]; 
    for (x_,y_) in RND[:300]:
        yy=np.hypot(X-x_,Y-y_)*ps<=1.0
        th=bm(1.0,x_,y_); ap_=th[np.hypot(X-x_,Y-y_)*ps<=RAD].sum()/ba; a_h=h0/ap_; a_o=c0/ap_
        mh_=(wH+a_h*th)[yy].max(); mo_=(wO+a_o*th)[yy].max(); Rr.append((mh_/a_h, mo_/a_o, mh_/mo_ if mo_>0 else np.nan))
    Rr=np.array(Rr); P(f'  true wing aperture fluxes HCN {h0} HCO+ {c0} (true R={h0/c0:.2f}): peak-pixel R median {np.nanmedian(Rr[:,2]):.2f}, 16-84% {np.nanpercentile(Rr[:,2],16):.2f}-{np.nanpercentile(Rr[:,2],84):.2f}; peak/true-peak HCN {np.median(Rr[:,0]):.2f} HCO+ {np.median(Rr[:,1]):.2f}')
# ---- T7 RADEX: abundance ratio needed
LO={'hcn':{3:265.8864},'hco+':{3:267.5576}}
def radex(mol,Tk,nH,N,dV=150.):
    inp=f"/opt/lamda/{mol}.dat\nout_{mol}.txt\n170 370\n{Tk}\n1\nH2\n{nH}\n2.73\n{N}\n{dV}\n0\n"
    os.makedirs('/workspace/out/radex_q1w',exist_ok=True); subprocess.run(['radex_lvg'],input=inp,text=True,capture_output=True,cwd='/workspace/out/radex_q1w')
    for l in open('/workspace/out/radex_q1w/'+f'out_{mol}.txt'):
        p_=l.split()
        if len(p_)>=12 and p_[1]=='--' and abs(float(p_[4])-LO[mol][3])<2e-3: return float(p_[11]),float(p_[7])
    return np.nan,np.nan
P('## T7 RADEX (LVG, dv=150): K-ratio R_K(3-2) per unit abundance ratio A=N(HCN)/N(HCO+), N(HCO+)=1e14 (thin) and 1e15; A needed for R_K=0.55 (Sblue core), 0.40 (S core) and 1.42 (wing; Jy 1.40 x1.0126)')
T7=[]
os.chdir('/workspace/out')
for Tk in (20,50,100,200):
    for nH in (1e5,1e6,1e7):
        for NO in (1e14,1e15):
            Io,to=radex('hco+',Tk,nH,NO); Ih,th=radex('hcn',Tk,nH,NO)
            r1=Ih/Io; T7.append(dict(T=Tk,n=nH,N=NO,R_A1=r1,tauHCN=th,tauHCO=to))
            P(f'  Tk={Tk:3d} n={nH:.0e} N(HCO+)={NO:.0e}: R_K(A=1)={r1:.3f} tau(HCN,HCO+)=({th:.2f},{to:.2f}) -> A(core 0.55)={0.55/r1:.2f}  A(S core 0.40)={0.40/r1:.2f}  A(wing 1.42)={1.42/r1:.2f}  (wing/core A factor {1.42/0.55:.2f} if same excitation)')
res['T7']=T7
json.dump(res,open('/workspace/out/s38_sblue2.json','w'),indent=1,default=lambda x:x.tolist() if hasattr(x,'tolist') else float(x))
open('/workspace/out/s38_sblue2.txt','w').write('\n'.join(txt))
