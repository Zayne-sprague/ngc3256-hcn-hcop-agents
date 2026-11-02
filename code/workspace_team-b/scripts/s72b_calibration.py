# Q11 round 2, fixes 2(centroid, Gaussian bound) and 3: trials/selection calibration of the Soutflow blue-wing HCN excess, blue-vs-red asymmetry along the outflow axis,
# HCN/HCO+ centroid sensitivity, Xb4d vs Xb4b consistency.  run: python3 -I scripts/s72b_calibration.py -> out/q11r2/s72b.txt/json
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter, label
OUT='/workspace/out/q11r2'; rng=np.random.default_rng(73); txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); RES={}
Spos=REGIONS['S'][:2]; Npos=REGIONS['N'][:2]; SO=sout_pos(); Brp=REGIONS['Bridge'][:2]
axPA=np.degrees(np.arctan2(*(lambda d:(d[0],d[1]))(((REGIONS['Sblue'][0]-Spos[0])*np.cos(np.deg2rad(Spos[1]))*3600,(REGIONS['Sblue'][1]-Spos[1])*3600))))%360
P(f'# S outflow axis (S -> Sblue blob of Q0) PA = {axPA:.0f} deg E of N; Soutflow offset (-0.2,-1.85) is at PA {np.degrees(np.arctan2(-0.2,-1.85))%360:.0f} deg')
PR={o:prep_pair(o) for o in ('Xb4d','Xb4b')}
def noise(o,wins,nr=200,rad=1.0):
    H,O=PR[o]['HCN'],PR[o]['HCO+']; RP=random_pos(H['c'],H['v'],H['ps'],H['pxy'],nr,rad,5.0,9.5,rng,excl=[H['pxy'](*p) for p in (Spos,Npos,SO)])
    out=np.zeros((nr,len(wins),2))
    for i,xy in enumerate(RP):
        for j,R in enumerate((H,O)):
            sp=aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],xy,rad))
            for k,(a,b) in enumerate(wins): out[i,k,j]=W(R,sp,a,b)
    return out
WINS=[(-400,-100),(100,400),(-405,-205),(205,405)]
NZ={o:noise(o,WINS) for o in PR}; SG={o:NZ[o].std(0) for o in PR}     # [win, mol]
for o in PR: P(f'{o}@1.05": random-aperture sigma (Jy km/s) blue-half HCN {SG[o][0,0]:.3f} HCO+ {SG[o][0,1]:.3f}; red-half {SG[o][1,0]:.3f}/{SG[o][1,1]:.3f}; excess window HCN {SG[o][2,0]:.3f} HCO+ {SG[o][2,1]:.3f}')
# ================= 1. trials: grid of aperture centres within 3" of S
o='Xb4d'; H,O=PR[o]['HCN'],PR[o]['HCO+']
def blue_red(pos,o='Xb4d',vloc=0.0):
    H,O=PR[o]['HCN'],PR[o]['HCO+']; sh=ap(H,pos); so=ap(O,pos)
    return dict(hb=W(H,sh,vloc-400,vloc-100),cb=W(O,so,vloc-400,vloc-100),hr=W(H,sh,vloc+100,vloc+400),cr=W(O,so,vloc+100,vloc+400))
P('\n## 1. Grid of aperture centres (r=1") within 3" of S and of N (0.5" step, v_loc=0 windows blue -400..-100, red +100..+400), Xb4d@1.05"; criterion HCN blue half >=10 sigma and R_Tb(blue)>=1.3')
def grid(cen,rmax=3.0,step=0.5):
    out=[]
    for dx in np.arange(-rmax,rmax+1e-6,step):
        for dy in np.arange(-rmax,rmax+1e-6,step):
            if np.hypot(dx,dy)<=rmax: out.append((dx,dy,blue_red(offset_pos(cen,dx,dy))))
    return out
for nm,cen in (('S',Spos),('N',Npos)):
    G=grid(cen); zs=np.array([g[2]['hb']/SG[o][0,0] for g in G]); zc=np.array([g[2]['cb']/SG[o][0,1] for g in G]); R=np.array([g[2]['hb']/g[2]['cb']*KT if g[2]['cb']>0 else np.nan for g in G])
    ok=(zs>=10)&(R>=1.3); ok5=(zs>=5)&(zc>=3)&(R>=1.3); ok0=zs>=10
    P(f'  grid around {nm}: {len(G)} centres (~{np.pi*9/ (np.pi*1.0**2):.0f} independent r=1" apertures); HCN blue>=10 sigma at {ok0.sum()} centres; >=10s & R>=1.3: {ok.sum()} centres; (HCN>=5s & HCO+>=3s & R>=1.3): {ok5.sum()}; max HCN blue {zs.max():.1f} sigma')
    for g,z,r in zip(G,zs,R):
        if z>=8: P(f'      centre offset (E {g[0]:+.1f}", N {g[1]:+.1f}") from {nm}: HCN blue {g[2]["hb"]:.3f} ({z:.1f}s) HCO+ {g[2]["cb"]:.3f} ({g[2]["cb"]/SG[o][0,1]:.1f}s) R_Tb {r:.2f}; red HCN {g[2]["hr"]:+.3f} HCO+ {g[2]["cr"]:+.3f}')
    RES[f'grid_{nm}']=dict(n=len(G),n10=int(ok0.sum()),n10_R13=int(ok.sum()))

# outer (skirt) window -405..-205 / +205..+405 about v=0: criterion HCN >=8 sigma, HCO+ >=3 sigma and R_Tb>=1.5 (core ratios are 0.4-0.7)
def outer(pos,o='Xb4d'):
    H,O=PR[o]['HCN'],PR[o]['HCO+']; sh=ap(H,pos); so=ap(O,pos)
    return dict(hb=W(H,sh,-405,-205),cb=W(O,so,-405,-205),hr=W(H,sh,205,405),cr=W(O,so,205,405))
P('\n## 1b. Same grid with the OUTER windows (-405..-205 blue, +205..+405 red; >2.4 sigma_v from line centre so the line cores do not leak in); criterion: HCN>=8 sigma & HCO+>=3 sigma & R_Tb>=1.5')
for nm,cen in (('S',Spos),('N',Npos)):
    n8=0;nR=0;ngrid=0;lst=[]
    for dx in np.arange(-3,3.01,0.5):
        for dy in np.arange(-3,3.01,0.5):
            if np.hypot(dx,dy)>3: continue
            b=outer(offset_pos(cen,dx,dy)); ngrid+=1; zh=b['hb']/SG['Xb4d'][2,0]; zc=b['cb']/SG['Xb4d'][2,1]; R=b['hb']/b['cb']*KT if b['cb']>0 else np.nan
            if zh>=8: n8+=1
            if zh>=8 and zc>=3 and R>=1.5: nR+=1; lst.append((dx,dy,b['hb'],zh,b['cb'],zc,R))
    P(f'  grid around {nm}: {ngrid} centres; HCN outer-blue >=8 sigma at {n8}; also HCO+>=3s & R>=1.5 at {nR}')
    for l in lst: P(f'      (E {l[0]:+.1f}", N {l[1]:+.1f}"): HCN {l[2]:.3f} ({l[3]:.1f}s) HCO+ {l[4]:.3f} ({l[5]:.1f}s) R_Tb {l[6]:.2f}')
    RES[f'grid_outer_{nm}']=dict(n=ngrid,n8=n8,n8_R15=nR)
# field-wide: local maxima of HCN blue-half map (-400..-100) above 5 sigma and the HCO+ and ratio there
kH=(H['v']>=-400)&(H['v']<=-100); kO=(O['v']>=-400)&(O['v']<=-100)
mH=H['c'][kH].sum(0)*H['dv']; mO=O['c'][kO].sum(0)*O['dv']
# convert to per-aperture (r=1") equivalent: aperture sum/beamarea = mean map * (aperture area / beam area) -> use direct: smooth by r=1" top-hat
from scipy.ndimage import uniform_filter
npx=int(round(2*1.0/H['ps'])); aH=uniform_filter(mH,npx)*(np.pi*1.0**2/H['ps']**2)/H['ba']; aO=uniform_filter(mO,npx)*(np.pi*1.0**2/H['ps']**2)/O['ba']
sgH,sgO=SG[o][0,0],SG[o][0,1]
pk=(aH==maximum_filter(aH,size=int(1.5/H['ps'])))&(aH>5*sgH); lab,nl=label(pk); yy,xx=np.where(pk)
sx,sy=H['pxy'](*Spos); nx,ny=H['pxy'](*Npos)
P(f'\n  field-wide local maxima (r=1" aperture sums over the whole 48" field) with HCN blue-half (-400..-100) >= 5 sigma ({sgH:.3f}): '+f'{len(yy)}')
for y,x in sorted(zip(yy,xx),key=lambda t:-aH[t]):
    P(f'     pixel offset from S (E {-(x-sx)*H["ps"]:+.1f}", N {(y-sy)*H["ps"]:+.1f}"): HCN {aH[y,x]:.3f} ({aH[y,x]/sgH:.1f}s) HCO+ {aO[y,x]:.3f} ({aO[y,x]/sgO:.1f}s) R_Tb {aH[y,x]/aO[y,x]*KT if aO[y,x]>0 else np.nan:.2f}')
# ================= 2. blue/red asymmetry along the outflow axis
P(f'\n## 2. Positions along the S outflow axis (PA {axPA:.0f}); t>0 = SSW (blue-lobe side), t<0 = NNE (toward the Bridge/N); r=1" apertures, v_loc=0, Jy km/s (sigma HCN {SG[o][0,0]:.3f}, HCO+ {SG[o][0,1]:.3f}); R_Tb blue/red only if both lines >2 sigma')
AX=[]
for t in np.arange(-4.0,4.01,0.5):
    pos=offset_pos(Spos,t*np.sin(np.deg2rad(axPA)),t*np.cos(np.deg2rad(axPA))); b=blue_red(pos)
    r1=b['hb']/b['cb']*KT if b['hb']>2*sgH and b['cb']>2*sgO else np.nan; r2=b['hr']/b['cr']*KT if b['hr']>2*SG[o][1,0] and b['cr']>2*SG[o][1,1] else np.nan
    AX.append(dict(t=t,**b,Rb=r1,Rr=r2))
    bo=outer(pos); ro=bo['hb']/bo['cb']*KT if bo['hb']>2*SG[o][2,0] and bo['cb']>2*SG[o][2,1] else np.nan; rr_=bo['hr']/bo['cr']*KT if bo['hr']>2*SG[o][3,0] and bo['cr']>2*SG[o][3,1] else np.nan
    AX[-1].update(outer_blue=bo['hb'],outer_blueCO=bo['cb'],outer_red=bo['hr'],outer_redCO=bo['cr'],R_outer_blue=ro,R_outer_red=rr_)
    P(f'   t={t:+.1f}": OUTER blue HCN {bo["hb"]:+.3f} ({bo["hb"]/SG[o][2,0]:+.1f}s) HCO+ {bo["cb"]:+.3f} ({bo["cb"]/SG[o][2,1]:+.1f}s) R {ro:.2f} | OUTER red HCN {bo["hr"]:+.3f} ({bo["hr"]/SG[o][3,0]:+.1f}s) HCO+ {bo["cr"]:+.3f} ({bo["cr"]/SG[o][3,1]:+.1f}s) R {rr_:.2f}')
    P(f'   t={t:+.1f}": blue HCN {b["hb"]:+.3f} ({b["hb"]/sgH:+.1f}s) HCO+ {b["cb"]:+.3f} ({b["cb"]/sgO:+.1f}s) R {r1:.2f} | red HCN {b["hr"]:+.3f} ({b["hr"]/SG[o][1,0]:+.1f}s) HCO+ {b["cr"]:+.3f} ({b["cr"]/SG[o][1,1]:+.1f}s) R {r2:.2f}')
RES['axis']=AX
# ================= 3. centroid sensitivity & Gaussian fits at Soutflow, Xb4d and Xb4b
P('\n## 3. Centroid sensitivity of the Soutflow excess; windows relative to the chosen frame centre. k = HCN/HCO+ (Jy) in the frame core -205..+95; E=HCN(-405..-205)-k*HCO+(-405..-205)')
for o in ('Xb4d','Xb4b'):
    H,O=PR[o]['HCN'],PR[o]['HCO+']; sh=ap(H,SO); so=ap(O,SO); wn=NZ[o]
    def fitc(R,s,lim):
        f=lambda x,a,m,w:a*np.exp(-0.5*((x-m)/w)**2); k=np.abs(R['v'])<lim
        try: p,cv=curve_fit(f,R['v'][k],s[k],p0=[s[k].max(),0,70],bounds=([0,-200,15],[np.inf,200,250])); return p[1],np.sqrt(cv[1,1])
        except Exception: return np.nan,np.nan
    cH250=fitc(H,sh,250); cO250=fitc(O,so,250); cH120=fitc(H,sh,120); cO120=fitc(O,so,120)
    P(f'  {o}: Gaussian centroid |v|<250: HCN {cH250[0]:+.0f}+-{cH250[1]:.0f}, HCO+ {cO250[0]:+.0f}+-{cO250[1]:.0f}; fit range |v|<120 (skirt excluded): HCN {cH120[0]:+.0f}+-{cH120[1]:.0f}, HCO+ {cO120[0]:+.0f}+-{cO120[1]:.0f}')
    for lab_,(vh,vo) in (('v_loc=0 both',(0,0)),('HCO+ centroid both (-29)',(-29,-29)),('HCN centroid both (-70)',(-70,-70)),('own centroids (HCN -70, HCO+ -29)',(-70,-29)),('own, skirt-free (|v|<120) centroids',(cH120[0],cO120[0]))):
        hk=W(H,sh,vh-205,vh+95); ck=W(O,so,vo-205,vo+95); k=hk/ck; hw=W(H,sh,vh-405,vh-205); cw=W(O,so,vo-405,vo-205)
        E=hw-k*cw; eE=np.hypot(SG[o][2,0],k*SG[o][2,1])
        P(f'    {lab_:36s}: k={k:.3f}  HCN win {hw:+.3f} HCO+ win {cw:+.3f}  E={E:+.3f}+-{eE:.3f} ({E/eE:+.1f}s);  blue-half R_Tb (vloc-400..-100) {W(H,sh,vh-400,vh-100)/W(O,so,vo-400,vo-100)*KT:.2f}')
    # free Gaussian fit to the residual spectrum
    k0=W(H,sh,-205,95)/W(O,so,-205,95); io=np.argsort(O['v']); r=sh-k0*np.interp(H['v'],O['v'][io],so[io]); v=H['v']   # NB: cube velocity axes are descending -> sort before np.interp (s70 did not: its Gaussian fit was invalid)
    lf=(v<-600)|(v>300); rms=r[lf].std(); lf=(v<-600)|(v>300); rms=r[lf].std()
    f=lambda x,a,m,w:a*np.exp(-0.5*((x-m)/w)**2)
    for lo,hi in ((-500,-150),(-700,100)):
        m=(v>-700)&(v<-100)
        try:
            p,cv=curve_fit(f,v[m],r[m],p0=[r[m].max(),-300,100],bounds=([0,lo,10],[1,hi,500]),sigma=np.full(m.sum(),rms)); pe=np.sqrt(np.diag(cv)); chi=(((r[m]-f(v[m],*p))/rms)**2).sum()
            chi0=((r[m]/rms)**2).sum(); P(f'    {o} free-ish Gaussian to residual (centre bound {lo}..{hi}, fit -700..-100, chan rms {rms*1e3:.2f} mJy): peak {p[0]*1e3:.2f}+-{pe[0]*1e3:.2f} mJy centre {p[1]:.0f}+-{pe[1]:.0f} FWHM {2.355*p[2]:.0f}+-{2.355*pe[2]:.0f}; chi2 {chi:.0f} vs null {chi0:.0f} for {m.sum()} chan (channels are correlated: smoothing, ~1.2 chan/ resolution element)')
        except Exception as ex: P('    fit failed '+str(ex))

# R(v) in 50 km/s bins at Soutflow (Xb4d), both lines, v_loc=0: continuity of the HCN/HCO+ spectra from line core into the skirt
o='Xb4d'; H,O=PR[o]['HCN'],PR[o]['HCO+']; sh=ap(H,SO); so=ap(O,SO)
P('\n## 3b. Soutflow Xb4d@1.05" spectra in 50 km/s bins (Jy km/s per bin; sigma per bin from random apertures), v_loc=0')
nzb=[]
rngb=np.random.default_rng(5); RPb=random_pos(H['c'],H['v'],H['ps'],H['pxy'],100,1.0,5.0,9.5,rngb,excl=[H['pxy'](*p_) for p_ in (Spos,Npos,SO)])
bins=[(a,a+50) for a in range(-550,150,50)]
rb=np.array([[ [W(R,aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],xy,1.0)),a,b) for (a,b) in bins] for R in (H,O)] for xy in RPb]); sgb=rb.std(0)
RV=[]
for j,(a,b) in enumerate(bins):
    h,c=W(H,sh,a,b),W(O,so,a,b); r_=h/c*KT if (h>2*sgb[0,j] and c>2*sgb[1,j]) else np.nan
    RV.append(dict(v=0.5*(a+b),HCN=h,HCO=c,eHCN=sgb[0,j],eHCO=sgb[1,j],R=r_))
    P(f'   {0.5*(a+b):+5.0f}: HCN {h:+.3f} ({h/sgb[0,j]:+.1f}s)  HCO+ {c:+.3f} ({c/sgb[1,j]:+.1f}s)  R_Tb {r_:.2f}')
RES['Rv']=RV
# ================= 4. Xb4d vs Xb4b at Soutflow (same recipe) and Bridge red
P('\n## 4. Xb4d vs Xb4b (both 1.05") Jy km/s, v_loc=0 windows')
cmp={}
for nm,pos in (('Soutflow',SO),('Bridge',Brp),('S',Spos),('N',Npos)):
    row={}
    for o in ('Xb4d','Xb4b'):
        H,O=PR[o]['HCN'],PR[o]['HCO+']; sh=ap(H,pos); so=ap(O,pos)
        row[o]={'HCN_blue':W(H,sh,-405,-205),'HCO_blue':W(O,so,-405,-205),'HCN_red':W(H,sh,205,405),'HCO_red':W(O,so,205,405),'HCN_b400':W(H,sh,-400,-100),'HCO_b400':W(O,so,-400,-100),'HCN_r400':W(H,sh,100,400),'HCO_r400':W(O,so,100,400)}
    cmp[nm]=row
    for key,wi,mol in (('HCN_blue',2,0),('HCO_blue',2,1),('HCN_b400',0,0),('HCO_b400',0,1),('HCN_r400',1,0),('HCO_r400',1,1)):
        a,b=row['Xb4d'][key],row['Xb4b'][key]; ea,eb=SG['Xb4d'][wi,mol],SG['Xb4b'][wi,mol]
        P(f'  {nm:8s} {key:9s}: Xb4d {a:+.3f}+-{ea:.3f}  Xb4b {b:+.3f}+-{eb:.3f}  diff(Xb4b-Xb4d) {b-a:+.3f}+-{np.hypot(ea,eb):.3f} ({(b-a)/np.hypot(ea,eb):+.1f}s)  Xb4b/Xb4d {b/a if a>0 else np.nan:.2f}')
RES['cmp']=cmp
# Xb4b ring at 1.85" from S: E blue per azimuth
o='Xb4b'; H,O=PR[o]['HCN'],PR[o]['HCO+']; P('  Xb4b ring 1.85" from S (E=HCN-0.534*HCO+ in -405..-205; sigma %.3f): '%np.hypot(SG[o][2,0],0.534*SG[o][2,1]))
for az in range(0,360,45):
    pos=offset_pos(Spos,1.85*np.sin(np.deg2rad(az)),1.85*np.cos(np.deg2rad(az))); h=W(H,ap(H,pos),-405,-205); c=W(O,ap(O,pos),-405,-205)
    P(f'     az{az:03d}: HCN {h:+.3f} HCO+ {c:+.3f}  E {h-0.534*c:+.3f}')
json.dump(RES,open(OUT+'/s72b.json','w'),indent=1,default=float); open(OUT+'/s72b.txt','w').write('\n'.join(txt))
