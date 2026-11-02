# Q11 round 2, fix 2: identity test of the 'blend' alternative for the Soutflow blue skirt (HCN frame -405..-205 km/s) with companion lines.
# Candidate carriers (Q0 line list out/linelist.csv, JPL/CDMS): HDO, CH3CHO, c-C3H, l-C3H (v=0, v4=1). For each species and T (60,100,200 K, optically thin LTE, partition function cancels):
#   F_i(Jy km/s in the same aperture) ∝ g_u A e^{-E_u/T} (Rayleigh-Jeans: W_K ∝ g A e^{-E/T}/nu^2, Jy ∝ W nu^2).  Template = sum of Gaussians (FWHM assumed 195 or 100 km/s) at v_loc offset;
#   normalised so that the template integral over the HCN-frame window -405..-205 equals the observed excess E (HCN - k*HCO+; Xb4d Soutflow 0.262 Jy km/s at v_loc=0).
#   Matched filter on the Soutflow aperture spectrum of each of the four Xb4d spws (spw25, 27(HCN), 29, 31): s = sum(d t)/sum(t^2); s=1 expected if the blend explains E; noise sigma_s from 200 random line-free apertures.
# Reference lines of a species = its catalogue lines inside the excess window (HCN frame -405..-205 km/s, +-40 km/s margin, centred at v_loc); companions = all other lines in the four spws except those within 400 km/s of the window centre (contaminated by HCN/the skirt itself).
# Spectra are high-pass filtered (running median, 600 km/s) before matching so that local baselines are used. Output f_max = 95% upper limit on the fraction of E that the species can carry.
# run: python3 -I scripts/s72c_blend_identity.py -> out/q11r2/s72c.txt/json
import sys, json, csv; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
OUT='/workspace/out/q11r2'; rng=np.random.default_rng(74); txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); RES={}
SO=sout_pos(); H0=6.62607015e-27; CC=2.99792458e10
def gu_(nu_ghz,smu2,A): nu=nu_ghz*1e9; return 64*np.pi**4*nu**3*smu2*1e-36/(3*H0*CC**3*A)
# ---- line list
SPECIES={'HDO':lambda l:l['sp']=='HDO','CH3CHO':lambda l:l['sp']=='CH3CHO','c-C3H':lambda l:l['sp']=='c-C3H','l-C3H':lambda l:l['sp']=='l-C3H'}
L=[]
for l in csv.DictReader(open('/workspace/out/linelist.csv')):
    try: f=float(l['freq_GHz']); E=float(l['upper_state_energy_K']); A=10**float(l['aij']); sm=float(l['sijmu2'])
    except: continue
    for sp,fn in SPECIES.items():
        if fn(l) and E<900 and sm>0: L.append(dict(sp=sp,f=f,Eu=E,A=A,g=gu_(f,sm,A),qn=l['resolved_QNs'].replace('\n',' ')[:28],nm=l['name'][:28]))
P(f'# {len(L)} catalogue lines (Eu<900 K) of {list(SPECIES)}')
# ---- cubes: all channels, 1.05" circular beam, per-pixel median baseline over the whole spw
SPW=('spw25','spw27','spw29','spw31'); CUB={}
for sp in SPW:
    c,v,ps,wc,x0,y0,b0=load_cube('Xb4d',sp,F0['HCN'][3],vmin=-1e5,vmax=1e5)
    nu=F0['HCN'][3]*(1-(v+VSYS)/C); c=conv(c,ps,b0,(1.05,1.05,0.0)); c=c-np.median(c,axis=0)
    CUB[sp]=dict(c=c,nu=nu,dv_hz=abs(np.median(np.diff(nu))),ps=ps,n=c.shape[1],ba=1.1331*1.05**2/ps**2,pxy=pxy_fn(wc,x0,y0),vhcn=v)
    P(f'  {sp}: obs {nu.min()/1e9:.3f}-{nu.max()/1e9:.3f} GHz -> rest (VSYS) {nu.min()/1e9/(1-VSYS/C):.3f}-{nu.max()/1e9/(1-VSYS/C):.3f}; shape {c.shape}')
R27=CUB['spw27']; Spos=REGIONS['S'][:2]; Npos=REGIONS['N'][:2]
excl=[R27['pxy'](*p) for p in (Spos,Npos,SO)]
RP=random_pos(R27['c'],R27['vhcn'],R27['ps'],R27['pxy'],200,1.0,5.0,9.5,rng,excl=excl)
def spec_ap(R,xy): return aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],xy,1.0))
for sp in SPW:
    R=CUB[sp]; R['d']=spec_ap(R,R['pxy'](*SO)); R['rand']=np.array([spec_ap(R,xy) for xy in RP])
    R['rms']=np.median(R['rand'].std(0))
    P(f'  {sp}: Soutflow aperture spectrum rms (random apertures, per channel, median) {R["rms"]*1e3:.2f} mJy; channel width {R["dv_hz"]/1e6:.3f} MHz')
for sp in SPW:   # baseline-A-consistent version of the aperture spectra for the E measurement (median of -900..-450 km/s in the HCN frame / own frame)
    R=CUB[sp]; f0=F0['HCO+'][3] if sp=='spw25' else F0['HCN'][3]; uu=C*(1-R['nu']/f0)-VSYS; kb=(uu>-900)&(uu<-450)
    R['dA']=R['d']-np.median(R['d'][kb]); R['randA']=R['rand']-np.median(R['rand'][:,kb],axis=1)[:,None]
# observed excess E (v_loc=0 frame) and k from s70/s72b: recompute here with spw27/spw25
H=R27; O=CUB['spw25']
def win(R,d,a,b,f0=F0['HCN'][3]):
    u=C*(1-R['nu']/f0)-VSYS; k=(u>=a)&(u<=b); return d[k].sum()*abs(np.median(np.diff(u)))
def uvec(R,f0): return C*(1-R['nu']/f0)-VSYS
# HCO+ spectrum on same velocity frame for k
hk=win(H,H['dA'],-205,95); ck=win(O,O['dA'],-205,95,F0['HCO+'][3]); k=hk/ck
hw=win(H,H['dA'],-405,-205); cw=win(O,O['dA'],-405,-205,F0['HCO+'][3]); E=hw-k*cw
sgE=np.std([win(H,r,-405,-205)-k*win(O,rr,-405,-205,F0['HCO+'][3]) for r,rr in zip(H['randA'],O['randA'])])
P(f'\nSoutflow Xb4d@1.05" v_loc=0: k={k:.3f}; HCN window {hw:.3f}, HCO+ window {cw:.3f} -> E={E:.3f}+-{sgE:.3f} Jy km/s (random-aperture sigma, {len(RP)} apertures)')
RES['E']=dict(k=k,HCN=hw,HCO=cw,E=E,sigE=sgE)
# ---- templates
from scipy.ndimage import median_filter
def hp(x,nch=137):
    return x-median_filter(x,size=nch,mode='nearest')
for sp in SPW:
    R=CUB[sp]; R['dhp']=hp(R['d']); R['randhp']=np.array([hp(r) for r in R['rand']])
def template(R,lines,scale,vloc,fwhm):
    nu=R['nu']; t=np.zeros_like(nu); sg=fwhm/2.3548
    for l in lines:
        nui=l['f']*1e9*(1-(VSYS+vloc)/C); dvv=C*(nui-nu)/nui; t+=scale*l['F']*np.exp(-0.5*(dvv/sg)**2)/(sg*np.sqrt(2*np.pi))
    return t        # Jy per channel
def uH_of(l,vloc): return C*(1-l['f']*1e9*(1-(VSYS+vloc)/C)/F0['HCN'][3])-VSYS
def run(sp_name,T,vloc,fwhm):
    ls=[dict(l) for l in L if l['sp']==sp_name]
    for l in ls: l['F']=l['g']*l['A']*np.exp(-l['Eu']/T)
    lo,hi=H['nu'].min(),H['nu'].max()
    ref=[l for l in ls if lo<l['f']*1e9*(1-(VSYS+vloc)/C)<hi and -445<uH_of(l,vloc)-vloc<-165]    # window -405..-205 +-40 (frame of the line centroid)
    if not ref: return None
    # contributions to the window: fraction of each Gaussian inside -405..-205 (HCN frame) -> scale so that template window flux = E
    t27=template(H,ref,1.0,vloc,fwhm); u=uvec(H,F0['HCN'][3]); kw=(u>=-405)&(u<=-205); w27=t27[kw].sum()*abs(np.median(np.diff(u)))
    scale=E/w27
    refset={id(l) for l in ref}
    comp=[l for l in ls if id(l) not in refset]
    out=dict(sp=sp_name,T=T,vloc=vloc,fwhm=fwhm,n_ref=len(ref),ref=[(l['f'],l['qn'],l['Eu'],round(uH_of(l,vloc))) for l in ref],scale=scale,mf={})
    res=[]
    for sp in SPW:
        R=CUB[sp]; 
        cl=[l for l in comp if R['nu'].min()<l['f']*1e9*(1-(VSYS+vloc)/C)<R['nu'].max()]
        if sp=='spw27': cl=[l for l in cl if abs(uH_of(l,vloc)+305)>400]
        cl=[l for l in cl if l['F']*scale>0.02*E]            # only lines predicted >= 2% of E matter
        if not cl: out['mf'][sp]=None; continue
        t=hp(template(R,cl,scale,vloc,fwhm)); n2=(t**2).sum()
        if n2<=0: out['mf'][sp]=None; continue
        s_=(R['dhp']*t).sum()/n2; sr=(R['randhp']@t)/n2; sig=sr.std()
        out['mf'][sp]=dict(s=float(s_),sig=float(sig),nlines=len(cl),maxpred_Jykms=float(max(l['F'] for l in cl)*scale))
        res.append((s_,sig))
    if res:
        w=np.array([1/r[1]**2 for r in res]); sv=np.array([r[0] for r in res]); sc=(w*sv).sum()/w.sum(); sg=w.sum()**-0.5
        out['comp']=dict(s=float(sc),sig=float(sg),fmax95=float(max(sc+2*sg,0)))
    return out
RUN=[]
for sp_name in SPECIES:
    for T in (60,100,200,500):
        for vloc in (0,-29):
            for fwhm in (100,195):
                o=run(sp_name,T,vloc,fwhm)
                if o: RUN.append(o)
RES['runs']=RUN
P('\n## Species lines falling in the excess window (T-independent list; vloc=-29)')
seen=set()
for o in RUN:
    if o['vloc']==-29 and o['sp'] not in seen: seen.add(o['sp']); P(f'  {o["sp"]:7s}: '+'; '.join(f'{f:.4f} GHz {q} Eu={e:.0f}K u={u:+d}' for f,q,e,u in o['ref']))
P('\n## Companion test: s = (measured)/(predicted if the species alone carries E), stacked over spw25/27/29/31 (spw27 excludes lines within 400 km/s of window), sigma from 200 random apertures; f_max = s+2sigma = 95% upper limit on the fraction of E the species can carry')
for o in RUN:
    c=o.get('comp'); mf=o['mf']; f=lambda sp:(f"{mf[sp]['s']:+.3f}+-{mf[sp]['sig']:.3f}(n={mf[sp]['nlines']})" if mf[sp] else 'n/a')
    if c: P(f'  {o["sp"]:7s} T={o["T"]:3d} vloc={o["vloc"]:+3d} FWHM={o["fwhm"]:3d}: spw25 {f("spw25")} spw27 {f("spw27")} spw29 {f("spw29")} spw31 {f("spw31")} | stack s={c["s"]:+.3f}+-{c["sig"]:.3f}  f_max(95%)={c["fmax95"]:.3f}')
    else: P(f'  {o["sp"]:7s} T={o["T"]:3d} vloc={o["vloc"]:+3d} FWHM={o["fwhm"]:3d}: no companion >=2% of E in the four spws -> untestable')
# ---- strongest companion lines (T=100, vloc=-29, FWHM=195), outside host contamination zone
P('\n## Strongest individual companion lines (T=100 K, v_loc=-29, FWHM 195; species alone carries E): predicted vs measured Jy km/s within +-100 km/s at Soutflow, high-pass baseline; sigma from random apertures')
TAB=[]
for sp_name in SPECIES:
    o=run(sp_name,100,-29,195)
    if not o: continue
    ls=[l for l in L if l['sp']==sp_name]
    for l in ls: l['F']=l['g']*l['A']*np.exp(-l['Eu']/100)
    refs={round(r[0],4) for r in o['ref']}
    cand=[]
    for l in ls:
        if round(l['f'],4) in refs: continue
        nui=l['f']*1e9*(1-(VSYS-29)/C)
        for sp in SPW:
            R=CUB[sp]
            if nui<R['nu'].min() or nui>R['nu'].max(): continue
            if sp=='spw27' and abs(uH_of(l,-29)+305)<=400: continue
            if l['F']*o['scale']>0.02*E: cand.append((l['F']*o['scale'],l,sp,nui))
    cand.sort(key=lambda t:-t[0]); done=[]
    for Pt,l,sp,nui in cand:
        if any(abs(l['f']-d)<0.02 for d in done): continue      # merge hyperfine components
        done.append(l['f']); R=CUB[sp]; dvv=C*(nui-R['nu'])/nui; kk=np.abs(dvv)<=100; dvch=abs(np.median(np.diff(dvv)))
        # predicted flux in window includes all components within 0.02 GHz: sum
        Ptot=sum(m['F'] for m in ls if abs(m['f']-l['f'])<0.02)*o['scale']; pred=Ptot*0.77
        meas=R['dhp'][kk].sum()*dvch; sg=(R['randhp'][:,kk].sum(1)*dvch).std()
        TAB.append(dict(sp=sp_name,f=l['f'],qn=l['qn'],Eu=l['Eu'],spw=sp,pred_win=float(pred),meas=float(meas),sig=float(sg)))
        P(f'  {sp_name:7s} {l["f"]:.4f} GHz {l["qn"]:26s} Eu={l["Eu"]:4.0f}K ({sp}): predicted {pred:9.3f}  measured {meas:+.3f}+-{sg:.3f} ({meas/sg:+.1f}s)  -> species can carry <= {max(meas+2*sg,0)/pred*100:.2f}% of E')
        if len(done)>=4: break

# ---- shape test in the HCN spw: residual r = HCN - k*HCO+ (baseline A) vs (i) 3-parameter Gaussian, (ii) reference-line templates of each species (amplitude fixed by E; and amplitude free)
from scipy.optimize import curve_fit
vv=uvec(H,F0['HCN'][3]); uo=uvec(O,F0['HCO+'][3]); io=np.argsort(uo); r=H['dA']-k*np.interp(vv,uo[io],O['dA'][io]); m=(vv>-700)&(vv<-100); rms=H['rms']
gf=lambda x,a,mu,w:a*np.exp(-0.5*((x-mu)/w)**2)
pg,_=curve_fit(gf,vv[m],r[m],p0=[0.003,-300,100],bounds=([0,-700,10],[1,100,500])); chi_g=(((r[m]-gf(vv[m],*pg))/rms)**2).sum(); chi_0=((r[m]/rms)**2).sum()
P(f'\n## Shape test, HCN spw27, residual HCN-{k:.3f}*HCO+ (channel rms {rms*1e3:.2f} mJy; {m.sum()} channels -700..-100; channels are correlated so chi2 are indicative): null chi2 {chi_0:.0f}; free Gaussian (peak {pg[0]*1e3:.2f} mJy, centre {pg[1]:.0f}, FWHM {2.355*pg[2]:.0f}) chi2 {chi_g:.0f}')
SH=[]
for o in RUN:
    if o['T']!=100: continue
    ls=[dict(l) for l in L if l['sp']==o['sp']]
    for l in ls: l['F']=l['g']*l['A']*np.exp(-l['Eu']/100)
    refs={round(q[0],4) for q in o['ref']}; rl=[l for l in ls if round(l['f'],4) in refs]
    t=template(H,rl,o['scale'],o['vloc'],o['fwhm']); chi_t=(((r[m]-t[m])/rms)**2).sum(); a=(r[m]*t[m]).sum()/(t[m]**2).sum(); chi_a=(((r[m]-a*t[m])/rms)**2).sum()
    SH.append(dict(sp=o['sp'],vloc=o['vloc'],fwhm=o['fwhm'],chi2_fixed=float(chi_t),chi2_free_amp=float(chi_a),amp=float(a)))
    P(f'  {o["sp"]:7s} vloc={o["vloc"]:+3d} FWHM={o["fwhm"]:3d}: template chi2 {chi_t:.0f} (amplitude fixed to E), best amplitude x{a:.2f} chi2 {chi_a:.0f}')
RES['shape']=dict(chi2_null=float(chi_0),chi2_gauss=float(chi_g),gauss=pg.tolist(),templates=SH)
RES['table']=TAB
json.dump(RES,open(OUT+'/s72c.json','w'),indent=1,default=float); open(OUT+'/s72c.txt','w').write('\n'.join(txt))
