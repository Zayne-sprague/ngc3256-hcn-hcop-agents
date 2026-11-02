# Q11 parts 3-4.  (3) identity test of the HCN excess over k*HCO+ at HCN-frame 2375-2575 km/s (= -405..-205 km/s about VSYS 2780) team-a p0006: (a) reproduce, (b) residual-spectrum shape and
# velocity (blue skirt vs fixed-frequency line), (c) does the excess exist at N/Bridge/S (a position-independent blend would), (d) catalogue lines in the window (out/linelist.csv, Q0) and the LTE column
# density a blend would need, relative to the HCN column the same aperture implies.  (4) N wing limits (3 sigma from random apertures).   Baseline: per-pixel median of -900<v<-450 (part 3) .
# Output: out/q11/s70.txt/.json   run: python3 -I scripts/s70_q11_blend_nlimit.py
import sys, json, os, csv; sys.path.insert(0,'/workspace/scripts')
os.environ.setdefault('Q11_BASE','A')
src=open('/workspace/scripts/s69_q11_wing_check.py').read().split('RES={}')[0]
exec(src)
OUT='/workspace/out/q11'
RES={}
H0=6.62607015e-27; CC=2.99792458e10; KB=1.380649e-16
def Wk(Jy,nu,rad_as,beam=None):
    """integrated Jy km/s over aperture of radius rad_as (source assumed to fill it) -> K km/s"""
    om=np.pi*(rad_as*np.pi/180/3600)**2; return Jy*1e-23*(CC**2)/(2*KB*(nu*1e9)**2*om)
def Nu_from_W(W,nu_ghz,A):   # W in K km/s, optically thin: N_u = 8 pi k nu^2 W /(h c^3 A)   (cgs, W in K cm/s)
    nu=nu_ghz*1e9; return 8*np.pi*KB*nu**2*(W*1e5)/(H0*CC**3*A)
def gu(nu_ghz,smu2,A):
    nu=nu_ghz*1e9; return 64*np.pi**4*nu**3*smu2*1e-36/(3*H0*CC**3*A)
# ---------------------------------------------------------------- Part 3
P('# Part 3: blend/identity test of HCN excess over k*HCO+ ; velocities relative to VSYS=2780 km/s (HCN-frame 2375-2575 = -405..-205; k fitted on 2575-2875 = -205..+95), aperture r=1.0", Xb4d at 1.05" and Xb4b at 1.05".')
P3={}
BIN=[(a,a+50) for a in range(-600,300,50)]
for ous in ('Xb4d','Xb4b'):
    o=prep(ous,1.05); H,O=o['HCN'],o['HCO+']
    RP=randpos(o,60,1.0,5.0,9.5,excl=[H['pxy'](*pos(p)) for p in POS])
    # noise on residual window: need paired noise for E=h-k*c per random aperture -> compute per-aperture window sums
    wl=[(-205,95),(-405,-205)]+BIN
    nz=noise_vec(o,RP,1.0,wl)
    P3[ous]={}
    P(f'## {ous}@1.05"')
    for pn in POS:
        xy=H['pxy'](*pos(pn)); sH=spec(H,apmask(H,xy,1.0)); sO=spec(O,apmask(O,xy,1.0))
        hk,_=win(H,sH,-205,95); ck,_=win(O,sO,-205,95); k=hk/ck
        hw,_=win(H,sH,-405,-205); cw,_=win(O,sO,-405,-205)
        E=hw-k*cw; eE=(nz[:,1,0]-k*nz[:,1,1]).std()      # scatter at random apertures of the same combination (k fixed)
        # also treat k uncertainty
        ek=k*np.hypot(nz[:,0,0].std()/hk,nz[:,0,1].std()/ck); eE_tot=np.hypot(eE,cw*ek)
        hc,_=win(H,sH,-100,100)
        P3[ous][pn]=dict(k=k,k_Tb=k*KT,HCN_w=hw,HCO_w=cw,E=E,eE=eE_tot,E_over_HCNw=E/hw if hw>0 else np.nan,HCN_core=hc)
        P(f'  {pn:8s} k={k:.3f} (R_Tb {k*KT:.3f}); window -405..-205: HCN {hw:+.3f} HCO+ {cw:+.3f}  excess E=HCN-k*HCO+ = {E:+.3f}+-{eE_tot:.3f} Jy km/s ({E/eE_tot:+.1f} sigma; {100*E/hw if hw>0 else float("nan"):.0f}% of HCN window flux; E/HCN_core(|v|<100)={E/hc:+.4f})')
    # residual spectrum at Soutflow and others in 50 km/s bins
    for pn in ('Soutflow','N','S'):
        xy=H['pxy'](*pos(pn)); sH=spec(H,apmask(H,xy,1.0)); sO=spec(O,apmask(O,xy,1.0)); k=P3[ous][pn]['k']
        rows=[]
        for j,(a,b) in enumerate(BIN):
            h,_=win(H,sH,a,b); c,_=win(O,sO,a,b); e=(nz[:,2+j,0]-k*nz[:,2+j,1]).std(); rows.append((0.5*(a+b),h,c,h-k*c,e))
        P3[ous][pn]['resid_bins']=rows
        P(f'    residual HCN-k*HCO+ (Jy km/s per 50 km/s) {pn}: '+' '.join(f'{v:+.0f}:{r:+.3f}({r/e:+.0f}s)' for v,h,c,r,e in rows if -550<v<100))
    # Gaussian fit to the Soutflow residual spectrum (channel level, fit range -650..50)
    xy=H['pxy'](*pos('Soutflow')); sH=spec(H,apmask(H,xy,1.0)); sO=spec(O,apmask(O,xy,1.0)); k=P3[ous]['Soutflow']['k']
    v=H['v']; r=sH-k*np.interp(v,O['v'],sO); m=(v>-520)&(v<-120)
    f=lambda x,a,mu,w:a*np.exp(-0.5*((x-mu)/w)**2)
    try:
        p,cv=curve_fit(f,v[m],r[m],p0=[r[m].max(),-300,50],bounds=([0,-500,10],[1,-150,300])); pe=np.sqrt(np.diag(cv))
        P(f'    Gaussian (centre bounded to -500..-150, fit range -520..-120) to Soutflow residual spectrum: peak {p[0]*1e3:.2f}+-{pe[0]*1e3:.2f} mJy, centre {p[1]:.0f}+-{pe[1]:.0f}, FWHM {2.355*p[2]:.0f}+-{2.355*pe[2]:.0f} km/s (channel rms ~{np.std(r[v<-800])*1e3:.2f} mJy); HDO 2(2,0)-3(1,3) 266.161 GHz expected at -307, c-C3H -323/-355, CH3CHO 12(6,6)-13(5,8) -365')
        P3[ous]['gauss_resid']=dict(p=p.tolist(),pe=pe.tolist(),fwhm=2.355*p[2])
    except Exception as ex: P(f'    gauss fit failed {ex}')
    # fixed-frequency vs local-velocity: excess window defined relative to vloc vs fixed (Soutflow vloc ~ -30)
    del o
RES['part3']=P3
# catalogue lines in window and required LTE columns
L=list(csv.DictReader(open('/workspace/out/linelist.csv'))); zf=(1-VSYS/C)
P('\n## catalogue lines (out/linelist.csv, Q0) falling at -405..-205 km/s of the HCN 3-2 frame, Eup<400 K, logA>-6 (galaxy-frame; a line of velocity u about systemic is observed at f*(1-VSYS/c))')
cands=[]
for l in L:
    try: f=float(l['freq_GHz']); Eu=float(l['upper_state_energy_K']); A=10**float(l['aij']); smu=float(l['sijmu2'])
    except: continue
    u=C*(1-f*zf/S.F0['HCN'][3]*1e9)-VSYS
    if -420<u<-190 and Eu<400 and A>1e-6: cands.append((u,l['name'],l['resolved_QNs'].replace('\n',' ')[:34],f,Eu,A,smu))
for c in sorted(cands): P(f'   u={c[0]:+.0f} {c[1]} {c[2]} f={c[3]:.4f} Eup={c[4]:.0f}K A={c[5]:.2e}')
# required column for HDO 2(2,0)-3(1,3) and CH3CHO 12(6,6)-13(5,8) at the measured Soutflow Xb4d excess
E=P3['Xb4d']['Soutflow']['E']; eE=P3['Xb4d']['Soutflow']['eE']; HW=P3['Xb4d']['Soutflow']['HCO_w']
P(f'\n## LTE column a blend would need to carry the Soutflow Xb4d excess E={E:.3f}+-{eE:.3f} Jy km/s (aperture r=1", source assumed to fill it: N scales as 1/filling); thin; compared with the thin-LTE HCO+ column for the HCO+ flux in the same window (a lower limit: HCO+ 3-2 is likely thick) {HW:.3f} Jy km/s (same assumptions)')
RQ={}
def Qrot_asym(T,A,B,C_,spin=1.0): return spin*np.sqrt(np.pi*(KB*T/H0/1e9)**3/(A*B*C_))/1.0
for nm,f,Eu,logA,smu,rot in (('HDO 2(2,0)-3(1,3)',266.16107,157.21292,-4.75649,None,(701.9,272.9,192.1)),('CH3CHO 12(6,6)-13(5,8) E',266.2132629,153.34802,-4.96295,None,(56.4,10.15,9.1))):
    # find catalogue sijmu2
    ll=[x for x in L if abs(float(x['freq_GHz'] or 0)-f)<2e-4 and x['name'].split()[0] in nm]
    smu=float(ll[0]['sijmu2']); A=10**logA; g=gu(f,smu,A); RQ[nm]=dict(g_u_from_Smu2=g)
    P(f' {nm}: nu={f} GHz Eup={Eu:.0f} K logA={logA}  g_u from catalogue S mu^2 = {g:.1f}')
    for T in (30,60,100,200):
        W=Wk(E,f,1.0); Nu=Nu_from_W(W,f,A)
        Q=Qrot_asym(T,*rot)   # classical rigid-rotor Q; catalogue S mu^2 gives g_u=5 (HDO, 2J+1, no nuclear spin) so no spin factor; CH3CHO order-of-magnitude only
        Ntot=Nu*Q/g*np.exp(Eu/T)
        # HCN thin column for the whole HCN window flux at same T (Eu=25.5 K, A=8.36e-4, g_u=21 (3 nuclear x 7), Q=T/2.127+1/3)
        Wh=Wk(HW,267.558,1.0); Nh=Nu_from_W(Wh,267.558,1.45e-3)*(T/2.14+1/3)/7*np.exp(25.7/T)
        RQ[nm][f'T{T}']=dict(W_Kkms=W,N_tot=Ntot,N_HCO_thin=Nh,ratio=Ntot/Nh)
        P(f'    T={T:3d} K: W_exc={W:.3f} K km/s (filling aperture)  N_tot(req) = {Ntot:.2e} cm-2 ; thin N(HCO+) for the HCO+ window flux = {Nh:.2e}; N_req/N(HCO+)_thin = {Ntot/Nh:.1f}')
RES['required_columns']=RQ
# ---------------------------------------------------------------- Part 4: N wing limits
P('\n# Part 4: dense-gas wing limits at N (aperture r=1.0", Xb4d at 1.05"), 3 sigma from random-aperture scatter of the same window; fractions of the line core flux (|v-v_loc|<225). Baselines: blue 225-350 from A (-900..-450); blue 350-750 and red 350-700 from (-1000..-780 U 500..700 for blue; -900..-450 for red).')
P4={}
for ous in ('Xb4d','Xb4b'):
    o=prep(ous,1.05); H,O=o['HCN'],o['HCO+']
    # alternative baseline D for deep blue window
    RP=randpos(o,60,1.0,5.0,9.5,excl=[H['pxy'](*pos(p)) for p in POS])
    xy=H['pxy'](*pos('N'))
    P4[ous]={}
    for mol,R in (('HCN',H),('HCO+',O)):
        sp=spec(R,apmask(R,xy,1.0)); v=R['v']
        core=win(R,sp,-225,225)[0]
        # baseline D shift: subtract extra per-pixel offset so that median over (-1000..-780)U(500..700) is zero -> apply to spectrum via difference between baselines
        cube=R['c']; refD=((v>-1000)&(v<-780))|((v>500)&(v<700)); off=lambda m_:np.median(cube[refD][:,m_],axis=0).mean() if False else None
        spD=sp-np.median(sp[refD])*np.ones_like(sp)    # aperture-summed spectrum with baseline D (median over D channels)
        outs={}
        for nm,(a,b),spx in (('blue225_350',(-350,-225),sp),('blue350_750',(-750,-350),spD),('red350_700',(350,700),sp)):
            val=win(R,spx,a,b)[0]
            nzv=[]
            for xyr in RP:
                s_=spec(R,apmask(R,xyr,1.0))
                s_=s_ if spx is sp else s_-np.median(s_[refD])
                nzv.append(win(R,s_,a,b)[0])
            nzv=np.array(nzv); sg=nzv.std()
            outs[nm]=dict(flux=val,sigma=sg,lim3=3*sg,lim3_frac_core=3*sg/core,flux_frac_core=val/core)
            P(f'  {ous} {mol:5s} {nm}: flux {val:+.3f}+-{sg:.3f} Jy km/s ({val/sg:+.1f} sigma), 3-sigma limit {3*sg:.3f} = {100*3*sg/core:.1f}% of core (|v|<225) flux {core:.2f}')
        P4[ous][mol]=dict(core=core,**outs)
    del o
RES['part4']=P4
json.dump(RES,open(OUT+'/s70.json','w'),indent=1,default=float)
open(OUT+'/s70.txt','w').write('\n'.join(txt))
