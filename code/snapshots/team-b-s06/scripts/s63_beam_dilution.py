# Q3 round 4, fix 1: beam-aware N vs S vs Bridge. Two compact cores (N, S, measured at 1.7" in Xb4d 3-2, team-a ours_1p7arcsec.csv) + a smooth component, convolved with the 10" Gaussian beam of the ladder data.
# (a) geometry: N-S separation, beam overlap weights at N, S, Bridge; (b) share of the 10" 3-2 flux that is in the two cores; (c) deblend: solve the three ladders (r_N, r_S, r_E; J/3-2 ratio in K km/s) from the three measured 10" ratios;
# (d) forward dilution: intrinsic N/S ratio-of-ratios -> 10" N/S ratio-of-ratios; (e) RADEX mapping of a Tk or n contrast between the cores; (f) the 6.5" (Xb47/Xb4f) T43/T32 check; (g) the detectability threshold.
# usage: python3 -I s63_beam_dilution.py -> out/q3d/s63_beam_dilution.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L, s62_lib as S
from regions import REGIONS
OUT='/workspace/out/q3d/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
res={}
cd=np.cos(np.radians(43.904))
def sep(a,b): return np.hypot((REGIONS[a][0]-REGIONS[b][0])*cd*3600,(REGIONS[a][1]-REGIONS[b][1])*3600)
FW=10.0; g=lambda d,fw=FW: np.exp(-4*np.log(2)*(d/fw)**2)
d_NS=sep('N','S'); d_NB=sep('N','Bridge'); d_SB=sep('S','Bridge')
P(f'# (a) geometry: N-S separation {d_NS:.2f}", N-Bridge {d_NB:.2f}", S-Bridge {d_SB:.2f}"  (10" FWHM Gaussian beam: overlap weight of a source at distance d, g=exp(-4ln2 (d/10)^2))')
P(f'   weight of the S core in the aperture at N: {g(d_NS):.3f}; of the N core at S: {g(d_NS):.3f}; Bridge sees N {g(d_NB):.3f}, S {g(d_SB):.3f}.  N-S is {d_NS/10:.2f} of the beam FWHM: the N, S and Bridge 10" apertures are not independent.')
P(f'   at the 6.5" beam (Xb47 4-3 vs Xb4f 3-2): weight {g(d_NS,6.5):.3f}; at 4.3" (Xb47 native major axis): {g(d_NS,4.3):.3f}')
# cores at 1.7" (K km/s; Xb4d 3-2 convolved to 1.7", moment 0): HCN, HCO+
CORE={'N':{'H':63.61,'O':85.06},'S':{'H':20.28,'O':53.72}}; OM17=1.1331*1.7**2; OM10=L.OM10
res['geometry']=dict(d_NS=d_NS,d_NB=d_NB,d_SB=d_SB,w_cross=g(d_NS),w_NB=g(d_NB),w_SB=g(d_SB))
def frac(pos,mol,jj,corescale=(1.,1.),wfun=g):
    """3-2 content (K km/s at 10") of: N core, S core, smooth component at position pos"""
    I,e=L.DATA[pos]; T=I[1] if mol=='H' else I[4]
    d={'N':{'N':0,'S':d_NS,'Bridge':d_NB},'S':{'N':d_NS,'S':0,'Bridge':d_SB}}
    cN=CORE['N'][mol]*OM17/OM10*corescale[0]*wfun(d['N'][pos]); cS=CORE['S'][mol]*OM17/OM10*corescale[1]*wfun(d['S'][pos])
    return cN,cS,max(T-cN-cS,0.),T
P('\n# (b) share of the measured 10" HCN/HCO+ 3-2 (Xb4f) in the two 1.7" cores (core flux = 1.7" beam intensity x 1.7" beam solid angle, i.e. a lower limit if the cores are larger than 1.7")')
FR={}
for mol,nm in (('H','HCN'),('O','HCO+')):
    for pos in ('N','Bridge','S'):
        cN,cS,e,T=frac(pos,mol,0); FR[(mol,pos)]=(cN/T,cS/T,e/T)
        P(f'   {nm:5s} 3-2 at {pos:6s}: N core {cN/T:.2f}  S core {cS/T:.2f}  smooth {e/T:.2f}  (T10={T:.2f} K km/s)')
def ratios(pos,mol): I=L.DATA[pos][0]; o=0 if mol=='H' else 3; return np.array([I[o]/I[o+1],I[o+2]/I[o+1]])
# (c) deblend ladders
P('\n# (c) deblend: observed 10" J/3-2 ratios at N, Bridge, S = sum_c f_c r_c with the fractions above (cores N,S + smooth). Ill-conditioned because the smooth fraction is nearly the same at all three positions; MC over core flux +-30% (each), beam FWHM 10+-0.3", position of cores +-0.3", stat errors; ratios a calibration uncertainty common to all three positions cancels.')
rng=np.random.default_rng(3); DB={}
for mol,nm in (('H','HCN'),('O','HCO+')):
    for jn,j in (('T21/T32',0),('T43/T32',1)):
        sol=[];
        for it in range(3000):
            cs=(1+0.3*rng.standard_normal(2)).clip(0.3); fw=10+0.3*rng.standard_normal(); dn=d_NS+0.3*rng.standard_normal()
            A=[];b=[]
            for pos in ('N','Bridge','S'):
                I,e=L.DATA[pos]; T=I[1] if mol=='H' else I[4]
                dd={'N':(0,dn),'Bridge':(d_NB,d_SB),'S':(dn,0)}[pos]
                cN=CORE['N'][mol]*OM17/OM10*cs[0]*g(dd[0],fw); cS=CORE['S'][mol]*OM17/OM10*cs[1]*g(dd[1],fw); ex=max(T-cN-cS,0.)
                A.append([cN/T,cS/T,ex/T]); o=0 if mol=='H' else 3
                rr=(I[o]/I[o+1],I[o+2]/I[o+1])[j]*(1+rng.standard_normal()*np.hypot(e[o+(0,2)[j]],e[o+1])); b.append(rr)
            try: sol.append(np.linalg.solve(np.array(A),np.array(b)))
            except Exception: pass
        sol=np.array(sol); q=np.percentile(sol,[16,50,84],axis=0)
        ro=[ratios(p,mol)[j] for p in ('N','Bridge','S')]
        DB[(mol,jn)]=sol
        P(f'   {nm:5s} {jn}: observed N/Bridge/S = {ro[0]:.3f}/{ro[1]:.3f}/{ro[2]:.3f};  core N {q[1,0]:.2f} [{q[0,0]:.2f},{q[2,0]:.2f}]  core S {q[1,1]:.2f} [{q[0,1]:.2f},{q[2,1]:.2f}]  smooth {q[1,2]:.2f} [{q[0,2]:.2f},{q[2,2]:.2f}]; core N/S {np.median(sol[:,0]/sol[:,1]):.2f} [{np.percentile(sol[:,0]/sol[:,1],16):.2f},{np.percentile(sol[:,0]/sol[:,1],84):.2f}]')
res['deblend']={f'{m}|{j}':dict(coreN=np.percentile(s_[:,0],[16,50,84]).tolist(),coreS=np.percentile(s_[:,1],[16,50,84]).tolist(),smooth=np.percentile(s_[:,2],[16,50,84]).tolist(),NoverS=np.percentile(s_[:,0]/s_[:,1],[16,50,84]).tolist()) for (m,j),s_ in DB.items()}
P('   -> the deblend of the cores from 10" data alone has errors of tens of per cent per core; it is a consistency screen, not a measurement of the nuclear ladders.')
# (d) forward dilution
P('\n# (d) forward dilution: intrinsic core ratio contrast rho = r_N/r_S for a ladder ratio, smooth component r_E equal to r_S ("S-like envelope") or to sqrt(r_N r_S); predicted 10" N/S ratio-of-ratios (position N / position S) and N/Bridge')
def dil(mol,rho,env='S'):
    out={}
    rN,rS=rho,1.; rE=rS if env=='S' else (np.sqrt(rN*rS) if env=='mid' else rN)
    r={}
    for pos in ('N','Bridge','S'):
        fN,fS,fE=FR[(mol,pos)]; r[pos]=fN*rN+fS*rS+fE*rE
    return r['N']/r['S'],r['N']/r['Bridge'],r['Bridge']/r['S']
for mol,nm in (('H','HCN'),('O','HCO+')):
    for env in ('S','mid'):
        P(f'   {nm} envelope ladder = {"S core" if env=="S" else "geometric mean"}:')
        for rho in (1.2,1.5,2.0,3.0,5.0):
            a,b,c=dil(mol,rho,env); P(f'      intrinsic N/S = {rho:.1f} -> 10" N/S = {a:.3f}   N/Bridge {b:.3f}   Bridge/S {c:.3f}')
        # contrast needed for +32 % / +8 %
        from scipy.optimize import brentq
        for thr in (1.08,1.32):
            try: rho=brentq(lambda r:dil(mol,r,env)[0]-thr,1.0,200.); P(f'      needed intrinsic N/S for a 10" N/S ratio-of-ratios of {thr:.2f}: {rho:.2f}')
            except Exception: P(f'      needed intrinsic N/S for {thr:.2f}: >200')
        res[f'dil|{nm}|{env}']={f'{r}':dil(mol,r,env)[0] for r in (1.2,1.5,2.0,3.0,5.0)}
# (e) RADEX mapping
F=S.Grid(); I0,e0=L.DATA['N']; pos=[dict(y=np.log(I0),rel=e0,cls=None)]
o=S.mcfit(F,pos,nd=64,keep_samples=False); Pp=S.Post(F,o['lnpost'],S.prior_w(F,'logT300'))
def run_ref(label,mask):
    global ref,rows
    rows=[]
    imax=int(np.argmax(np.where(mask,Pp.p,0))); ref=dict(Tk=10**F.lT[imax],ln=F.ln[imax],lNO=F.lNO[imax],lA=F.lA[imax],tauH=10**F.tauH[imax,1])
    P(f'\n# (e) [{label}] RADEX mapping from the single-position (N) posterior mode: Tk={ref["Tk"]:.0f} K, log n={ref["ln"]:.2f}, log N(HCO+)={ref["lNO"]:.1f}, log A={ref["lA"]:.1f} (tau_HCN32={ref["tauH"]:.1f}); ladder from the RADEX grid, filling cancels in J/3-2')
    def ladder_at(lT,ln):
        """T21/T32, T43/T32 for HCN and HCO+ at (logT, log n, ref columns), via the grid interpolation"""
        k=np.where((np.abs(F.lT-lT)<0.026)&(np.abs(F.ln-ln)<0.07)&(np.abs(F.lNO-ref['lNO'])<0.01)&(np.abs(F.lA-ref['lA'])<0.01))[0]
        if len(k)==0: return None
        m=F.M[k[0]]; return np.array([np.exp(m[0]-m[1]),np.exp(m[2]-m[1]),np.exp(m[3]-m[4]),np.exp(m[5]-m[4])])
    lab=['HCN T21/T32','HCN T43/T32','HCO+ T21/T32','HCO+ T43/T32']
    rows=[]
    def scen(lTS,lnS,lTN,lnN):
        rN=ladder_at(lTN,lnN); rS=ladder_at(lTS,lnS)
        if rN is None or rS is None: return None
        out=[]
        for k,(mol,j) in enumerate((('H',0),('H',1),('O',0),('O',1))):
            rr=rN[k]/rS[k]
            # 10" N/S with smooth component = S-like and the mean-weighted envelope
            fN,fS,fE=FR[(mol,'N')]; gN,gS,gE=FR[(mol,'S')]
            a=(fN*rN[k]+fS*rS[k]+fE*rS[k])/(gN*rN[k]+gS*rS[k]+gE*rS[k]); out.append((rr,a))
        return out
    lT0=np.log10(ref['Tk']); lT0=LT=S.LT[np.argmin(abs(S.LT-lT0))]
    P('   (intrinsic core N/S ratio-of-ratios -> diluted 10" N/S; smooth envelope with S-core ladder; threshold 2 sigma = 1.32)')
    for nm,lTN,lnN,lTS,lnS in (('Tk_N=1.5 Tk_S',lT0+np.log10(1.5),ref['ln'],lT0,ref['ln']),('Tk_N=2 Tk_S',lT0+np.log10(2),ref['ln'],lT0,ref['ln']),('Tk_N=3 Tk_S',lT0+np.log10(3),ref['ln'],lT0,ref['ln']),
                               ('Tk_S=Tk_N/1.5',lT0,ref['ln'],lT0-np.log10(1.5),ref['ln']),('Tk_S=Tk_N/2',lT0,ref['ln'],lT0-np.log10(2),ref['ln']),('Tk_S=Tk_N/3',lT0,ref['ln'],lT0-np.log10(3),ref['ln']),
                               ('n_N=3 n_S',lT0,ref['ln']+0.477,lT0,ref['ln']),('n_N=10 n_S',lT0,ref['ln']+1.0,lT0,ref['ln']),('n_S=n_N/3',lT0,ref['ln'],lT0,ref['ln']-0.477),('n_S=n_N/10',lT0,ref['ln'],lT0,ref['ln']-1.0)):
        lTN=S.LT[np.argmin(abs(S.LT-lTN))]; lTS=S.LT[np.argmin(abs(S.LT-lTS))]; lnN=S.LN[np.argmin(abs(S.LN-lnN))]; lnS=S.LN[np.argmin(abs(S.LN-lnS))]
        o_=scen(lTS,lnS,lTN,lnN)
        if o_ is None: P(f'   {nm}: grid hole'); continue
        P(f'   {nm:14s} (Tk_N,Tk_S={10**lTN:.0f},{10**lTS:.0f} K; log n_N,n_S={lnN:.2f},{lnS:.2f}): '+'; '.join(f'{lab[k]} intrinsic N/S {a:.2f} -> 10" {b:.2f}' for k,(a,b) in enumerate(o_)))
        rows.append(dict(name=nm,vals=[(float(a),float(b)) for a,b in o_]))

    res.setdefault('radex_map',{})[label]=rows; res.setdefault('ref',{})[label]={k:float(v) for k,v in ref.items()}
run_ref('mode of the single-position posterior',np.ones(F.nm,bool))
run_ref('mode with log n>=4.5 (dense, less thick branch)',F.ln>=4.5)
run_ref('mode with Tk<=60 K (cool branch)',10**F.lT<=60)
# measured 10" N/S
P('\n# measured 10" N/S ratio-of-ratios (position N / position S):')
for mol,nm in (('H','HCN'),('O','HCO+')):
    rn=ratios('N',mol); rs=ratios('S',mol); P(f'   {nm}: T21/T32 N/S = {rn[0]/rs[0]:.3f}; T43/T32 N/S = {rn[1]/rs[1]:.3f}  (stat <~1.5%)')
    res[f'meas|{nm}']=[float(rn[0]/rs[0]),float(rn[1]/rs[1])]
# calibration-limited 2 sigma: for ratio-of-ratios between positions in the same bands the calibration cancels; the 32% quoted by the reviewer is the detectability for any ladder ratio against absolute (independent bands: sqrt(2*0.1^2+...)) 
P('   detectability: 32% (2 sigma of an *absolute* ladder-ratio difference between two datasets with independent 10%+10% band calibrations and 7% per-line terms). For positions in the same cubes the shared calibration cancels, so the positional ratio-of-ratios is limited by stat (<2%) and by the 3-5% beam/pointing terms: the measured 0.94-1.06 is then a robust *10" statement*; it is nonetheless not informative about the cores because of the dilution above.')
# (f) 6.5" check: T43/T32 at 6.5" (Xb47 J4 6.5as / Xb4f J3 6.5as)
P('\n# (f) T43/T32 at the 6.5" common beam (Xb47 4-3 / Xb4f 3-2; weights 0.18 for the other nucleus) for N, Bridge, S')
def tk65(tag,row,m,k,ek):
    r=L.R30[tag]['rows'][row]; om=1.1331*6.5**2; return L.tk(r[k],m,{'Xb47_J4_6.5as':4,'Xb4f_J3_6.5as':3}[tag],om),L.tk(r[ek],m,{'Xb47_J4_6.5as':4,'Xb4f_J3_6.5as':3}[tag],om)
d65={}
for pos in ('N_pos','Bridge_pos','S_pos'):
    o_=[]
    for m,k,ek in (('HCN','HCN','eHCN'),('HCO+','HCOp','eHCOp')):
        a,ea=tk65('Xb47_J4_6.5as',pos,m,k,ek); b,eb=tk65('Xb4f_J3_6.5as',pos,m,k,ek); o_.append((a/b,a/b*np.hypot(ea/a,eb/b)))
    d65[pos]=o_
    P(f'   {pos[:-4]:6s}: HCN T43/T32 {o_[0][0]:.3f}+-{o_[0][1]:.3f}  HCO+ {o_[1][0]:.3f}+-{o_[1][1]:.3f}')
for mi,nm in ((0,'HCN'),(1,'HCO+')):
    P(f'   {nm}: N/S at 6.5" = {d65["N_pos"][mi][0]/d65["S_pos"][mi][0]:.3f} (stat {np.hypot(d65["N_pos"][mi][1]/d65["N_pos"][mi][0],d65["S_pos"][mi][1]/d65["S_pos"][mi][0]):.3f}); at 10": {res["meas|"+nm][1]:.3f}')
    res[f'ratio65|{nm}']=float(d65["N_pos"][mi][0]/d65["S_pos"][mi][0])
open(OUT+'s63_beam_dilution.txt','w').write('\n'.join(log)); json.dump(res,open(OUT+'s63_beam_dilution.json','w'),indent=1,default=float)
