# Q2 round 2, fix 3: ONE kinematic decomposition for HCN and HCO+ 3-2 in the S lobe aperture (Xb4d 1.05", r=1.0" and 1.5", per-pixel baseline median(-900..-450)):
#  narrow disk-like Gaussian + broad Gaussian, velocity and width of each component SHARED between the lines, amplitudes free (model C0);
#  separate kinematics per line (model Cs, test of the shared assumption);  C0 + HDO Gaussian at the HDO 2(2,0)-3(1,3) 266.1611 GHz position (HCN-frame offset -336-vloc, FWHM 100-200 fixed) on HCN only (model C0+HDO).
# Component fluxes (Gaussian integrals), window fluxes of the broad component (-405..-150 and alternatives), extrapolated |v|<150 part, flux-weighted velocities. Errors: 200 noise-injection refits.
# Output out/q2r2/s77c.txt/json
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.optimize import least_squares
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); rng=np.random.default_rng(5)
g=lambda v,a,c,s: a*np.exp(-0.5*((v-c)/s)**2)
def spectra(rad,prod='std_Xb4d'):
    o={}
    for mol in ('HCN','HCO+'):
        R=load_any(prod,mol,1.05,12.0); v=R['v']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
        sp=aspec(c,R['ba'],apmask(R['n'],R['ps'],R['pxy'](*LOBE),rad)); lf=(v>-900)&(v<-450)
        o[mol]=(v,sp,sp[lf].std(),R['dv'])
    return o
SIG=lambda f:f/2.355
def model(p,v,kind,mol,vloc_hdo=None):
    # p layout depends on kind
    if kind=='C0':   # shared: vn,sn,vb,sb, An_H,Ab_H,An_O,Ab_O
        vn,sn,vb,sb,anH,abH,anO,abO=p[:8]; a,b=(anH,abH) if mol=='HCN' else (anO,abO); m=g(v,a,vn,sn)+g(v,b,vb,sb)
        if kind=='C0' and len(p)>8 and mol=='HCN': m=m+g(v,p[8],p[9],p[10])
        return m
    if kind=='Cs':   # per line: An,vn,sn,Ab,vb,sb (6 each)
        q=p[:6] if mol=='HCN' else p[6:12]; return g(v,q[0],q[1],q[2])+g(v,q[3],q[4],q[5])
def fit(S,kind,hdo=None,seeds=6):
    vs={m:S[m][0] for m in S}; use={m:(S[m][0]>-800)&(S[m][0]<650) for m in S}
    def res(p):
        r=[]
        for m in ('HCN','HCO+'):
            v=vs[m][use[m]]; y=S[m][1][use[m]]; e=S[m][2]
            if kind=='C0' and hdo is not None and m=='HCN': mm=g(v,p[4],p[0],p[1])+g(v,p[5],p[2],p[3])+g(v,p[8],hdo[0],hdo[1])
            elif kind=='C0': mm=g(v,p[4 if m=='HCN' else 6],p[0],p[1])+g(v,p[5 if m=='HCN' else 7],p[2],p[3])
            else: q=p[:6] if m=='HCN' else p[6:12]; mm=g(v,q[0],q[1],q[2])+g(v,q[3],q[4],q[5])
            r.append((y-mm)/e)
        return np.concatenate(r)
    best=None; a0=S['HCO+'][1].max()
    for k in range(seeds):
        if kind=='C0':
            x0=[35,60,rng.uniform(-250,-50),rng.uniform(80,180),a0*0.5,a0*0.15,a0,a0*0.15]; lo=[0,15,-450,50,0,0,0,0]; hi=[100,120,50,300,np.inf,np.inf,np.inf,np.inf]
            if hdo is not None: x0+= [a0*0.1]; lo+=[0]; hi+=[np.inf]
        else:
            x0=[a0*0.5,35,60,a0*0.15,rng.uniform(-250,-50),rng.uniform(80,180)]*2; lo=[0,0,15,0,-450,50]*2; hi=[np.inf,100,120,np.inf,50,300]*2
        try: r=least_squares(res,x0,bounds=(lo,hi))
        except Exception: continue
        if best is None or r.cost<best.cost: best=r
    return best,res
def comp(p,kind,hdo):
    out={}
    for mol in ('HCN','HCO+'):
        if kind=='C0': an,ab,vn,sn,vb,sb=(p[4],p[5],p[0],p[1],p[2],p[3]) if mol=='HCN' else (p[6],p[7],p[0],p[1],p[2],p[3])
        else: q=p[:6] if mol=='HCN' else p[6:12]; an,vn,sn,ab,vb,sb=q
        Fb=np.sqrt(2*np.pi)*ab*sb; Fn=np.sqrt(2*np.pi)*an*sn
        from scipy.special import erf
        def win(a,b): return Fb*0.5*(erf((b-vb)/(np.sqrt(2)*sb))-erf((a-vb)/(np.sqrt(2)*sb)))
        out[mol]=dict(Fbroad=Fb,Fnarrow=Fn,vb=vb,sb=sb,vn=vn,sn=sn,w_405_150=win(-405,-150),w_350_100=win(-350,-100),w_405_205=win(-405,-205),w_inner=win(-150,150),w_blue_all=win(-1000,0),
                      Fhdo=(np.sqrt(2*np.pi)*p[8]*hdo[1] if (hdo is not None and mol=='HCN' and kind=='C0') else 0.0))
    return out
RES={}
for rad in (1.0,1.5):
    S=spectra(rad); dv=S['HCN'][3]
    P(f'### lobe aperture r={rad}": rms/chan HCN {S["HCN"][2]*1e3:.2f}  HCO+ {S["HCO+"][2]*1e3:.2f} mJy (beam-averaged Jy per channel; 4.4 km/s channels)')
    fits={}
    for name,kind,hdo in (('C0','C0',None),('Cs','Cs',None),('C0+HDO','C0',(-365.0,100/2.355)),('C0+HDO195','C0',(-365.0,195/2.355))):
        r,resf=fit(S,kind,hdo); chi=2*r.cost; ndat=len(resf(r.x)); npar=len(r.x)
        # noise-injection errors: refit simulated spectra (best model + gaussian noise, noise correlated? white assumed; chi2/dof inflation used)
        dof=ndat-npar; infl=max(1.0,chi/dof)
        P(f'  {name}: chi2={chi:.0f}/{dof} (red {chi/dof:.2f}) params {np.round(r.x,3).tolist()}')
        fits[name]=(r,chi,dof)
        c=comp(r.x,kind,hdo)
        # MC
        draws=[]
        v1,v2=S['HCN'][0],S['HCO+'][0]
        for it in range(60):
            S2={}
            for m in S:
                v=S[m][0]; sp=S[m][1]; e=S[m][2]
                # model spectrum on full axis from fit
                if kind=='C0':
                    p=r.x; a,b=(p[4],p[5]) if m=='HCN' else (p[6],p[7]); mm=g(v,a,p[0],p[1])+g(v,b,p[2],p[3])
                    if hdo is not None and m=='HCN': mm=mm+g(v,p[8],hdo[0],hdo[1])
                else: q=r.x[:6] if m=='HCN' else r.x[6:12]; mm=g(v,q[0],q[1],q[2])+g(v,q[3],q[4],q[5])
                S2[m]=(v,mm+rng.normal(0,e*np.sqrt(infl),len(v)),e,S[m][3])
            try:
                r2,_=fit(S2,kind,hdo,seeds=2); draws.append(comp(r2.x,kind,hdo))
            except Exception: pass
        for mol in ('HCN','HCO+'):
            d={k:(c[mol][k],float(np.std([x[mol][k] for x in draws]))) for k in c[mol]}
            RES[f'{rad}|{name}|{mol}']=d
            P(f'    {mol:5s}: broad flux {d["Fbroad"][0]:.3f}+-{d["Fbroad"][1]:.3f} Jy km/s; centroid {d["vb"][0]:.0f}+-{d["vb"][1]:.0f}, sigma {d["sb"][0]:.0f}+-{d["sb"][1]:.0f} km/s | narrow F {d["Fnarrow"][0]:.2f} v={d["vn"][0]:.0f} sig={d["sn"][0]:.0f} | broad in -405..-150: {d["w_405_150"][0]:.3f}+-{d["w_405_150"][1]:.3f}; -350..-100: {d["w_350_100"][0]:.3f}; -405..-205: {d["w_405_205"][0]:.3f}; |v|<150: {d["w_inner"][0]:.3f}; v<0 total {d["w_blue_all"][0]:.3f}'+(f'; HDO line flux {d["Fhdo"][0]:.3f}+-{d["Fhdo"][1]:.3f}' if mol=='HCN' and 'HDO' in name else ''))
    dchi=fits['C0'][1]-fits['Cs'][1]; P(f'  shared-kinematics test: chi2(C0)-chi2(Cs) = {dchi:.1f} for {fits["Cs"][0].x.size-fits["C0"][0].x.size} fewer parameters (C0 has {len(fits["C0"][0].x)} vs Cs {len(fits["Cs"][0].x)})')
    P(f'  HDO test: chi2(C0)-chi2(C0+HDO100) = {fits["C0"][1]-fits["C0+HDO"][1]:.1f}; HDO195: {fits["C0"][1]-fits["C0+HDO195"][1]:.1f} (1 extra amplitude parameter)')
    RES[f'{rad}|chi']={k:(v[1],v[2]) for k,v in fits.items()}
json.dump(RES,open(OUT+'/s77c.json','w'),indent=1,default=float); open(OUT+'/s77c.txt','w').write('\n'.join(txt))
