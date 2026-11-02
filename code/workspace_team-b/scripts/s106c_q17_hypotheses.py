# Q17 round 2, step c: likelihood-type comparison of three origins for the red wing (+180..320) over the whole nuclear/Bridge region (reviewer fix 3).
# Data: s106a_maps_*.npz (r=1" aperture-mean maps on a 1" grid within 10" of the S-N midpoint; Xb4d, 1.05", baseline A). Gaussian noise sigma = sd of control apertures (field+matched-core scatter, below).
# Null M0: red = a*core + b + point-like wing terms at N and S nuclei (free amps >=0)  [generic core-scaled leakage + nuclear wings].
# H_A (S counter-lobe): + one component on the S outflow axis: PA of the lobe (183 +- 10 deg from S, antiparallel to the blue lobe), radius r free 0.5..6", intrinsic sigma 0..1".
# H_B (N red-only flow): + one component within 2.5" of N, position free (ANY direction); H_B' restricts to the counter-direction of M18 North-OF-Blue (PA 111+-20 from N).
# H_C (Bridge streaming): + a uniform-amplitude strip along the straight N-S line (width sigma 0..1.5", free amplitude; extends 0..1 of the way N<->S).
# Component template = Gaussian of sigma^2 = sigma_beam^2 + sigma_ap^2 + sigma_int^2 (beam 1.05" FWHM, r=1" aperture ~ sigma 0.5").  Compare chi2, BIC; also give per-hypothesis predictions (ridge flux, position).
# usage: python3 -I s106c_q17_hypotheses.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.optimize import least_squares
OUT='/workspace/out/q17r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
cosd=np.cos(np.deg2rad(S0[1])); rng=np.random.default_rng(3)
sb=1.05/2.355; sap=0.5
def rel(pos): return np.array([(pos[0]-S0[0])*3600*cosd,(pos[1]-S0[1])*3600])   # (E,N) arcsec from S
NR=rel(N0); print('N rel S',NR)
for mol in ('HCN','HCO+'):
    d=np.load(f'{OUT}/s106a_maps_{mol}.npz'); ps=float(d['ps']); sel=d['sel']&d['grid']; core=d['core']; red=d['red']
    R=load_any('std_Xb4d',mol,1.05,24.0); n=R['n']; Y,X=np.mgrid[:n,:n]; x0,y0=R['pxy'](*S0)
    E=-(X-x0)*ps; Nn=(Y-y0)*ps     # arcsec E, N rel S
    # restrict to sel & grid and core>0.1 (inside emission) ... keep all sel points (core-free pixels constrain b)
    m=sel&(core>-1); Ed=E[m]; Nd=Nn[m]; cd=core[m]; yd=red[m]; sig=0.027 if mol=='HCN' else 0.026
    def tmpl(pE,pN,si): s2=sb**2+sap**2+si**2; return np.exp(-0.5*((Ed-pE)**2+(Nd-pN)**2)/s2)
    def strip(si,f0,f1):
        s2=sb**2+sap**2+si**2; t=np.zeros_like(Ed)
        for f in np.linspace(f0,f1,15): t+=np.exp(-0.5*((Ed-f*NR[0])**2+(Nd-f*NR[1])**2)/s2)
        return t/15
    base=lambda p: p[0]*cd+p[1]+p[2]*tmpl(NR[0],NR[1],0)+p[3]*tmpl(0,0,0)
    def fit(model,p0,lo,hi):
        best=None
        for trial in range(12):
            q0=np.array(p0,float)*(1 if trial==0 else rng.uniform(0.5,1.5,len(p0)))
            q0=np.clip(q0,np.array(lo)+1e-6,np.array(hi)-1e-6)
            r=least_squares(lambda p:(model(p)-yd)/sig,q0,bounds=(lo,hi))
            if best is None or r.cost<best.cost: best=r
        return best
    ndat=len(yd); out={}
    r0=fit(base,[0.01,0,0.1,0.03],[-1,-0.1,0,0],[1,0.1,2,2]); chi0=2*r0.cost
    P(f'### {mol}: {ndat} grid points, sigma={sig}; M0 (core*a+b + N,S nuclear wings): chi2={chi0:.1f}, a={r0.x[0]:.4f} b={r0.x[1]:+.3f} N wing amp {r0.x[2]:.3f} S {r0.x[3]:.3f}')
    out['M0']=dict(chi2=chi0,k=4)
    # H_A: S axis
    def mA(p):
        a,b,aN,aS,A,r,si,pa=p; pE=r*np.sin(np.deg2rad(pa)); pN=r*np.cos(np.deg2rad(pa)); return a*cd+b+aN*tmpl(NR[0],NR[1],0)+aS*tmpl(0,0,0)+A*tmpl(pE,pN,si)
    rA=fit(mA,[0.01,0,0.1,0.03,0.2,3.3,0.3,0],[-1,-0.1,0,0,0,0.5,0,-10],[1,0.1,2,2,2,6,1,10]); 
    # H_B free-position near N
    def mB(p):
        a,b,aN,aS,A,dE,dN,si=p; return a*cd+b+aN*tmpl(NR[0],NR[1],0)+aS*tmpl(0,0,0)+A*tmpl(NR[0]+dE,NR[1]+dN,si)
    rB=fit(mB,[0.01,0,0.1,0.03,0.2,0,-1.8,0.3],[-1,-0.1,0,0,0,-2.5,-2.5,0],[1,0.1,2,2,2,2.5,2.5,1])
    # H_B' counter-direction of N blue lobe (PA 111 +-20 from N)
    def mBp(p):
        a,b,aN,aS,A,r,pa,si=p; return a*cd+b+aN*tmpl(NR[0],NR[1],0)+aS*tmpl(0,0,0)+A*tmpl(NR[0]+r*np.sin(np.deg2rad(pa)),NR[1]+r*np.cos(np.deg2rad(pa)),si)
    rBp=fit(mBp,[0.01,0,0.1,0.03,0.1,2,111,0.3],[-1,-0.1,0,0,0,0.3,91,0],[1,0.1,2,2,2,5,131,1])
    # H_C strip
    def mC(p):
        a,b,aN,aS,A,si=p; return a*cd+b+aN*tmpl(NR[0],NR[1],0)+aS*tmpl(0,0,0)+A*strip(si,0.0,1.0)
    rC=fit(mC,[0.01,0,0.1,0.03,0.1,0.5],[-1,-0.1,0,0,0,0],[1,0.1,2,2,2,1.5])
    # H_C2 strip only the part between the nuclei, excluding the nuclei themselves (0.15..0.85)
    def mC2(p):
        a,b,aN,aS,A,si=p; return a*cd+b+aN*tmpl(NR[0],NR[1],0)+aS*tmpl(0,0,0)+A*strip(si,0.15,0.85)
    rC2=fit(mC2,[0.01,0,0.1,0.03,0.1,0.5],[-1,-0.1,0,0,0,0],[1,0.1,2,2,2,1.5])
    for nm,r_,k,desc in (('H_A S counter-lobe (on axis PA 173-193)',rA,8,lambda x:f'A={x[4]:.3f} at r={x[5]:.2f}" PA={x[7]:+.0f} int sigma {x[6]:.2f}"'),
                         ('H_B N red-only flow (free position within 2.5" of N)',rB,8,lambda x:f'A={x[4]:.3f} at dE={x[5]:+.2f} dN={x[6]:+.2f} from N (=> {NR[0]+x[5]:+.2f},{NR[1]+x[6]:+.2f} from S) int sigma {x[7]:.2f}"'),
                         ('H_B2 N counter-lobe direction (PA 111+-20)',rBp,8,lambda x:f'A={x[4]:.3f} at r={x[5]:.2f}" PA {x[6]:.0f} from N'),
                         ('H_C Bridge strip N-S (0..1)',rC,6,lambda x:f'A={x[4]:.3f} width sigma {x[5]:.2f}"'),
                         ('H_C2 Bridge strip between nuclei (0.15-0.85)',rC2,6,lambda x:f'A={x[4]:.3f} width sigma {x[5]:.2f}"')):
        chi=2*r_.cost; dbic=(chi-chi0)+(k-4)*np.log(ndat)
        P(f'  {nm}: chi2={chi:.1f} (dchi2={chi-chi0:+.1f}, extra params {k-4}, dBIC={dbic:+.1f}); {desc(r_.x)}; core slope a={r_.x[0]:.4f}')
        out[nm]=dict(chi2=chi,k=k,dchi2=chi-chi0,dBIC=dbic,x=[float(t) for t in r_.x])
    # profile likelihood: chi2 vs position of a single component along the S-N axis (dE=0.15), amplitude free (>=0), width 0.3"
    P('  profile of chi2 vs position of a point-like comp on the S-N axis (dN rel S): '+', '.join(f'{pn:+.1f}:{2*fit(lambda p,pn=pn:p[0]*cd+p[1]+p[2]*tmpl(NR[0],NR[1],0)+p[3]*tmpl(0,0,0)+p[4]*tmpl(0.15,pn,0.3),[0.01,0,0.1,0.03,0.1],[-1,-0.1,0,0,0],[1,0.1,2,2,2]).cost-chi0:+.0f}' for pn in (-1,0.5,1.5,2.5,3.3,4.2,5.1,6.0,7.0)))
    res[mol]=out
    # model-independent: red wing at positions north of N (far side) and mirrors
    for nm,(pE,pN) in {'N far side (+1.85" N of N)':(NR[0],NR[1]+1.85),'N far side (+3" N of N)':(NR[0],NR[1]+3.0),'N counter-lobe dir (PA111, 2.4")':(NR[0]+2.4*np.sin(np.deg2rad(111)),NR[1]+2.4*np.cos(np.deg2rad(111))),'ridge':(0.15,3.3),'Bridge midpoint (+2.55")':(0.1,2.55),'Bridge +4.2"':(0.13,4.2)}.items():
        dd=np.hypot(E-pE,Nn-pN); ii=np.unravel_index(np.argmin(dd),dd.shape); P(f'     red wing {nm}: {red[ii]:.3f}, core {core[ii]:.2f}, blue wing {d["blue"][ii]:.3f}')
json.dump(res,open(OUT+'/s106c.json','w'),indent=1); open(OUT+'/s106c.txt','w').write('\n'.join(txt))
