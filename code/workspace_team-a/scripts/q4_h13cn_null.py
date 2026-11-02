# Q4b-4: calibration of the H13CN 4-3 'hint' at N: (i) the same pipeline (same wing models, windows, template = HCN 4-3 profile at the nucleus) run at line-free CONTROL positions
# (annulus around the field centre, >=4.5" from every catalogued source, 4.6x3.1" common beam) -> empirical distribution of the fitted template flux (expressed as N-equivalent alpha = flux/I_HCN(N));
# (ii) injection-recovery: HCN-shaped line with alpha_inj = 0.03, 0.075 (of the local nucleus HCN 4-3 flux) injected into N, S, Bridge and control spectra, refit with every wing model.
# Fit: joint nonlinear least squares of alpha + wing parameters (multi-start), error from the Jacobian; inadequate models (zero CO, Lorentzian) dropped as in q4_h13cn.py.
# Output: out/q4b/h13cn_null_controls.csv, h13cn_null_injection.csv, h13cn_null_summary.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q4_h13cn import *
from q4_lib import *
from scipy.optimize import least_squares
OUTB='/workspace/out/q4b/'
MODELS=['exp wing','log-quadratic wing','log-cubic wing','power-law @CO','CO=broadened HCO+43 template']
def fit_alpha(P,y,t,nm):
    fn,npar,starts=models(P)[nm]; sg=P['sg']; best=None
    for p0 in starts:
        try: rr=least_squares(lambda p:(y-p[0]*t-fn(p[1:]))/sg,[0.02]+list(p0),max_nfev=300,bounds=(([-1]+[-np.inf]*npar,[1]+[np.inf]*npar) if nm!='CO=broadened HCO+43 template' else ([-1,0,5],[1,1e3,300])))
        except Exception: continue
        if best is None or rr.cost<best.cost: best=rr
    if best is None: return np.nan,np.nan,np.nan
    J=best.jac; 
    try: cov=np.linalg.inv(J.T@J); e=np.sqrt(cov[0,0])
    except Exception: e=np.nan
    return best.x[0],e,2*best.cost
def posdata(refP,ra,de,fmin,fmax):
    """control spectrum at (ra,de); template and HCO+ profile are those of the reference region (N)."""
    f20,s20,_,_=spectrum('Xb47','spw20','N',TB,pos=(ra,de)); s20=s20*1e3
    P=dict(refP); m=(f20>fmin)&(f20<fmax); P['f']=f20[m]; P['y']=s20[m]; P['sg']=robust_sigma(s20[(f20>340.7)&(f20<341.7)])
    return P
if __name__=='__main__':
    rng=np.random.default_rng(3)
    refs={r:prep(r,341.80,342.28) for r in ('N','S','Bridge')}
    cen=(0.5*(REG.loc['N','RA_deg']+REG.loc['S','RA_deg']),0.5*(REG.loc['N','Dec_deg']+REG.loc['S','Dec_deg']))
    cosd=np.cos(np.radians(cen[1])); src=[(REG.loc[r,'RA_deg'],REG.loc[r,'Dec_deg']) for r in REG.index if r!='Disk']
    pts=[]
    for dx in np.arange(-15,15.1,3.5):
        for dy in np.arange(-15,15.1,3.5):
            r=np.hypot(dx,dy)
            if r<7 or r>15: continue
            ra=cen[0]+dx/3600/cosd; de=cen[1]+dy/3600
            if min(np.hypot((ra-a)*cosd*3600,(de-b)*3600) for a,b in src)<4.5: continue
            pts.append((ra,de,dx,dy))
    print('control positions',len(pts))
    rows=[]; IT=refs['N']['IT']
    for fx in (342.28,342.20):
        for (ra,de,dx,dy) in pts:
            P=posdata(refs['N'],ra,de,341.80,fx); t=refs['N']['T'](((1-P['f']/F13)*C)-P['v0'])
            # local CO-wing proxy: mean flux in 342.30-342.39 GHz (just below the band edge)
            for nm in MODELS:
                a,e,chi=fit_alpha(P,P['y'],t,nm); rows.append(dict(kind='control',dx=dx,dy=dy,fx=fx,model=nm,alpha=a,err=e,chi2=chi,ndata=len(P['f']),sg=P['sg']))
    C_=pd.DataFrame(rows); C_.round(4).to_csv(OUTB+'h13cn_null_controls.csv',index=False)
    # injection
    inj=[]; 
    sel=[pts[i] for i in rng.choice(len(pts),min(8,len(pts)),replace=False)]
    for fx in (342.28,342.20):
        for reg in ('N','S','Bridge','ctrl'):
            for pos in ([(None,None)] if reg!='ctrl' else [(p[0],p[1]) for p in sel]):
                if reg=='ctrl': P=posdata(refs['N'],pos[0],pos[1],341.80,fx); tref=refs['N']
                else: P=prep(reg,341.80,fx); tref=refs[reg]
                t=tref['T'](((1-P['f']/F13)*C)-P['v0']); P0=dict(P)
                for ai in (0.,0.03,0.075):
                    y=P0['y']+ai*t
                    for nm in MODELS:
                        a,e,chi=fit_alpha(P0,y,t,nm); inj.append(dict(region=reg,pos=str(pos),fx=fx,alpha_inj=ai,model=nm,alpha_rec=a,err=e,chi2=chi,ndata=len(P['f'])))
    I_=pd.DataFrame(inj); I_.round(4).to_csv(OUTB+'h13cn_null_injection.csv',index=False)
    print(C_.groupby(['fx','model']).alpha.describe().round(3).to_string())
