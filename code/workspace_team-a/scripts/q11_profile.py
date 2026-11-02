# Q11 forward-model profile likelihood on the intrinsic disk gradient eps of HCN/HCO+, per band/array map.
# For each map (2-1 Xb55; 3-2 Xb4f; 3-2 Xb4d): free sky morphology (h, core FWHM, core fraction), transfer-function scale f (u_a=f*Bmin/lambda, Bmin=8.7 m, ramp 2x),
# free R0; eps FIXED on a grid; chi2 = sum over 5 annuli x 2 lines (0-4,4-8,8-12,12-16,16-24"), sigma = 5%|F|+1 Jy km/s (conservative: zero-level / edge systematics).
# True-ratio convention: R(r)=R0*(1-eps*clip((r-4)/6,0,1)); the annulus ratio of the TRUE (unfiltered, 10" beam) sky is also reported.
import sys; sys.path.insert(0,'/workspace/scripts')
from q11_forward import *
from scipy.optimize import differential_evolution
from multiprocessing import Pool
MAPS={'21 Xb55':('HCN 2-1','HCO+ 2-1',177.3,55.),'32 Xb4f':('HCN 3-2@Xb4f','HCO+ 3-2@Xb4f',265.9,37.),'32 Xb4d':('HCN 3-2@Xb4d','HCO+ 3-2@Xb4d',265.9,22.)}
EPS=(-0.1,0.,0.1,0.2,0.3,0.4)
OBS={k:(obs_annuli(v[0])[0],obs_annuli(v[1])[0]) for k,v in MAPS.items()}
def mod(p,eps,fr,pbw,Bm=8.7):
    h,sc,fc,f=p; T=Tf(fr,Bm,f); sky=sky_hco(h,sc,fc)
    return annuli(observe(sky*Rr(Rc,1.,eps),T,pbw)),annuli(observe(sky,T,pbw))
def run(args):
    mp,eps=args; kh,ko,fr,pbw=MAPS[mp]; yH,yO=OBS[mp]; sH=0.05*abs(yH)+1; sO=0.05*abs(yO)+1
    def chi(p):
        a,b=mod(p,eps,fr,pbw); A=(a*yH/sH**2).sum()/(a*a/sH**2).sum()   # amplitude shared: HCN = A*R0*..., HCO+ = B*..
        B=(b*yO/sO**2).sum()/(b*b/sO**2).sum()
        return (((A*a-yH)/sH)**2).sum()+(((B*b-yO)/sO)**2).sum()
    r=differential_evolution(chi,[(3,40),(1.5,10),(0.02,0.95),(0.05,3)],seed=2,maxiter=50,popsize=10,tol=1e-7)
    a,b=mod(r.x,eps,fr,pbw); A=(a*yH/sH**2).sum()/(a*a/sH**2).sum(); B=(b*yO/sO**2).sum()/(b*b/sO**2).sum()
    # true (unfiltered) 10"-beam annulus ratio of the best-fit sky
    sky=sky_hco(*r.x[:3]); T1=np.ones_like(U); th=annuli(observe(sky*Rr(Rc,1.,eps),T1,1e4)); to=annuli(observe(sky,T1,1e4))
    return mp,eps,r.fun,r.x.tolist(),(A*a/(B*b)).tolist(),(th/to*A/B).tolist()
if __name__=='__main__':
    jobs=[(m,e) for m in MAPS for e in EPS]
    with Pool(6) as P: res=P.map(run,jobs)
    json.dump(dict(obs={k:[v[0].tolist(),v[1].tolist()] for k,v in OBS.items()},res=res),open('/workspace/out/q11/q11_profile.json','w'),indent=1)
    for mp,eps,chi,x,rm,rt in res:
        print(f'{mp} eps={eps:+.1f} chi2={chi:7.1f} h={x[0]:.1f} sc={x[1]:.1f} fc={x[2]:.2f} f={x[3]:.2f} | filtered-model R {np.round(rm[:3],3)} | true R {np.round(rt[:3],3)}  obs R {np.round(OBS[mp][0][:3]/OBS[mp][1][:3],3)}')
