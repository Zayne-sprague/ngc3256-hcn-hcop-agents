# Q16r2 (4): injection-recovery of a team-a-like template deblend (NNLS, templates from the Xb4d 3-2 map, same template used at J=3 and J=4; components N, S, Bridge, Outer).
# Truth model (known): J3 map = the observed (plane-subtracted) Xb4f map; J4 = kappa*[ conv(X*rho) + rho_h*R ], X = Xb4d template (>4 sigma), R = Xb4f - conv(X) (extended gas),
# rho(x) = rho_h + a_N g_N + a_S g_S (g = gaussian of width s_true around each nucleus). Correlated noise (beam-smoothed, rms = 0.3 Jy/beam km/s) added to J4 (and J3).
# Output: recovered T43/T32 of N, S (median and 16-84% over noise realisations) vs input rho(N), rho(S) and vs the 'aperture' ratio of the truth map at the N, S positions.
# usage: python3 -I s79_q16_injection.py -> out/q16b/s79.txt/json
import sys,json; sys.path.insert(0,'/workspace/scripts')
from s78_q16_ratiofield import *
from scipy.ndimage import gaussian_filter
OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
def deblend(D,d3,d4,rc=1.5,rcB=1.0,Rfit=10,noise=0.3,floor=0.05,outer=True,thr=4.0,mol_t='HCN'):
    mol=D['mol']; ec,nc,e,n=D['ec'],D['nc'],D['e'],D['n']; X=D['X']
    dd=lambda p: np.hypot(e-POS[p][0],n-POS[p][1])
    mN=dd('N')<rc; mS=dd('S')<rc; mB=(dd('Bridge')<rcB)&~mN&~mS; comps={'N':mN,'S':mS,'Bridge':mB}
    if outer: comps['Outer']=~(mN|mS|mB)
    cols=[to_coarse(conv_map(np.where(mk,X,0),D['ps'],D['b0'],B0),'4d3',D['mol_t'],ec,nc) for mk in comps.values()]
    sel=np.hypot(ec,nc)<Rfit; A=np.array([c[sel] for c in cols]).T; out={}
    for nm,d in (('3',d3),('4',d4)):
        y=d[sel]; err=np.sqrt(noise**2+(floor*np.abs(y))**2); a,_=nnls(A/err[:,None],y/err); out[nm]=a
    return {k:(ratio43(out['4'][i],out['3'][i],mol) if out['3'][i]>1e-6 else np.nan) for i,k in enumerate(comps)}
def truth(D,rho_h,rho_N,rho_S,rho_B,s):
    mol=D['mol']; kap=(NU[mol][4]/NU[mol][3])**2
    g=lambda c: np.exp(-((D['e']-POS[c][0])**2+(D['n']-POS[c][1])**2)/(2*s*s))
    rho=rho_h+(rho_N-rho_h)*g('N')+(rho_S-rho_h)*g('S')+(rho_B-rho_h)*g('Bridge')
    j4=kap*(to_coarse(conv_map(D['X']*rho,D['ps'],D['b0'],B0),'4d3',D['mol_t'],D['ec'],D['nc'])+rho_h*D['R'])
    return D['d3'].copy(),j4
def aper(d3,d4,D,p,mol):
    x=lambda d: float(sample_arr(d,D,p)); return ratio43(x(d4),x(d3),mol)
def sample_arr(d,D,p):
    from scipy.ndimage import map_coordinates
    ec,nc=D['ec'],D['nc']; Y,X=np.mgrid[:d.shape[0],:d.shape[1]]; A=np.c_[X.ravel(),Y.ravel(),np.ones(X.size)]
    ce=np.linalg.lstsq(A,ec.ravel(),rcond=None)[0]; cn=np.linalg.lstsq(A,nc.ravel(),rcond=None)[0]
    Mx=np.array([[ce[0],ce[1]],[cn[0],cn[1]]]); b=np.array([ce[2],cn[2]]); xy=np.linalg.solve(Mx,np.array([POS[p][0]-b[0],POS[p][1]-b[1]]))
    return map_coordinates(d,[[xy[1]],[xy[0]]],order=1)[0]
def noise_map(rng,shape,rms=0.3):
    z=gaussian_filter(rng.standard_normal(shape),2.0); return z/z.std()*rms
if __name__=='__main__':
    mol='HCN'; D=setup(mol,mol,ring=None)       # no plane removal: truth is built on the plane-free maps, so noise-free recovery isolates the method
    D=setup(mol,mol,ring=(12,22))
    P('# observed-data deblend with this code (plane-subtracted maps, noise 0.6 = empirical ring rms): '+str({k:round(v,3) for k,v in deblend(D,D['d3'],D['d4'],noise=0.6).items()}))
    P('# team-a p0012: N 0.395 S 0.268;  matched-beam aperture: N 0.450 S 0.357 Bridge 0.362')
    # scenarios: (rho_h, rho_N, rho_S, rho_B, s)
    SC={'A uniform halo/S/Bridge 0.36, N 0.45':(0.36,0.45,0.36,0.36,1.5),'B same, cores sharper (s=1)':(0.36,0.45,0.36,0.36,1.0),'C S core 0.268 (team-a), halo 0.36':(0.36,0.45,0.268,0.36,1.5),
        'D S core 0.45, halo 0.36':(0.36,0.45,0.45,0.36,1.5),'E halo 0.20 (extended gas poor in J4), S 0.36':(0.20,0.45,0.36,0.30,1.5),'F all 0.40':(0.40,0.40,0.40,0.40,1.5),
        'G halo 0.36, N 0.55, S 0.40, Br 0.36, s=2.5':(0.36,0.55,0.40,0.36,2.5),'H S 0.36 halo 0.30 N 0.5 s=1.5':(0.30,0.50,0.36,0.33,1.5)}
    res={}
    for nm,(rh,rN,rS,rB,s) in SC.items():
        d3,d4=truth(D,rh,rN,rS,rB,s)
        ap={p:aper(d3,d4,D,p,mol) for p in ('N','S','Bridge')}
        rng=np.random.default_rng(3); rec={k:[] for k in ('N','S','Bridge')}; rec0=deblend(D,d3,d4,noise=0.3)
        for it in range(20):
            r=deblend(D,d3+noise_map(rng,d3.shape),d4+noise_map(rng,d4.shape),noise=0.3)
            for k in rec: rec[k].append(r.get(k,np.nan))
        q={k:np.nanpercentile(v,[16,50,84]) for k,v in rec.items()}
        P(f'{nm:48s} input rho_N {rN:.3f} rho_S {rS:.3f} | truth aperture (5.7x3.5") N {ap["N"]:.3f} S {ap["S"]:.3f} Br {ap["Bridge"]:.3f} | deblend noiseless N {rec0["N"]:.3f} S {rec0["S"]:.3f} | recovered N {q["N"][1]:.3f} [{q["N"][0]:.3f},{q["N"][2]:.3f}]  S {q["S"][1]:.3f} [{q["S"][0]:.3f},{q["S"][2]:.3f}]  -> S bias {q["S"][1]-rS:+.3f}')
        res[nm]=dict(inp=[rh,rN,rS,rB,s],aperture=ap,noiseless={k:float(v) for k,v in rec0.items()},rec={k:q[k].tolist() for k in q})
    json.dump(res,open(OUT+'s79.json','w'),indent=1); open(OUT+'s79.txt','w').write('\n'.join(log))
