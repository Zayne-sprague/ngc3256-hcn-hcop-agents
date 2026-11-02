# Q16r2 (3): forward-model "ratio field" deblend (replaces the failed template deblend of s71/s77, whose chi2/dof is 20-300).
# J=3 is reproduced EXACTLY by construction: Xb4f J3 map (beam B0) = conv(Xb4d J3 map -> B0) + R, R = residual (extended gas missed by Xb4d, any sign kept).
# J=4 model = conv( Xb4d_J3(x) * rho(x) -> B0 ) * kappa + rho_h * R * kappa,  kappa=(nu4/nu3)^2 (flux units),  rho(x) = T43/T32 ratio field
#   rho(x) = rho_h + a_N g_N(x) + a_S g_S(x) [+ a_B g_B(x)],  g_c = exp(-d_c^2/2 s^2) (s = transition radius of the core excess), linear in (rho_h,a_N,a_S,a_B) for fixed s.
# Compared with the observed Xb47 J4 map convolved to B0. Core ratio of S = rho(S) = rho_h + a_S. Errors: map rms plus a 5% model floor; covariance propagated.
# Variants: s, Bridge term on/off, Xb4d template HCN or HCO+, fit radius, positivity of the model.   usage: python3 -I s78_q16_ratiofield.py -> out/q16b/s78.txt/json
import sys,json,itertools; sys.path.insert(0,'/workspace/scripts')
from q16b_lib import *
OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
def to_coarse(tm,tag_fine,mol_f,ec,nc):
    e,n=grid(tag_fine,mol_f); Y,X=np.mgrid[:tm.shape[0],:tm.shape[1]]; A_=np.c_[X.ravel(),Y.ravel(),np.ones(X.size)]
    ce=np.linalg.lstsq(A_,e.ravel(),rcond=None)[0]; cn=np.linalg.lstsq(A_,n.ravel(),rcond=None)[0]
    Mx=np.array([[ce[0],ce[1]],[cn[0],cn[1]]]); b=np.array([ce[2],cn[2]]); xy=np.linalg.solve(Mx,np.vstack([ec.ravel()-b[0],nc.ravel()-b[1]]))
    return map_coordinates(tm,[xy[1],xy[0]],order=1,mode='nearest').reshape(ec.shape)
def plane(ec,nc,d,r0,r1):
    m=(np.hypot(ec,nc)>=r0)&(np.hypot(ec,nc)<r1); A=np.c_[np.ones(m.sum()),ec[m],nc[m]]
    p=np.linalg.lstsq(A,d[m],rcond=None)[0]; return p[0]+p[1]*ec+p[2]*nc, float(np.std(d[m]-A@p))
def setup(mol,mol_t,thr=4.0,ring=(12,22),sigfloor=None):
    """thr: Xb4d template pixels kept only above thr*sigma (no clipping of noise: a clipped noisy template gives a large positive bias after convolution to 5.7x3.5");
       ring: local plane baseline fitted to both coarse maps in this radius range (None = no baseline); sigfloor: map error floor = empirical rms of the ring after plane removal"""
    ec,nc=grid('4f3',mol); e,n=grid('4d3',mol_t)
    m4=M[f'4d3|{mol_t}|m0']; X=np.clip(m4,0,None) if thr is None else np.where(m4>thr*M[f'4d3|{mol_t}|e'],m4,0.0); ps=float(M[f'4d3|{mol_t}|ps']); b0=tuple(M[f'4d3|{mol_t}|beam'])
    cX=to_coarse(conv_map(X,ps,b0,B0),'4d3',mol_t,ec,nc)
    d3=sample('4f3',mol,ec,nc); d4=sample('47_4',mol,ec,nc); e4=sample('47_4',mol,ec,nc,'e')
    sf=0.0
    if ring is not None:
        p3,s3=plane(ec,nc,d3,*ring); p4,s4=plane(ec,nc,d4,*ring); d3=d3-p3; d4=d4-p4; sf=max(s3*0.0,s4)
    if sigfloor is not None: sf=sigfloor
    e4=np.hypot(e4,sf)
    R=d3-cX
    return dict(ec=ec,nc=nc,e=e,n=n,X=X,ps=ps,b0=b0,cX=cX,d3=d3,d4=d4,e4=e4,R=R,mol=mol,mol_t=mol_t,sf=sf)
def design(D,s,bridge,Rfit):
    mol=D['mol']; kap=(NU[mol][4]/NU[mol][3])**2; cols=[]
    cols.append(kap*(D['cX']+D['R']))     # rho_h
    for c in (('N','S','Bridge') if bridge else ('N','S')):
        g=np.exp(-((D['e']-POS[c][0])**2+(D['n']-POS[c][1])**2)/(2*s*s))
        cols.append(kap*to_coarse(conv_map(D['X']*g,D['ps'],D['b0'],B0),'4d3',D['mol_t'],D['ec'],D['nc']))
    sel=np.hypot(D['ec'],D['nc'])<Rfit
    return np.array([c[sel] for c in cols]).T,sel
def solve(D,A,sel,floor):
    y=D['d4'][sel]; err=np.sqrt(D['e4'][sel]**2+(floor*np.abs(y))**2)
    Aw=A/err[:,None]; yw=y/err; cov=np.linalg.inv(Aw.T@Aw); p=cov@Aw.T@yw; chi=float(np.sum((Aw@p-yw)**2))
    return p,cov,chi,len(y)-len(p)
if __name__=='__main__':
    P('# T43/T32 (Tb) at the nuclei from a forward-modelled ratio field. rho(S)=rho_h+a_S, rho(N)=rho_h+a_N, rho_h = ratio of everything else (halo/Bridge/lobe).  matched-beam aperture (s72): HCN N .450 S .357 Br .362; HCO+ N .565 S .486 Br .494;  team-a core: HCN N .395 S .268; HCO+ N .482 S .442')
    res={}
    for mol in ('HCN','HCO+'):
        for mol_t in ('HCN','HCO+'):
            D=setup(mol,mol_t)
            for s,bridge,Rfit,floor in itertools.product((1.0,1.5,2.0,3.0),(False,True),(8,12),(0.05,)):
                A,sel=design(D,s,bridge,Rfit); p,cov,chi,dof=solve(D,A,sel,floor)
                rh=p[0]; rN=rh+p[1]; rS=rh+p[2]
                eN=np.sqrt(cov[0,0]+cov[1,1]+2*cov[0,1]); eS=np.sqrt(cov[0,0]+cov[2,2]+2*cov[0,2]); eh=np.sqrt(cov[0,0])
                txt=f'{mol:5s} tmpl {mol_t:5s} s={s} bridge={bridge!s:5} Rfit={Rfit}: rho_h {rh:.3f}+-{eh:.3f}  rho_N {rN:.3f}+-{eN:.3f}  rho_S {rS:.3f}+-{eS:.3f}'
                if bridge: rB=rh+p[3]; eB=np.sqrt(cov[0,0]+cov[3,3]+2*cov[0,3]); txt+=f'  rho_Br {rB:.3f}+-{eB:.3f}'
                txt+=f'  chi2 {chi:.0f} dof {dof} (pixels are ~%d-fold oversampled)'%46
                P(txt); res[f'{mol}|{mol_t}|{s}|{bridge}|{Rfit}']=dict(p=p.tolist(),cov=cov.tolist(),chi2=chi,dof=dof)
    json.dump(res,open(OUT+'s78.json','w'),indent=1); open(OUT+'s78.txt','w').write('\n'.join(log))
