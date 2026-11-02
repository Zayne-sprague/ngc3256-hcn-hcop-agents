# Q16r2 (2): re-implement a team-a-like template deblend of the Xb4f 3-2 and (common-beam) Xb47 4-3 maps and test what it returns for S and N.
# Components (masks on the Xb4d 3-2 moment-0 map at 0.94x0.74"): N, S, Bridge (circles of radius rc/rcB), Outer = everything else inside 26". Each component map is
# convolved to the Xb4f beam (5.72x3.54 PA82) and sampled on the Xb4f grid. The Xb4f J3 map and the Xb47 J4 map (convolved to the same beam) are fit with NNLS
# amplitudes a_c(J) (same template for both J). T43/T32 of component c = a_c(4)/a_c(3) * (nu3/nu4)^2 (RJ, flux->Tb).
# Variants: template (Xb4d HCN or HCO+), radius rc, fit radius, outer component on/off, template clipping. -> out/q16b/s77.txt/json
import sys,json,itertools; sys.path.insert(0,'/workspace/scripts')
from q16b_lib import *
OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
def templates(mol_t,rc,rcB,outer=True,clip=True):
    m=M[f'4d3|{mol_t}|m0'].copy(); ps=float(M[f'4d3|{mol_t}|ps']); b0=tuple(M[f'4d3|{mol_t}|beam'])
    if clip: m=np.clip(m,0,None)
    e,n=grid('4d3',mol_t); C={}
    d=lambda p: np.hypot(e-POS[p][0],n-POS[p][1])
    mN=d('N')<rc; mS=d('S')<rc; mB=(d('Bridge')<rcB)&~mN&~mS
    comps={'N':mN,'S':mS,'Bridge':mB}
    comps['Outer']=~(mN|mS|mB)&(np.hypot(e,n)<26)
    out={}
    for k,mk in comps.items():
        if k=='Outer' and not outer: continue
        out[k]=conv_map(np.where(mk,m,0),ps,b0,B0)
    return out,e,n
def fit(mol,mol_t,rc=1.5,rcB=1.0,Rfit=14,outer=True,clip=True,floor=0.0,signs=None):
    T,e,n=templates(mol_t,rc,rcB,outer,clip)
    ec,nc=grid('4f3',mol)
    sel=np.hypot(ec,nc)<Rfit
    A=np.array([sample('4d3',mol_t,ec,nc,'m0')*0 for _ in T]) # placeholder shape
    # sample each convolved component onto the Xb4f grid (templates live on Xb4d grid)
    cols=[]
    ps=float(M[f'4d3|{mol_t}|ps'])
    for k,tm in T.items():
        # bilinear sample of tm on 4d grid
        Y,X=np.mgrid[:tm.shape[0],:tm.shape[1]]; A_=np.c_[X.ravel(),Y.ravel(),np.ones(X.size)]
        ce=np.linalg.lstsq(A_,e.ravel(),rcond=None)[0]; cn=np.linalg.lstsq(A_,n.ravel(),rcond=None)[0]
        Mx=np.array([[ce[0],ce[1]],[cn[0],cn[1]]]); b=np.array([ce[2],cn[2]]); xy=np.linalg.solve(Mx,np.vstack([ec.ravel()-b[0],nc.ravel()-b[1]]))
        cols.append(map_coordinates(tm,[xy[1],xy[0]],order=1,mode='nearest'))
    A=np.array(cols).T[sel.ravel()]
    res={}
    for tag in ('4f3','47_4'):
        # J4 map resampled onto the Xb4f grid
        d=sample(tag,mol,ec,nc,'m0').ravel()[sel.ravel()]
        err=np.sqrt(sample(tag,mol,ec,nc,'e').ravel()[sel.ravel()]**2+(floor*np.abs(d))**2)
        a,rn=nnls(A/err[:,None],d/err); res[tag]=(a,float(np.sum(((A@a-d)/err)**2)),len(d)-len(a))
    names=list(T)
    r={k:ratio43(res['47_4'][0][i],res['4f3'][0][i],mol) if res['4f3'][0][i]>0 else np.nan for i,k in enumerate(names)}
    return r,res,names
if __name__=='__main__':
    P('# team-a p0012 deblend core ladders T43/T32: HCN N 0.395 S 0.268 ; HCO+ N 0.482 S 0.442.   matched-beam aperture (s72): HCN N 0.450 S 0.357 Bridge 0.362; HCO+ N 0.565 S 0.486 Bridge 0.494')
    allr={}
    for mol in ('HCN','HCO+'):
        for mol_t in ('HCN','HCO+'):
            for rc,rcB,outer,Rfit in itertools.product((1.0,1.5,2.0,2.5),(1.0,),(True,False),(10,14,20)):
                r,res,names=fit(mol,mol_t,rc,rcB,Rfit,outer)
                allr[f'{mol}|{mol_t}|{rc}|{rcB}|{outer}|{Rfit}']=r
                P(f'{mol:5s} tmpl {mol_t:5s} rc={rc} outer={outer!s:5} Rfit={Rfit:2d}: '+' '.join(f'{k} {v:.3f}' for k,v in r.items())+f' | chi2/dof J3 {res["4f3"][1]:.0f}/{res["4f3"][2]} J4 {res["47_4"][1]:.0f}/{res["47_4"][2]}')
    json.dump(allr,open(OUT+'s77.json','w'),indent=1); open(OUT+'s77.txt','w').write('\n'.join(log))
