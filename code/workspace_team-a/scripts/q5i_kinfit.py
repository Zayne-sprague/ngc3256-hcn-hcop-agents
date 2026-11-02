# Q5i: forward-model kinematics. Intrinsic elliptical-Gaussian emission + rotation curve v=vsys+Vmax*(r.u)/sqrt(r^2+Rt^2) + sigma_int(R)=so+(sc-so)exp(-R/0.25"),
# model cube -> clean-beam convolution -> the SAME per-pixel Gaussian fit as the data (q5e_pix.npz) -> chi2 on (v0, sigma, peak) of the fitted pixels.
# python3 -I scripts/q5i_kinfit.py [N|S] -> out/q5h/q5i_<src>_<line>.json/.log, q5i_pv_<src>.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
import json
from scipy.signal import fftconvolve
from scipy.optimize import least_squares
log=None; MODE='gauss'
def pr(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True)
    if log: log.write(s+'\n'); log.flush()
d=load(); DX=d['dx']; DY=d['dy']; pix=float(d['pix'])
PIXD=np.load(O+'q5e_pix.npz')
SRC={'N':(-0.05,2.63),'S':(-0.15,-2.60)}
RFIT=1.0; HB=2.2       # fitted pixels r<RFIT; model box half-size HB (arcsec)
def setup(src,k):
    x,y=SRC[src]; nm='HCN' if k=='H' else 'HCO+'
    ix0=np.argmin(np.abs(DX[0]-x))   # x index
    iy0=np.argmin(np.abs(DY[:,0]-y)); h=int(HB/pix)
    sl=(slice(iy0-h,iy0+h+1),slice(ix0-h,ix0+h+1)); X=DX[sl]-x*0; Y=DY[sl]
    v=d['v'+k]; vm=(v>2230)&(v<3330); vv=v[vm]
    P=PIXD[f'{nm}_{src}']; P=P[np.hypot(P[:,0]-x,P[:,1]-y)<RFIT]
    # indices of the data pixels in the box
    ii=np.array([[np.argmin(np.abs(DY[sl][:,0]-p[1])),np.argmin(np.abs(DX[sl][0]-p[0]))] for p in P])
    bc=beamcov(*d['beam'+k]); hk=15; yy,xx=np.mgrid[-hk:hk+1,-hk:hk+1]*pix; Ci=np.linalg.inv(bc)
    ker=np.exp(-0.5*(Ci[0,0]*xx**2+2*Ci[0,1]*xx*yy+Ci[1,1]*yy**2)); ker/=ker.sum()
    m0,_=win_map(d,k,*VW); Im=np.clip(m0[sl],0,None); Im=Im/Im.max()
    return dict(Imap=Im,src=src,k=k,nm=nm,X=X,Y=Y,vv=vv,P=P,ii=ii,ker=ker,x=x,y=y,sl=sl,vm=vm)
def cubemodel(p,S):
    x0,y0,fa,fb,pai,pak,Vm,Rt,sc,so,vs,amp=p
    X=S['X']-x0; Y=S['Y']-y0
    t=np.deg2rad(pai); u=np.array([np.sin(t),np.cos(t)]); a=(X*u[0]+Y*u[1]); b=(-X*u[1]+Y*u[0])
    I=np.exp(-0.5*((a/(fa/2.3548))**2+(b/(fb/2.3548))**2)) if MODE=='gauss' else S['Imap']   # 'map': measured (already beam-smoothed) m0 map as the intrinsic intensity -> upper bound on the smear
    tk=np.deg2rad(pak); r=np.hypot(X,Y); s_=(X*np.sin(tk)+Y*np.cos(tk))
    vl=vs+Vm*s_/np.sqrt(r**2+Rt**2); sg=so+(sc-so)*np.exp(-r/0.25)
    vv=S['vv'][:,None,None]
    cube=I[None]*np.exp(-0.5*((vv-vl[None])/sg[None])**2)/(sg[None]*np.sqrt(2*np.pi))
    cube=fftconvolve(cube,S['ker'][None],mode='same',axes=(1,2))
    return cube*amp
def gfit(spec,vv,niter=25):
    """vectorised Gauss-Newton/LM fit of A exp(-(v-v0)^2/2s^2) to spec[npix,nch]; init from the peak channel and the intensity-weighted moments"""
    n=spec.shape[0]; ip=spec.argmax(1); A=spec.max(1); v0=vv[ip]; w=np.clip(spec,0,None); s=np.sqrt(np.clip((w*(vv[None]-v0[:,None])**2).sum(1)/np.clip(w.sum(1),1e-12,None),100,None)); s=np.clip(s,10,150)
    lam=np.full(n,1e-2); par=np.c_[A,v0,s]
    def res(par): return par[:,0:1]*np.exp(-0.5*((vv[None]-par[:,1:2])/par[:,2:3])**2)-spec
    r=res(par); c=(r**2).sum(1)
    for _ in range(niter):
        g=np.exp(-0.5*((vv[None]-par[:,1:2])/par[:,2:3])**2); dv_=(vv[None]-par[:,1:2])
        J=np.stack([g,par[:,0:1]*g*dv_/par[:,2:3]**2,par[:,0:1]*g*dv_**2/par[:,2:3]**3],-1)   # n,nch,3
        JTJ=np.einsum('nci,ncj->nij',J,J); JTr=np.einsum('nci,nc->ni',J,r)
        for _l in range(4):
            M=JTJ+lam[:,None,None]*np.einsum('nii->ni',JTJ)[:,:,None]*np.eye(3)[None]+1e-18*np.eye(3)[None]
            step=-np.linalg.solve(M,JTr[:,:,None])[:,:,0]; pn=par+step; pn[:,2]=np.clip(pn[:,2],3,300); rn=res(pn); cn=(rn**2).sum(1)
            good=cn<c; par=np.where(good[:,None],pn,par); r=np.where(good[:,None],rn,r); c=np.where(good,cn,c); lam=np.where(good,lam*0.3,lam*5)
    return par
def obsmodel(p,S):
    cube=cubemodel(p,S); sp=cube[:,S['ii'][:,0],S['ii'][:,1]].T       # npix,nch
    # restrict to the window used for the data fits (2280-3280)
    m=(S['vv']>2280)&(S['vv']<3280); return gfit(sp[:,m],S['vv'][m])
def resid(p,S):
    g=obsmodel(p,S); P=S['P']; ev=np.maximum(P[:,4],5.); es=np.maximum(P[:,6],5.); eA=np.maximum(0.1*P[:,2],3e-4)
    rv=(g[:,1]-P[:,3])/ev; rs=(g[:,2]-P[:,5])/es; ra=(g[:,0]-P[:,2])/eA
    return np.r_[rv,rs,ra]
def run(src,k,fixed_int=None):
    S=setup(src,k); x,y=S['x'],S['y']
    if src=='N': starts=[(85,0.4),(70,0.2),(100,0.7)]; vm0=150.
    else: starts=[(120,0.4),(100,0.2),(140,0.7)]; vm0=100.
    amp0=np.median(S['P'][:,2])*np.sqrt(2*np.pi)*40.
    best=None
    for pak,Rt in starts:
        p0=[x,y,0.45,0.25,110,pak,vm0,Rt,50,35,2780,amp0*3]
        lo=[x-.3,y-.3,0.05,0.03,0,pak-60,5,0.03,10,10,2600,amp0*0.05]; hi=[x+.3,y+.3,1.5,1.5,180,pak+60,600,3,150,100,2950,amp0*60]
        if MODE=='map': p0[2:5]=[0.45,0.25,110]; lo[2:5]=[0.449,0.249,109.9]; hi[2:5]=[0.451,0.251,110.1]
        try: r=least_squares(resid,p0,args=(S,),bounds=(lo,hi),x_scale=[.05,.05,.1,.1,20,20,30,.2,15,15,10,amp0],max_nfev=200,diff_step=1e-3)
        except Exception as e: pr('fail',e); continue
        if best is None or r.cost<best.cost: best=r
    r=best; n=len(r.fun); dof=n-12; chi2r=2*r.cost/dof
    J=r.jac; cov_=np.linalg.pinv(J.T@J)*chi2r
    Nb=1.13*d['beam'+k][0]*d['beam'+k][1]/(2*pix)**2; infl=np.sqrt(Nb)   # beam-correlated pixels: inflate errors
    err=np.sqrt(np.diag(cov_))*infl
    names=['x0','y0','fwhm_maj','fwhm_min','pa_int','pa_kin','Vmax_los','Rt','sig_c','sig_out','vsys','amp']
    pr(f'{src} {S["nm"]}: n_resid={n} chi2_red={chi2r:.2f} (v,sig,peak residual rms: {np.std(r.fun[:n//3]):.2f},{np.std(r.fun[n//3:2*n//3]):.2f},{np.std(r.fun[2*n//3:]):.2f} sigma-units); errors inflated by sqrt(Nb)={infl:.2f}')
    for nm_,v_,e_ in zip(names,r.x,err): pr(f'   {nm_:9s} {v_:9.3f} +- {e_:.3f}')
    # intrinsic sigma(R), v(R) with MC from the covariance
    rng=np.random.default_rng(1); draws=rng.multivariate_normal(r.x,cov_*Nb,400)
    Rg=np.array([0.0,0.12,0.25,0.5,0.75,1.0]); sig=[];vr=[]
    for q in draws:
        so,sc=q[9],q[8]; Rt=q[7]; sig.append(so+(sc-so)*np.exp(-Rg/0.25)); vr.append(abs(q[6])*Rg/np.sqrt(Rg**2+Rt**2))
    sig=np.array(sig); vr=np.array(vr)
    out=dict(src=src,line=S['nm'],chi2_red=chi2r,n=n,params=dict(zip(names,r.x)),err=dict(zip(names,err)),
             R=Rg.tolist(),sig_int=np.percentile(sig,[16,50,84],0).tolist(),vlos=np.percentile(vr,[16,50,84],0).tolist())
    for i_,R_ in enumerate(Rg): pr(f'   R={R_:.2f}": sigma_int={out["sig_int"][1][i_]:.0f} [{out["sig_int"][0][i_]:.0f},{out["sig_int"][2][i_]:.0f}]  v_los(R)={out["vlos"][1][i_]:.0f} [{out["vlos"][0][i_]:.0f},{out["vlos"][2][i_]:.0f}] km/s')
    # observed vs model per-pixel sigma near the peak, and a smearing diagnostic: model sigma with sigma_int forced to 0
    p0_=r.x.copy(); p0_[8]=p0_[9]=8.0   # small sigma (subtracted in quadrature below) -> pure rotation smear
    g0=obsmodel(p0_,S); g0[:,2]=np.sqrt(np.clip(g0[:,2]**2-64.,0,None)); gm=obsmodel(r.x,S); rr=np.hypot(S['P'][:,0]-S['x'],S['P'][:,1]-S['y'])
    for lo_,hi_ in ((0,0.2),(0.2,0.4),(0.4,0.7),(0.7,1.0)):
        m=(rr>=lo_)&(rr<hi_)
        if m.sum(): pr(f'   r in [{lo_},{hi_}): n={m.sum()} sigma_obs(data,mean)={S["P"][m,5].mean():.1f} model={gm[m,2].mean():.1f} pure-rotation-smear(sigma_int~0)={g0[m,2].mean():.1f}')
        out[f'smear_{lo_}_{hi_}']=float(g0[m,2].mean()) if m.sum() else None
    out['sig_obs_by_r']={f'{lo_}-{hi_}':float(S['P'][(rr>=lo_)&(rr<hi_),5].mean()) if ((rr>=lo_)&(rr<hi_)).sum() else None for lo_,hi_ in ((0,0.2),(0.2,0.4),(0.4,0.7),(0.7,1.0))}
    json.dump(out,open(O.replace('q5','q5h')+f'q5i_{MODE}_{src}_{S["nm"]}.json','w'),indent=1,default=float)
    return out
if __name__=='__main__':
    args=sys.argv[1:]
    if 'map' in args: MODE='map'; args.remove('map')
    log=open('/workspace/out/q5h/q5i_'+MODE+'_'+'_'.join(args or ['all'])+'.log','w')
    for src in (args or ['N','S']):
        for k in ('O','H'): run(src,k)
