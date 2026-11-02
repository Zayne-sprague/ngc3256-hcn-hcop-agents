# Q5 library: Xb4b native cubes, empirical noise fields, beam-convolved elliptical-Gaussian model fits with noise-injection errors.
import numpy as np
from scipy.optimize import least_squares
D40=40.; PC=194.0   # pc per arcsec (D=40 Mpc, Pereira-Santaella+2024)
O='/workspace/out/q5/'
VW=(2430.,3130.)    # km/s moment window (2780+-350)
def load(ous='Xb4b'):
    return np.load(O+f'cube_native_{ous}.npz')
def win_map(d,k,lo,hi):
    v=d['v'+k]; m=(v>=lo)&(v<hi); return d[k][m].sum(0)*float(d['dv'+k]),m.sum()
def null_maps(d,k,nch_target,nblk=None):
    """independent line-free sums (same #channels as the measurement window) -> list of noise maps; uses FREE ranges of q3_lib (line-free after baseline)"""
    free={'H':[(1380,2100)],'O':[(1280,1680),(1880,2100)]}[k]; v=d['v'+k]; dv=float(d['dv'+k]); out=[]
    for a,b in free:
        idx=np.nonzero((v>=a)&(v<=b))[0]; n=len(idx)
        nb=max(1,n//nch_target) if n>=nch_target else 0
        if nb==0: continue
        step=n//nb
        for i in range(nb): out.append(d[k][idx[i*step:i*step+nch_target]].sum(0)*dv)
    return out
def null_maps_scaled(d,k,nch_target,nsub=70):
    """if the free ranges are shorter than the window, build blocks of nsub channels and rescale by sqrt(n_target/nsub)"""
    free={'H':[(1380,2100)],'O':[(1280,1680),(1880,2100)]}[k]; v=d['v'+k]; dv=float(d['dv'+k]); out=[]
    for a,b in free:
        idx=np.nonzero((v>=a)&(v<=b))[0]
        for i in range(len(idx)//nsub): out.append(d[k][idx[i*nsub:(i+1)*nsub]].sum(0)*dv*np.sqrt(nch_target/nsub))
    return out
def cov(smaj,smin,pa):
    t=np.deg2rad(pa); R=np.array([[np.sin(t),np.cos(t)],[np.cos(t),-np.sin(t)]])   # x=East(+RA offset), y=North; PA from north through east: major axis direction (sin PA, cos PA)
    return R@np.diag([smaj**2,smin**2])@R.T
def model(par,X,Y,Bcov,ncomp,bg):
    """sum of components (F, x0, y0, smaj, smin, pa[deg]) [intrinsic sigma, arcsec], beam covariance Bcov; output Jy/beam-units (F in units*beam)"""
    z=np.zeros_like(X)
    for i in range(ncomp):
        F,x0,y0,sa,sb,pa=par[6*i:6*i+6]; C=cov(sa,sb,pa)+Bcov; Ci=np.linalg.inv(C)
        dx=X-x0; dy=Y-y0; q=Ci[0,0]*dx*dx+2*Ci[0,1]*dx*dy+Ci[1,1]*dy*dy
        z+=F*np.sqrt(np.linalg.det(Bcov)/np.linalg.det(C))*np.exp(-0.5*q)
    if bg: z=z+par[-3]+par[-2]*X+par[-1]*Y
    return z
def beamcov(bmaj,bmin,bpa): return cov(bmaj/2.3548,bmin/2.3548,bpa)
def fit(img,X,Y,mask,p0,lo,hi,Bcov,ncomp,bg=True,sig=None):
    w=1./sig if sig is not None else 1.
    f=lambda p:((model(p,X,Y,Bcov,ncomp,bg)-img)*w)[mask]
    return least_squares(f,p0,bounds=(lo,hi),x_scale='jac',max_nfev=2000)
def decon(par):
    """intrinsic FWHM (maj,min) in arcsec from sigma params, with the major axis first"""
    sa,sb=abs(par[3]),abs(par[4]); return 2.3548*max(sa,sb),2.3548*min(sa,sb)
