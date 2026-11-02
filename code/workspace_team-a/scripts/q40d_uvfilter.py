# Q40d: matched-uv-range decomposition from the image products. Every map (1-0 Harada, 1-0 Michiyama, 2-1 Xb55, 3-2 Xb4f) is regridded to a common 0.4" grid, FT'd, multiplied by a common uv high-pass
# T(u)=0.5(1+erf((u-umin)/1.5kl)), converted to the common 8.64x6.0 beam (Fourier-space convolution kernel), and fitted with the SAME filtered beam-shaped templates (N,S,E1,E2,W; NNLS) -> linear filtering is
# forward-modelled exactly, so the missing-short-spacing bowl is absorbed in the templates. Noise: filtered null-window maps added and refit. umin=0 is the control (free-sign constant included).
import sys; sys.path.insert(0,'/workspace/scripts')
from q40c_uvcmp import *
from scipy.special import erf
from scipy.optimize import nnls
import q40b_fit as B
TB=B.TBS['b86']; KFJ=B.KF
SETS=[(1,'Harada'),(1,'Michiyama'),(2,'Xb55'),(3,'Xb4f')]
def OMf(b): return np.pi*b[0]*b[1]/(4*np.log(2))
gx=(np.arange(NG)-NG/2)*PIXG; XG,YG=np.meshgrid(gx,gx)
u=np.fft.fftshift(np.fft.fftfreq(NG,PIXG)); UU,VV=np.meshgrid(u,u); RK=np.hypot(UU,VV)*KL
def quad(S): return S[0,0]*UU**2+2*S[0,1]*UU*VV+S[1,1]*VV**2
def filt(d,arr,umin,tb):
    g=regrid(d,arr); F=np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(g)))
    S=B.cov2(*tb)-cov2(*d['beam'])
    if np.linalg.eigvalsh(S).min()<-1e-3: raise ValueError('beam')
    K=np.exp(-2*np.pi**2*quad(S))*OMf(tb)/OMf(d['beam'])
    T=0.5*(1+erf((RK-umin)/1.5)) if umin>0 else 1.0
    return np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F*K*T)))),T
def templ(names,umin,shift={},size=0.):
    out=[]
    for n in names:
        ex,ey=B.EXY[n]; sx,sy=shift.get(n,(0,0)); Si=np.linalg.inv(B.cov2(*TB,s=size))
        X=np.stack([XG-ex-sx,YG-ey-sy],-1); g=np.exp(-0.5*np.einsum('...i,ij,...j->...',X,Si,X))
        F=np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(g))); T=0.5*(1+erf((RK-umin)/1.5)) if umin>0 else 1.0
        out.append(np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F*T)))))
    return out
RMAX=14.
sel=np.hypot(XG,YG)<RMAX
def fitmap(y,cols,umin):
    A=np.array([c[sel] for c in cols]).T
    if umin==0: A=np.hstack([A,np.ones((A.shape[0],1)),-np.ones((A.shape[0],1))])
    x,_=nnls(A,y[sel]); return x[:len(cols)]
if __name__=='__main__':
    NAMES=('N','S','E1','E2','W'); out={}; log=[]
    def P(s): print(s,flush=True); log.append(s)
    DS={}
    for J,nm in SETS:
        for ln,lab in ((0,'HCN'),(1,'HCO+')): DS[(nm,lab)]=B.getmaps(nm,ln)
    for umin in (0,3,5,6,8,10):
        cols=templ(NAMES,umin); res={}
        for J,nm in SETS:
            for lab in ('HCN','HCO+'):
                d,m,ns=DS[(nm,lab)]; y,_=filt(d,m,umin,TB); a=fitmap(y,cols,umin)
                dd=[]
                for n in ns:
                    n=np.nan_to_num(n); yn,_=filt(d,n,umin,TB); dd.append(fitmap(y+yn,cols,umin)-a)
                res[(nm,lab)]=(a,np.std(dd,axis=0) if len(dd)>=3 else np.full(5,np.nan))
        out[umin]=res
        def R(nm,J,reg):
            def v(lab):
                a,e=res[(nm,lab)]
                if reg=='N+S': return a[0]+a[1],np.hypot(e[0],e[1])
                i=NAMES.index(reg); return a[i],e[i]
            h,q=v('HCN'),v('HCO+')
            if h[0]<=0 or q[0]<=0: return np.nan,np.nan
            r=h[0]/q[0]*KFJ[J]; return r,r*np.hypot(h[1]/h[0],q[1]/q[0])
        P(f'## umin={umin} kl, common beam 8.64x6.0 (R_K = HCN/HCO+)')
        for reg in ('N','S','N+S','E1','E2','W'):
            P(f' {reg:4s} '+'  '.join(f'{nm}: {R(nm,J,reg)[0]:.3f}+-{R(nm,J,reg)[1]:.3f}' for J,nm in SETS))
        for c in ('E1','W'):
            P(f' {c}/(N+S): '+'  '.join(f'{nm}: {R(nm,J,c)[0]/R(nm,J,"N+S")[0]:.2f}' for J,nm in SETS))
        P(' HCN E1 / HCO+ E1 flux (Jy km/s @b86): '+'  '.join(f'{nm}:{res[(nm,"HCN")][0][2]:.2f}/{res[(nm,"HCO+")][0][2]:.2f}' for J,nm in SETS))
    open('/workspace/out/q40b/q40d_uvfilter.log','w').write('\n'.join(log))
    import pickle; pickle.dump(out,open('/workspace/out/q40b/q40d_uvfilter.pkl','wb'))
