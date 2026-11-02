# Q5e: kinematics at native 0.25": per-pixel Gaussian line fits (v, sigma), planar velocity gradients, peak-pixel and aperture line widths (HCN, HCO+ 3-2), beam-smearing correction.
# python3 -I scripts/q5e_kin.py -> out/q5/q5e.log, q5e_pix.npz, q5e_grad.csv, q5e_width.csv, q5e_kin.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=load(); dx=d['dx']; dy=d['dy']; pix=float(d['pix']); log=open(O+'q5e.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
SRC={'S':(-0.15,-2.60),'N':(0.0,2.60),'Bridge':(0.3,0.15)}
def gauss(v,A,v0,s): return A*np.exp(-0.5*((v-v0)/s)**2)
def fit_spec(v,s,e,v00=2780.,win=(2280,3280)):
    m=(v>win[0])&(v<win[1]); vv=v[m]; ss=s[m]
    i=np.argmax(ss); p0=[ss[i],vv[i],40.]
    try:
        p,c=curve_fit(gauss,vv,ss,p0=p0,sigma=np.full(m.sum(),e),absolute_sigma=True,bounds=([0,vv.min(),5],[10,vv.max(),300]),maxfev=4000)
    except Exception: return None
    return p,np.sqrt(np.diag(c)),(ss-gauss(vv,*p))/e
# ---- per-pixel fits on a 2-pixel grid (0.094") within r<3.3" of N, S, Bridge
PIX={}
for k,nm in (('H','HCN'),('O','HCO+')):
    a=d[k]; v=d['v'+k]; rms=d['rms'+k]; out={}
    for src,(x,y) in SRC.items():
        R=3.3 if src!='Bridge' else 1.6
        iy,ix=np.nonzero(np.hypot(dx-x,dy-y)<R); sel=(iy%2==0)&(ix%2==0); iy,ix=iy[sel],ix[sel]
        rows=[]
        for a_,b_ in zip(iy,ix):
            s=a[:,a_,b_]; e=rms[a_,b_]
            # 3-channel boxcar peak S/N test avoids fitting noise
            sm=np.convolve(s,np.ones(3)/3,'same')[(v>2430)&(v<3130)]
            if sm.max()<4*e/np.sqrt(3): continue
            f=fit_spec(v,s,e)
            if f is None: continue
            p,pe,res=f
            if p[0]/pe[0]<5 or pe[1]>25 or p[2]>250: continue
            rows.append((dx[a_,b_],dy[a_,b_],p[0],p[1],pe[1],p[2],pe[2],p[0]*p[2]*np.sqrt(2*np.pi)))
        out[src]=np.array(rows)
        pr(f'{nm} {src}: {len(rows)} fitted pixels (peak S/N>5, sigma_v err cut)')
    PIX[nm]=out
np.savez(O+'q5e_pix.npz',**{f'{nm}_{s}':PIX[nm][s] for nm in PIX for s in PIX[nm]})
# ---- planar gradient fits v = v0 + gx*dx + gy*dy (weights 1/err_v, floor 3 km/s) in r<R
rows=[]
Nb=1.13*0.245*0.235/pix**2/4      # beam area in sampled (2-px) pixels -> effective independent samples
def gradfit(P,x,y,R):
    m=np.hypot(P[:,0]-x,P[:,1]-y)<R
    if m.sum()<8: return None
    Q=P[m]; w=1/np.maximum(Q[:,4],3.)
    A=np.c_[np.ones(len(Q)),Q[:,0]-x,Q[:,1]-y]
    c,*_=np.linalg.lstsq(A*w[:,None],Q[:,3]*w,rcond=None)
    res=Q[:,3]-A@c; sres=np.std(res); neff=max(len(Q)/Nb,3.)
    cov_=np.linalg.inv((A*w[:,None]).T@(A*w[:,None]))*np.mean((res*w)**2)
    sc=np.sqrt(len(Q)/neff)    # inflate for beam-correlated pixels
    g=np.hypot(c[1],c[2]); pa=np.degrees(np.arctan2(c[1],c[2]))%360   # direction of increasing velocity, PA east of north
    eg=np.sqrt(cov_[1,1]*c[1]**2+cov_[2,2]*c[2]**2)/g*sc
    return dict(n=len(Q),v0=c[0],gx=c[1],gy=c[2],grad=g,e_grad=eg,pa_red=pa,pa_axis=pa%180,res_rms=sres,mean_sigma=np.average(Q[:,5],weights=Q[:,2]),n_eff=neff)
for nm in PIX:
    for src,(x,y) in SRC.items():
        for R in (0.4,0.6,0.8,1.0,1.4,2.0):
            r=gradfit(PIX[nm][src],x,y,R)
            if r: rows.append(dict(line=nm,src=src,R=R,**r))
gd=pd.DataFrame(rows); gd.to_csv(O+'q5e_grad.csv',index=False)
pd.set_option('display.width',220)
pr('\n# planar velocity gradients (km/s per arcsec); pa_red = PA (E of N) of the direction of increasing velocity; mean_sigma = intensity-weighted Gaussian sigma of the pixel fits (km/s)')
pr(gd.round(2).to_string())
# ---- line widths: peak pixel (one beam, native) and apertures r=0.12", 0.5", 1.0"; single Gaussian
rows=[]
from scipy import ndimage as ndi
for k,nm in (('H','HCN'),('O','HCO+')):
    a=d[k]; v=d['v'+k]; rms=d['rms'+k]; m0,_=win_map(d,k,*VW)
    for src,(x,y) in SRC.items():
        for rad in (0.0,0.12,0.25,0.5,1.0):
            if rad==0:      # brightest pixel of the 0.1"-smoothed m0 within 0.4" of the nominal position (the peak pixel, as used for 'one beam' values)
                msk=np.hypot(dx-x,dy-y)<0.4; z=np.where(msk,ndi.gaussian_filter(m0,2),-1e9); iy,ix=np.unravel_index(np.argmax(z),z.shape)
                s=a[:,iy,ix]; e=rms[iy,ix]; cx,cy=dx[iy,ix],dy[iy,ix]
            else:
                msk=np.hypot(dx-cx_ref(src) if False else dx-x,dy-y)<rad; s=a[:,msk].mean(1); e=np.median(rms[msk])/np.sqrt(max(msk.sum()/29.,1)); cx,cy=x,y
            f=fit_spec(v,s,e)
            if f is None: continue
            p,pe,res=f; fw=2.3548*p[2]
            # non-parametric: FWHM from the measured profile around peak (half maximum crossings, interpolated)
            m_=(v>2280)&(v<3280); vv=v[m_]; ss=ndi.uniform_filter1d(s[m_],3); ip=np.argmax(ss); h=ss[ip]/2
            lo=ip
            while lo>0 and ss[lo]>h: lo-=1
            hi=ip
            while hi<len(ss)-1 and ss[hi]>h: hi+=1
            rows.append(dict(line=nm,src=src,aperture_r=rad,x=cx,y=cy,peak_mJy=p[0]*1e3,e_peak=pe[0]*1e3,v0=p[1],e_v0=pe[1],sigma=p[2],e_sigma=pe[2],fwhm=fw,e_fwhm=2.3548*pe[2],fwhm_halfmax=abs(vv[hi]-vv[lo]),rchi2=np.mean(res**2)))
wd=pd.DataFrame(rows); wd.to_csv(O+'q5e_width.csv',index=False)
pr('\n# line widths (single Gaussian fit +-500 km/s about 2780; aperture_r=0 -> peak pixel, native beam)')
pr(wd.round(1).to_string())
# ---- figure: v and sigma maps from the pixel fits
fig,ax=plt.subplots(2,4,figsize=(16,9))
for i,nm in enumerate(('HCN','HCO+')):
    for j,(src,(x,y)) in enumerate((('N',SRC['N']),('S',SRC['S']))):
        P=PIX[nm][src]
        sc=ax[i,2*j].scatter(P[:,0],P[:,1],c=P[:,3],s=14,cmap='RdBu_r',vmin=2680,vmax=2880); ax[i,2*j].set_title(f'{nm} {src} v0'); plt.colorbar(sc,ax=ax[i,2*j])
        sc=ax[i,2*j+1].scatter(P[:,0],P[:,1],c=P[:,5],s=14,cmap='viridis',vmin=10,vmax=100); ax[i,2*j+1].set_title(f'{nm} {src} sigma'); plt.colorbar(sc,ax=ax[i,2*j+1])
        for q in (2*j,2*j+1): ax[i,q].set_xlim(x+2.5,x-2.5); ax[i,q].set_ylim(y-2.5,y+2.5); ax[i,q].set_aspect('equal')
plt.savefig(O+'q5e_kin.png',dpi=80)
