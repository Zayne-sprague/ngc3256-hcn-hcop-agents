# Q3 round 2: two-component (narrow disk-like + broad blue) Gaussian decomposition of the HCN and HCO+ 3-2 spectra at the fixed Soutflow apertures (shared v, sigma between lines); ratio of each component. Output: out/q3b/twocomp.log/csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
import pandas as pd
from scipy.optimize import least_squares
cT=(F0['HCO+']/F0['HCN'])**2
log=open(OUT+'twocomp.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
rng=np.random.default_rng(1)
def data(x,y,rad,ous='Xb4d',fw=1.05,nb=3):
    out={}
    for sp in ('HCN','HCO+'):
        c=vcube(sp,ous,fw,crop=9. if ous=='Xb4d' else 7.); s=aper(c,x,y,rad)*1e3; s,e=baseline(c['v'],s,FREE[sp],1)
        k=(c['v']>2000)&(c['v']<3300); v=c['v'][k]; s=s[k]; n=len(v)//nb*nb
        vv=v[:n].reshape(-1,nb).mean(1); ss=s[:n].reshape(-1,nb).mean(1)
        # empirical noise of nb-channel means
        sm=np.convolve(s_all:=np.array(baseline(c['v'],aper(c,x,y,rad)*1e3,FREE[sp],1)[0]),np.ones(nb)/nb,'valid'); m=np.zeros(len(c['v']),bool)
        for a,b in FREE[sp]: m|=(c['v']>=a)&(c['v']<=b)
        mv=np.convolve(m.astype(float),np.ones(nb)/nb,'valid')>0.999; es=robust_sigma(sm[mv])
        out[sp]=(vv,ss,es)
    return out
def model(p,v,shared=True):
    vA,sA,vB,sB,hA,hB,oA,oB=p
    g=lambda a,mu,s: a*np.exp(-0.5*((v-mu)/s)**2)
    return g(hA,vA,sA)+g(hB,vB,sB),g(oA,vA,sA)+g(oB,vB,sB)
def fit(d,p0=(2800,60,2570,110,1.5,1.5,3.0,1.0),shared=True):
    (vh,sh,eh),(vo,so,eo)=d['HCN'],d['HCO+']
    def res(p):
        mh,_=model(p,vh); _,mo=model(p,vo); return np.concatenate([(sh-mh)/eh,(so-mo)/eo])
    lo=[2700,25,2300,60,0,0,0,0]; hi=[2900,120,2750,300,50,50,50,50]
    r=least_squares(res,p0,bounds=(lo,hi)); J=r.jac; cov=np.linalg.pinv(J.T@J); 
    return r.x,np.sqrt(np.diag(cov)),np.sum(r.fun**2),len(r.fun)-len(p0),cov
rows=[]
for ous in ('Xb4d','Xb4b'):
  for nm,(x,y,rad) in {'Soutflow line-integrated centroid, 2"':(-0.72,-4.79,1.0),'M18 South-OF-Blue, 2"':(-0.65,-4.83,1.0),'Soutflow wing peak, 2"':(-0.15,-4.5,1.0),'Soutflow wing peak, 1" (point)':(-0.15,-4.5,0.5)}.items():
    d=data(x,y,rad,ous); p,e,chi2,dof,cov=fit(d)
    RA=p[4]/p[6]*cT; RB=p[5]/p[7]*cT
    # MC for ratio errors from the covariance (amplitudes only, positive)
    S=rng.multivariate_normal(p,cov,5000); S=S[(S[:,6]>0)&(S[:,7]>0)]; rA=S[:,4]/S[:,6]*cT; rB=S[:,5]/S[:,7]*cT
    fB_h=p[5]*p[3]/(p[4]*p[1]+p[5]*p[3]); fB_o=p[7]*p[3]/(p[6]*p[1]+p[7]*p[3])
    # one-component fit for comparison
    def res1(q):
        (vh,sh,eh),(vo,so,eo)=d['HCN'],d['HCO+']; return np.concatenate([(sh-q[2]*np.exp(-0.5*((vh-q[0])/q[1])**2))/eh,(so-q[3]*np.exp(-0.5*((vo-q[0])/q[1])**2))/eo])
    r1=least_squares(res1,[2750,100,2,3]); chi1=np.sum(r1.fun**2)
    rows.append(dict(ous=ous,ap=nm,vA=p[0],e_vA=e[0],sA=p[1],e_sA=e[1],vB=p[2],e_vB=e[2],sB=p[3],e_sB=e[3],R_A=RA,R_A_lo=np.percentile(rA,16),R_A_hi=np.percentile(rA,84),R_B=RB,R_B_lo=np.percentile(rB,16),R_B_hi=np.percentile(rB,84),
         frac_broad_HCN=fB_h,frac_broad_HCO=fB_o,chi2=chi2,dof=dof,chi2_1comp=chi1,noise_HCN=d['HCN'][2],noise_HCO=d['HCO+'][2]))
D=pd.DataFrame(rows); D.round(3).to_csv(OUT+'twocomp.csv',index=False); pd.set_option('display.width',250); pr(D.round(2).T.to_string())
# ---- independent fits per line (free kinematics): integrated flux of the broad component and its ratio
pr('\n== per-line two-Gaussian fits (free v, sigma per line); narrow component constrained to 2790-2830 km/s, sigma 20-70')
def fit1(v,s,e):
    g=lambda a,mu,sg: a*np.exp(-0.5*((v-mu)/sg)**2)
    def res(p): return (s-g(p[0],p[1],p[2])-g(p[3],p[4],p[5]))/e
    r=least_squares(res,[2,2805,35,1.5,2650,110],bounds=([0,2790,20,0,2300,60],[50,2830,70,50,2800,350]))
    cov=np.linalg.pinv(r.jac.T@r.jac); return r.x,np.sqrt(np.diag(cov)),np.sum(r.fun**2),len(v)-6
rows=[]
for ous in ('Xb4d','Xb4b'):
  for nm,(x,y,rad) in {'Soutflow line-integrated centroid, 2"':(-0.72,-4.79,1.0),'M18 South-OF-Blue, 2"':(-0.65,-4.83,1.0)}.items():
    d=data(x,y,rad,ous); P={}
    for sp in ('HCN','HCO+'):
        v,s,e=d[sp]; p,pe,chi,dof=fit1(v,s,e); P[sp]=(p,pe,chi,dof)
        FA=p[0]*p[2]*np.sqrt(2*np.pi)/1e3*1; FB=p[3]*p[5]*np.sqrt(2*np.pi)/1e3
        pr(f'{ous} {nm} {sp}: narrow v={p[1]:.0f} sig={p[2]:.0f} amp={p[0]:.2f}mJy; broad v={p[4]:.0f}+-{pe[4]:.0f} sig={p[5]:.0f}+-{pe[5]:.0f} amp={p[3]:.2f}+-{pe[3]:.2f}mJy; chi2/dof {chi:.0f}/{dof}')
        rows.append(dict(ous=ous,ap=nm,line=sp,vN=p[1],sN=p[2],aN=p[0],vB=p[4],sB=p[5],aB=p[3],FB_Jy=FB*abs(np.median(np.diff(v)))*0+FB*1.0,chi2=chi,dof=dof))
    # integrated broad flux (mJy km/s ->): amp*sigma*sqrt(2pi) in mJy km/s
    FBH=P['HCN'][0][3]*P['HCN'][0][5]*np.sqrt(2*np.pi); FBO=P['HCO+'][0][3]*P['HCO+'][0][5]*np.sqrt(2*np.pi)
    FNH=P['HCN'][0][0]*P['HCN'][0][2]*np.sqrt(2*np.pi); FNO=P['HCO+'][0][0]*P['HCO+'][0][2]*np.sqrt(2*np.pi)
    pr(f'   broad-component flux ratio HCN/HCO+ (Tb) = {FBH/FBO*cT:.2f}; narrow-component = {FNH/FNO*cT:.2f}; broad fraction of the line flux: HCN {FBH/(FBH+FNH):.2f}, HCO+ {FBO/(FBO+FNO):.2f}')
pd.DataFrame(rows).round(3).to_csv(OUT+'twocomp_perline.csv',index=False)
