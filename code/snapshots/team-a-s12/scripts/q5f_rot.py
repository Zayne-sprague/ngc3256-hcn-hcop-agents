# Q5f: rotation curves along the kinematic axes (pixel-fit centroids binned by projected offset), enclosed dynamical masses and virial masses.
# python3 -I scripts/q5f_rot.py -> out/q5/q5f.log, q5f_rot.csv, q5f_rot.png
import sys; sys.path.insert(0,'/workspace/scripts')
from q5lib import *
import pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
Z=np.load(O+'q5e_pix.npz'); log=open(O+'q5f.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
G=4.30091e-3   # pc (km/s)^2 / Msun
SRC={'S':(-0.15,-2.60),'N':(0.0,2.60)}
INC={'N':(30,(20,40)),'S':(80,(70,90))}   # Sakamoto+2014 (ApJ 797,90) i_N~30 deg, i_S~80 deg (70-90); assumptions
PAK={'N':75.,'S':100.}                    # kinematic axes: measured gradient directions (q5e_grad.csv): N ~80 (77-85), S ~120 (108-124) at R<=1"; Sakamoto: N 75, S 90
rows=[]; fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for ia,src in enumerate(('N','S')):
    x0,y0=SRC[src]
    for nm in ('HCN','HCO+'):
        P=Z[f'{nm}_{src}']
        for pa in (PAK[src]-15,PAK[src],PAK[src]+15):
            t=np.deg2rad(pa); ex,ey=np.sin(t),np.cos(t)   # unit vector along axis (E of N)
            X=(P[:,0]-x0)*ex+(P[:,1]-y0)*ey; Yp=-(P[:,0]-x0)*ey+(P[:,1]-y0)*ex
            m=np.abs(Yp)<0.25
            for R in (0.25,0.5,0.75,1.0,1.5):
                # intensity-weighted mean velocity in bins [R-0.125,R+0.125] on each side (signed X)
                vp=[];vn=[]
                for sgn,l in ((1,vp),(-1,vn)):
                    mm=m&(sgn*X>R-0.125)&(sgn*X<R+0.125)
                    if mm.sum()<3: l.append(np.nan); continue
                    w=P[mm,2]; l.append(np.average(P[mm,3],weights=w)); 
                v_p,v_n=vp[0],vn[0]
                if np.isnan(v_p) or np.isnan(v_n): continue
                vrot=abs(v_p-v_n)/2; vc=(v_p+v_n)/2
                # sigma of the bins (intensity weighted)
                sg=np.average(P[m&(np.abs(X)>R-0.125)&(np.abs(X)<R+0.125),5],weights=P[m&(np.abs(X)>R-0.125)&(np.abs(X)<R+0.125),2])
                i0,(ilo,ihi)=INC[src]
                Rpc=R*PC
                M=lambda i: Rpc*(vrot/np.sin(np.deg2rad(i)))**2/G
                rows.append(dict(src=src,line=nm,pa=pa,R_arcsec=R,R_pc=Rpc,v_plus=v_p,v_minus=v_n,vrot_proj=vrot,vcen=vc,sigma=sg,i=i0,vrot_deproj=vrot/np.sin(np.deg2rad(i0)),Mdyn=M(i0),Mdyn_ilo=M(ilo),Mdyn_ihi=M(ihi),
                    Mpress=Rpc*(vrot**2/np.sin(np.deg2rad(i0))**2+2*sg**2)/G))
df=pd.DataFrame(rows); df.to_csv(O+'q5f_rot.csv',index=False); pd.set_option('display.width',220)
pr(df[df.pa==df.pa.where(df.src=='N',0).max() if False else True].round(2).to_string() if False else '')
for src in ('N','S'):
    d_=df[(df.src==src)&(df.pa==PAK[src])]
    pr(f'\n# {src}: rotation (kinematic axis PA {PAK[src]}, slit half-width 0.25", bins +-0.125"), i={INC[src][0]} deg; Mdyn=R vrot^2/(G sin^2 i) [enclosed, circular]; Mpress adds 2 sigma^2 (asymmetric-drift-like pressure support)')
    pr(d_[['line','R_arcsec','R_pc','v_plus','v_minus','vrot_proj','vcen','sigma','vrot_deproj','Mdyn','Mdyn_ilo','Mdyn_ihi','Mpress']].round(1).to_string(float_format=lambda x:f'{x:.3g}'))
    pr('PA variation (+-15 deg) of vrot_proj at R=0.5, 1.0:')
    pr(df[(df.src==src)&(df.R_arcsec.isin([0.5,1.0]))][['line','pa','R_arcsec','vrot_proj','vcen']].round(1).to_string())
    for nm,c in (('HCN','C0'),('HCO+','C1')):
        s_=df[(df.src==src)&(df.pa==PAK[src])&(df.line==nm)]; ax[ia].plot(s_.R_arcsec,s_.vrot_proj,'o-',c=c,label=nm); ax[ia].set_title(src+' projected rotation (half velocity difference)'); ax[ia].legend(); ax[ia].set_xlabel('R (")'); ax[ia].set_ylabel('km/s')
plt.savefig(O+'q5f_rot.png',dpi=80)
