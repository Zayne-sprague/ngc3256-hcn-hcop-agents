# Q9: velocity-resolved HCN/HCO+ 3-2 ratio (Jy ratio at a common beam; x1.0126 for Tb) in annuli about N, S, Bridge.
# Xb4d @1.05" circular and Xb4b @0.30"; 40-km/s bins (9 channels); noise from random line-free bins/positions; constancy tests; widths.
# python3 -I q9_vel.py  (needs out/q7/cube_*.npz from q7_cache.py)
import sys; sys.path.insert(0,'/workspace/scripts')
from q7_lib import *
import pandas as pd
rng=np.random.default_rng(9)
TBF=(265.886434/267.557626)**-2     # Tb_HCN/Tb_HCO+ = (S_HCN/S_HCO+)*(nu_O/nu_H)^2
NB=9
cb={(o,l):Cube(o,l) for o in ('Xb4b','Xb4d') for l in ('HCN','HCO+')}
CASES={'Xb4d@1.05':('Xb4d',(1.05,1.05,0.),{'N':[(0,.5),(.5,1),(1,2)],'S':[(0,.5),(.5,1),(1,2)],'Bridge':[(0,1)]}),
       'Xb4b@0.30':('Xb4b',(0.30,0.30,0.),{'N':[(0,.25),(.25,.5),(.5,1)],'S':[(0,.25),(.25,.5),(.5,1)]})}
REG={'N':(0.0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15)}
def smcube(c,tb):
    out=np.empty_like(c.cube)
    k=None
    for i in range(c.cube.shape[0]): out[i]=smooth(c.cube[i],c,tb)
    return out
def apspec(c,cube,x0,y0,r1,r2,tb):
    m=(np.hypot(c.X-x0,c.Y-y0)>=r1)&(np.hypot(c.X-x0,c.Y-y0)<r2)
    return cube[:,m].sum(1)/(OM(tb)/c.pix**2)     # Jy per channel (flux density)
def rebin(v,s,nb=NB):
    n=len(v)//nb*nb; return v[:n].reshape(-1,nb).mean(1), s[:n].reshape(-1,nb).sum(1)*(abs(np.median(np.diff(v))))   # Jy km/s per bin
rows=[];prof=[];tests=[]
for case,(ous,tb,regs) in CASES.items():
    H,O=cb[(ous,'HCN')],cb[(ous,'HCO+')]
    sH,sO=smcube(H,tb),smcube(O,tb)
    imO=smooth(O.m0(-350,350),O,tb)
    for reg,ann in regs.items():
        x0,y0=peak_near(O,imO,*REG[reg])
        for r1,r2 in ann:
            a=apspec(H,sH,x0,y0,r1,r2,tb); b=apspec(O,sO,x0,y0,r1,r2,tb)
            vb,A=rebin(H.v,a); vb2,B=rebin(O.v,b)
            # velocity grids differ by sub-channel shift between HCN and HCO+ cubes: interpolate HCO+ to HCN velocity grid
            # (cubes are in the same frame via their own rest freqs; both channel widths equal)
            # noise: random centres, line-free bins (|v|>550), same aperture
            nh=[];no=[]
            sel=np.where((np.abs(H.v)>550))[0]; sel2=np.where((np.abs(O.v)>550))[0]
            for _ in range(400):
                xr,yr=rng.uniform(-5,5,2)
                s1=apspec(H,sH,xr,yr,r1,r2,tb); s2=apspec(O,sO,xr,yr,r1,r2,tb)
                for ss,vv,sl,lst in ((s1,H.v,sel,nh),(s2,O.v,sel2,no)):
                    i0=rng.integers(0,len(sl)-NB)
                    idx=sl[i0:i0+NB]
                    if np.all(np.diff(idx)==1): lst.append(ss[idx].sum()*abs(np.median(np.diff(vv))))
            eH=1.4826*np.median(np.abs(np.array(nh)-np.median(nh))); eO=1.4826*np.median(np.abs(np.array(no)-np.median(no)))
            vO_on_H=np.interp(vb,vb2[::-1],B[::-1]) if vb2[0]>vb2[-1] else np.interp(vb,vb2,B)
            for v,Ah,Bo in zip(vb,A,vO_on_H):
                if abs(v)>450: continue
                R=Ah/Bo if Bo!=0 else np.nan
                ok=(Ah>3*eH)&(Bo>3*eO)
                rows.append(dict(case=case,reg=reg,r1=r1,r2=r2,x0=round(x0,2),y0=round(y0,2),v=round(v,1),HCN=Ah,eHCN=eH,HCO=Bo,eHCO=eO,R=R if ok else np.nan,
                                 eR=abs(R)*np.hypot(eH/Ah,eO/Bo) if ok else np.nan,snrH=Ah/eH,snrO=Bo/eO,ok=ok,Rlim=(3*eH/Bo if (Bo>3*eO and not Ah>3*eH) else np.nan)))
            # profile properties from native-res channel spectra within +-400
            for ln,s,c in (('HCN',a,H),('HCO+',b,O)):
                m=np.abs(c.v)<400; v=c.v[m]; y=s[m]; pk=y.max(); half=np.where(y>pk/2)[0]
                w=np.clip(y,0,None); mom1=(v*w).sum()/w.sum(); sig=np.sqrt(((v-mom1)**2*w).sum()/w.sum())
                # FWHM by linear interpolation of half-max crossings
                o=np.argsort(v); vv,yy=v[o],y[o]; ip=yy.argmax(); 
                lo=ip
                while lo>0 and yy[lo]>pk/2: lo-=1
                hi=ip
                while hi<len(yy)-1 and yy[hi]>pk/2: hi+=1
                vlo=np.interp(pk/2,[yy[lo],yy[lo+1]],[vv[lo],vv[lo+1]]); vhi=np.interp(pk/2,[yy[hi],yy[hi-1]],[vv[hi],vv[hi-1]])
                prof.append(dict(case=case,reg=reg,r1=r1,r2=r2,line=ln,peak_Jy=pk,v_peak=vv[ip],centroid=mom1,fwhm=vhi-vlo,sigma_mom2=sig))
df=pd.DataFrame(rows); pr=pd.DataFrame(prof)
df.to_csv('/workspace/out/q9/q9_vel_ratio.csv',index=False); pr.to_csv('/workspace/out/q9/q9_profiles.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500)
print('Ratio = Jy ratio (multiply by %.4f for Tb); 40 km/s bins; S/N>3 in both lines'%TBF)
for (case,reg,r1,r2),g in df.groupby(['case','reg','r1','r2'],sort=False):
    g=g.sort_values('v'); print('\n',case,reg,f'r={r1}-{r2}"')
    print(g[['v','HCN','eHCN','HCO','eHCO','R','eR','Rlim']].round(3).to_string(index=False))
    gg=g[g.ok]
    if len(gg)>=3:
        w=1/gg.eR**2; Rm=(gg.R*w).sum()/w.sum(); chi=((gg.R-Rm)**2*w).sum()
        p=np.polyfit(gg.v,gg.R,1,w=1/gg.eR); cov=np.linalg.inv(np.vstack([gg.v/gg.eR,1/gg.eR]).dot(np.vstack([gg.v/gg.eR,1/gg.eR]).T))
        core=gg[np.abs(gg.v)<=60]; wing=gg[np.abs(gg.v)>60]
        print(f'  nbins={len(gg)} weighted mean R={Rm:.3f}  chi2/dof(const)={chi:.1f}/{len(gg)-1}  slope={p[0]*100:.3f}+-{np.sqrt(cov[0,0])*100:.3f} per 100 km/s')
        tests.append(dict(case=case,reg=reg,r1=r1,r2=r2,nbins=len(gg),Rmean=Rm,chi2=chi,dof=len(gg)-1,slope100=p[0]*100,eslope100=np.sqrt(cov[0,0])*100,
                          Rcore=np.average(core.R,weights=1/core.eR**2) if len(core) else np.nan,Rwing_blue=np.average(gg[gg.v<-60].R,weights=1/gg[gg.v<-60].eR**2) if (gg.v<-60).any() else np.nan,
                          Rwing_red=np.average(gg[gg.v>60].R,weights=1/gg[gg.v>60].eR**2) if (gg.v>60).any() else np.nan))
pd.DataFrame(tests).to_csv('/workspace/out/q9/q9_constancy.csv',index=False)
print('\nPROFILES'); print(pr.round(2).to_string(index=False))
print('\nTESTS'); print(pd.DataFrame(tests).round(3).to_string(index=False))
