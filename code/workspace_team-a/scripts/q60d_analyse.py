"""Q60 step 3 (plain python): analyse uv-plane component spectra from q60c (out/q60/vis/*.npz).
Per EB/line/uvset/model: baseline (linear, line-free 470<|dv|<edge km/s), integral +-HW, curve of growth, N/S, HCN/HCO+.
Writes out/q60/uvfit.csv, out/q60/cog.csv. Run: python3 -I q60d_analyse.py"""
import numpy as np, glob, os, csv
C=299792.458; VC=2767.5
EBPROJ={'Xb7d4':'Xb47','X2cef':'Xb47','X3bb9':'Xb47','X178e':'Xb47','X11ec':'Xb47','X8300':'Xb47','X8d14':'Xb47','Xabe':'Xb47','X164':'Xb47',
 'Xdcb5':'00493','X70b8':'00493','X52e3':'00493','X5fca':'00223','Xdd41':'00223','X3cae':'01246-12m'}
def integ(v,F,E,hw,bl=(470.,1e4),vc=VC):
    d=v-vc; wb=(abs(d)>bl[0])&(abs(d)<bl[1]); wl=abs(d)<hw
    if wb.sum()<8 or wl.sum()<3: return np.nan,np.nan,wb.sum()
    p=np.polyfit(d[wb]/1e3,F[wb],1); r=F-np.polyval(p,d/1e3)
    dv=abs(np.median(np.diff(v))); I=r[wl].sum()*dv
    rms=np.std(r[wb]); return I,rms*dv*np.sqrt(wl.sum()),wb.sum()
rows=[];cog=[]
for f in sorted(glob.glob('/workspace/out/q60/vis/*.npz')):
    z=np.load(f); b=os.path.basename(f)[:-4]; eb=b.rsplit('_',2)[1]; ln=b.rsplit('_',1)[1]; proj=EBPROJ[eb]
    fr=z['freq']; v=(1-fr/float(z['rest']))*C; o=np.argsort(v); v=v[o]
    for k in z.keys():
        if not k.startswith('F_'): continue
        _,us,mod=k.split('_'); F=z[k][o]; E=z['E_'+us+'_'+mod][o]
        comps=['N','S','E'][:F.shape[1]]
        for hw in (250.,350.,500.):
            for ic,cn in enumerate(comps):
                I,e,nb=integ(v,F[:,ic],E[:,ic],hw)
                rows.append(dict(proj=proj,eb=eb,line=ln,uvset=us,model=mod,comp=cn,hw=hw,I=I,err=e,nbase=nb,vmin=v.min()-VC,vmax=v.max()-VC,r50=float(z['r50'])))
        if us=='all' and mod=='2':
            for ic,cn in enumerate(comps):
                for hw in (100,150,200,250,300,350,400,450,500):
                    I,e,nb=integ(v,F[:,ic],E[:,ic],float(hw)); cog.append(dict(proj=proj,eb=eb,line=ln,comp=cn,hw=hw,I=I,err=e))
for n,R in (('uvfit',rows),('cog',cog)):
    with open(f'/workspace/out/q60/{n}.csv','w') as fh:
        w=csv.DictWriter(fh,fieldnames=list(R[0].keys())); w.writeheader(); w.writerows(R)
print(len(rows),len(cog))
