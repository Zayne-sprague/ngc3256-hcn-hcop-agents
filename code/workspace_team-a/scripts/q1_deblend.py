# Q1 fix 1: deblend coarse-beam maps at their NATIVE elliptical beams (Band 5 2-1 9.4x6.2", Band 7 4-3 / CS 7-6 4.3x2.8", Xb4f 3-2 6.1x3.7" as a recovery test; the beams are
# elongated E-W so the N-S direction (5.1" apart) is better resolved than a circular smoothing to 6.5/10" would give) into N, S, Bridge, Outer clumps and diffuse
# components, using the Xb4d/Xb4b 3-2 moment-0 maps (~1" / 0.24") split into components and smoothed to the coarse beam as spatial templates.
# Model: y(pix) = sum_c a_c T_c(pix), all in Tb (K km/s); a_c = Tb(line)/Tb(3-2) of component c (non-negative, NNLS), weights sigma_eff^2 = sigma_null^2 + (0.03 y)^2
# (3% template-mismatch floor; without it chi2_r ~100-250 and the formal errors are meaningless). Errors: block bootstrap (blocks = one beam) in quadrature
# with the spread over an ensemble of template variants (component radius, Xb4d/Xb4b templates, Wclump assignment, +missing-flux template).
# Inject-and-recover uses truth built from one template set and recovery with another (template mismatch is the dominant systematic).
import sys,warnings,json,itertools; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
from scipy.ndimage import map_coordinates
from scipy.optimize import nnls
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15),'NEclump':(6.75,3.9),'Eclump':(5.25,-2.85),'Wclump':(-1.8,-2.4),'Sext':(-1.05,-5.4),'SWarm':(-3.75,-4.5)}
GROUP={'N':'N','S':'S','Bridge':'Bridge','NEclump':'Outer','Eclump':'Outer','Sext':'Outer','SWarm':'Outer','Wclump':'Outer'}
GROUP_WS={**GROUP,'Wclump':'S'}
COMPS=['N','S','Bridge','Outer','diffuse']
NU=lambda key: LINES[key][2]/(1+VSYS/C)
CACHE={}
def get(key,crop=None):
    if key not in CACHE: CACHE[key]=load(key,crop=crop)
    return CACHE[key]
def resample(img,ds,dc):
    w=WCS(ds['hd']).celestial
    px,py=w.world_to_pixel_values(RA0+dc['dx']/3600/np.cos(np.deg2rad(DE0)),DE0+dc['dy']/3600)
    return map_coordinates(img,[py-ds['y0'],px-ds['x0']],order=1,mode='constant',cval=0.)
TCACHE={}
def TbK(freq,b): return 1.222e6/(freq**2*b[0]*b[1])       # K per Jy/beam, elliptical Gaussian beam (maj,min arcsec)
def smooth_to(d,tb,arr):
    k=kernel(d['hd'],tb,d['beam']); 
    if k is None: raise ValueError('target beam smaller than native')
    return fftconvolve(np.nan_to_num(arr),k,mode='same')*OM(tb)/OM(d['beam'])
def templates(sp,tsrc,fw,dc_key,rad,grp,missing):
    key=(sp,tsrc,fw,dc_key,rad,id(grp),missing)
    if key in TCACHE: return TCACHE[key]
    dc=get(dc_key); tb=dc['beam']; tk=f'{sp} 3-2@{tsrc}'; ds=get(tk,crop=14. if tsrc=='Xb4b' else None)
    names=list(POS); dist=np.array([np.hypot(ds['dx']-POS[n][0],ds['dy']-POS[n][1]) for n in names]); near=dist.argmin(0); ok=dist.min(0)<rad
    m0=np.nan_to_num(ds['m0']); noise=np.nanstd([n for _,n in ds['nulls']]) if ds['nulls'] else 0
    k=TbK(NU(tk),tb); T={}
    for c in COMPS[:-1]:
        msk=np.zeros(m0.shape,bool)
        for i,n in enumerate(names):
            if grp[n]==c: msk|=(near==i)&ok
        T[c]=resample(smooth_to(ds,tb,np.where(msk,m0,0.)),ds,dc)*k
    # diffuse: Xb4d emission outside the compact apertures, clipped at 3 sigma of the native-resolution null noise (avoids fitting noise)
    sg=np.median([np.nanstd(n) for _,n in ds['nulls']]); dif=np.where((~ok)&(m0>3*sg),m0,0.)
    T['diffuse']=resample(smooth_to(ds,tb,dif),ds,dc)*k
    if missing:
        d4f=get(f'{sp} 3-2@Xb4f'); s4f=resample(smooth_to(d4f,tb,d4f['m0']),d4f,dc)*TbK(NU(f'{sp} 3-2@Xb4f'),tb)
        T['missing']=np.clip(s4f-sum(T.values()),0,None)
    TCACHE[key]=T; return T
def fitnn(y,sig,T,mask,floor=0.03):
    names=list(T); A=np.stack([T[n][mask] for n in names],1); s=np.sqrt(sig[mask]**2+(floor*y[mask])**2); a,_=nnls(A/s[:,None],y[mask]/s)
    res=(y[mask]-A@a)/s; return names,a,float(np.mean(res**2)),(A,y[mask],s)
def bootstrap(pack,fw,pix,mask,nb=60,seed=0):
    A,y,s=pack; rng=np.random.default_rng(seed)
    # block ids from pixel coordinates of masked pixels
    return None
def run(lkey,sp,tsrc,fw,rfit,rad,grp,missing,boot=0,seed=0):
    dc=get(lkey); sm=np.nan_to_num(dc['m0']); prof=noise_profile(dc,None); r=np.hypot(dc['dx'],dc['dy'])
    sig=np.vectorize(lambda q:sigma_at(prof,q))(r); k=TbK(NU(lkey),dc['beam']); y=sm*k; s=sig*k
    T=templates(sp,tsrc,fw,lkey,rad,grp,missing); mask=(r<rfit)&np.isfinite(y)&(sm!=0)
    names,a,chi2,pack=fitnn(y,s,T,mask)
    eb=np.full(len(a),np.nan)
    if boot:
        rng=np.random.default_rng(seed); A,yy,ss=pack; X=dc['dx'][mask]; Y=dc['dy'][mask]; bs=dc['beam'][0]
        bid=(np.floor(X/bs).astype(int)*1000+np.floor(Y/bs).astype(int)); ub=np.unique(bid); sel=[np.where(bid==u)[0] for u in ub]
        fits=[]
        for _ in range(boot):
            idx=np.concatenate([sel[i] for i in rng.integers(0,len(sel),len(sel))])
            fits.append(nnls(A[idx]/ss[idx,None],yy[idx]/ss[idx])[0])
        eb=np.std(fits,0)
    return dict(names=names,a=a,eb=eb,chi2=chi2,T=T,y=y,s=s,mask=mask,dc=dc,k=k)
SPECS=[('HCN 2-1','HCN',10.0,22.),('HCO+ 2-1','HCO+',10.0,22.),('HCN 3-2@Xb4f','HCN',10.0,22.),('HCO+ 3-2@Xb4f','HCO+',10.0,22.),
       ('HCN 4-3','HCN',6.5,14.),('HCO+ 4-3','HCO+',6.5,14.),('CS 7-6','HCO+',6.5,14.),('CS 7-6','HCN',6.5,14.)]
VARS=[dict(tsrc=t,rad=r,grp=g,missing=m,tag=f'{t}/r{r}/{"Wc->S" if g is GROUP_WS else "Wc->Outer"}/{"+miss" if m else "nomiss"}')
      for t,r,g,m in itertools.product(('Xb4d','Xb4b'),(1.0,1.5,2.0,2.5),(GROUP,GROUP_WS),(False,True))]
BASE=dict(tsrc='Xb4d',rad=1.5,grp=GROUP,missing=False)
rows=[]; ENS={}
for lkey,sp,fw,rf in SPECS:
    for v in VARS:
        if v['missing'] and 'Xb4f' in lkey: continue
        try: R=run(lkey,sp,v['tsrc'],fw,rf,v['rad'],v['grp'],v['missing'],boot=40 if (v['tag']==f"Xb4d/r1.5/Wc->Outer/nomiss") else 0)
        except ValueError: continue
        ENS[(lkey,sp,v['tag'])]=R
        for i,n in enumerate(R['names']):
            rows.append(dict(line=lkey,template=sp,fw=fw,variant=v['tag'],comp=n,a=R['a'][i],a_boot=R['eb'][i],chi2r=R['chi2'],Tsum=float(R['T'][n].sum())))
    b=ENS[(lkey,sp,'Xb4d/r1.5/Wc->Outer/nomiss')]; print(lkey,sp,'base chi2r %.2f'%b['chi2'],' '.join(f"{n}:{b['a'][i]:.3f}±{b['eb'][i]:.3f}" for i,n in enumerate(b['names'])),flush=True)
A=pd.DataFrame(rows); A.to_csv('/workspace/out/q1/q1_deblend_amplitudes.csv',index=False)

# ---- ensemble statistics of derived quantities ----
def get_a(l,sp,tag,c):
    R=ENS.get((l,sp,tag))
    if R is None or c not in R['names']: return None
    i=R['names'].index(c); return R['a'][i],R['eb'][i],float(R['T'][c].sum())
TAGS=[v['tag'] for v in VARS if not v['missing']]; TAGSM=[v['tag'] for v in VARS if v['missing'] and (('HCN 2-1','HCN',v['tag']) in ENS)]
BT='Xb4d/r1.5/Wc->Outer/nomiss'
def summarize(vals,boots):
    vals=np.array([x for x in vals if np.isfinite(x)]); 
    if len(vals)==0: return dict(value=np.nan)
    med=np.median(vals); lo,hi=np.percentile(vals,[16,84]); eb=np.nanmedian(boots) if len(boots) else np.nan
    return dict(value=med,spread_lo=lo,spread_hi=hi,boot_err=eb,n_var=len(vals))
out=[]
def q_ladder(sp,line,c,tags,kind):
    v=[];b=[]
    for t in tags:
        x=get_a(line,sp,t,c); 
        if x is None or x[0]<=0: v.append(np.nan); continue
        a,e,_=x; v.append(1/a if kind=='inv' else a); b.append(e/a**2 if kind=='inv' else e)
    return summarize(v,b)
for sel,tags in (('nomissing',TAGS),('+missing-template',TAGSM)):
    for c in ['N','S','Bridge','Outer','diffuse']:
        for sp in ('HCN','HCO+'):
            for nm,line,kind in (('3-2/2-1 (Tb)',f'{sp} 2-1','inv'),('4-3/3-2 (Tb)',f'{sp} 4-3','dir'),('3-2(Xb4f)/3-2(Xb4d) (recovery)',f'{sp} 3-2@Xb4f','dir')):
                r=q_ladder(sp,line,c,tags,kind); r.update(model=sel,comp=c,species=sp,quantity=nm); out.append(r)
        for sp in ('HCN','HCO+'):
            v=[];
            for t in tags:
                x=get_a('CS 7-6',sp,t,c); v.append(x[0] if x else np.nan)
            r=summarize(v,[]); r.update(model=sel,comp=c,species=sp,quantity=f'CS 7-6/{sp} 3-2 (Tb)'); out.append(r)
        for J,fw,lh,lo in (('2-1',10,'HCN 2-1','HCO+ 2-1'),('4-3',6.5,'HCN 4-3','HCO+ 4-3'),('3-2 (via Xb4f)',10,'HCN 3-2@Xb4f','HCO+ 3-2@Xb4f')):
            v=[];b=[]
            for t in tags:
                x,y=get_a(lh,'HCN',t,c),get_a(lo,'HCO+',t,c)
                if x and y and y[0]>0 and x[0]>=0: 
                    v.append(x[0]*x[2]/(y[0]*y[2])); 
                    if t==BT and np.isfinite(x[1]) and x[0]>0: b.append(v[-1]*np.hypot(x[1]/x[0],y[1]/y[0]))
                else: v.append(np.nan)
            r=summarize(v,b); r.update(model=sel,comp=c,species='HCN/HCO+',quantity=f'HCN/HCO+ {J} (Tb, component-integrated)'); out.append(r)
        # 1" 3-2 ratio of the component from the templates themselves (matched 3-2)
        v=[]
        for t in tags:
            x=ENS[('HCN 2-1','HCN',t)]['T'][c].sum(); y=ENS[('HCO+ 2-1','HCO+',t)]['T'][c].sum(); v.append(x/y)  # 10"x6" Band-5 beam (Tb)
        r=summarize(v,[]); r.update(model=sel,comp=c,species='HCN/HCO+',quantity='HCN/HCO+ 3-2 (Tb) from 1" templates, same aperture sets'); out.append(r)
E=pd.DataFrame(out); E.to_csv('/workspace/out/q1/q1_deblend_summary.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',500)
print(E[E.model=='nomissing'].round(3).to_string())
# N/S contrast per quantity (ensemble): ratio of component values within each variant
ns=[]
for qn,fn in (('HCN/HCO+ 2-1',None),):
    pass
for J,fw,lh,lo in (('2-1',10,'HCN 2-1','HCO+ 2-1'),('4-3',6.5,'HCN 4-3','HCO+ 4-3')):
    for sel,tags in (('nomissing',TAGS),('+missing-template',TAGSM)):
        v=[]
        for t in tags:
            xn,yn,xs,ys=get_a(lh,'HCN',t,'N'),get_a(lo,'HCO+',t,'N'),get_a(lh,'HCN',t,'S'),get_a(lo,'HCO+',t,'S')
            if None in (xn,yn,xs,ys): continue
            if min(yn[0],ys[0],xs[0])>0: v.append((xn[0]*xn[2]/(yn[0]*yn[2]))/(xs[0]*xs[2]/(ys[0]*ys[2])))
        v=np.array(v)
        if len(v)==0: continue
        ns.append(dict(quantity=f'(HCN/HCO+ {J})_N / (HCN/HCO+ {J})_S',model=sel,median=np.median(v),p16=np.percentile(v,16),p84=np.percentile(v,84),min=v.min(),max=v.max(),n=len(v)))
pd.DataFrame(ns).to_csv('/workspace/out/q1/q1_deblend_NS_contrast.csv',index=False); print(pd.DataFrame(ns).round(3).to_string())


# ---- ladders against the more complete Xb4f 3-2 (a_L / a_3-2(Xb4f); a(3-2 Xb4f) = Tb(Xb4f 3-2)/Tb(Xb4d 3-2) of the component) ----
pr=[]
for c in ['N','S','Bridge','Outer']:
    for sp in ('HCN','HCO+'):
        for nm,l1,l2 in (('4-3/3-2 [Xb4f 3-2] (Tb)',f'{sp} 4-3',f'{sp} 3-2@Xb4f'),('3-2 [Xb4f]/2-1 (Tb)',f'{sp} 3-2@Xb4f',f'{sp} 2-1')):
            v=[]
            for t in TAGS:
                x,y=get_a(l1,sp,t,c),get_a(l2,sp,t,c)
                v.append(x[0]/y[0] if (x and y and y[0]>0 and x[0]>0) else np.nan)
            r=summarize(v,[]); r.update(comp=c,species=sp,quantity=nm); pr.append(r)
PR=pd.DataFrame(pr); PR.to_csv('/workspace/out/q1/q1_deblend_ladders_xb4f.csv',index=False); print(PR.round(3).to_string())

# ---- inject-and-recover with template mismatch: truth from Xb4b/r2.0 templates (from its own fit, N/S amplitude ratio forced to f), recovered with the Xb4d/r1.5 base set ----
inj=[]
TRUTH=dict(tsrc='Xb4b',rad=2.0,grp=GROUP_WS,missing=False); TR_TAG='Xb4b/r2.0/Wc->S/nomiss'
for lkey,sp,fw,rf in (('HCN 2-1','HCN',10.0,22.),('HCO+ 2-1','HCO+',10.0,22.),('HCN 4-3','HCN',6.5,14.),('HCO+ 4-3','HCO+',6.5,14.),('CS 7-6','HCO+',6.5,14.)):
    Rt=ENS[(lkey,sp,TR_TAG)]; names=Rt['names']; at=Rt['a'].copy(); gm=np.sqrt(max(at[0],1e-6)*max(at[1],1e-6)) if at[0]>0 and at[1]>0 else max(at[:2].max(),1e-3)
    dc=Rt['dc']; k=Rt['k']; nulls=[np.nan_to_num(nm)*k for _,nm in dc['nulls']]
    for f in (1.0,1.25,1.5,2.0,3.0):
        a2=at.copy(); a2[0]=gm*np.sqrt(f); a2[1]=gm/np.sqrt(f)
        model=sum(a2[i]*Rt['T'][n] for i,n in enumerate(names)); rec=[]
        for z in nulls:
            for fl in (lambda q:q,lambda q:-q,lambda q:q[::-1],lambda q:-q[:,::-1]):
                zz=np.nan_to_num(fl(z)); yy=model+zz; Tb=ENS[(lkey,sp,BT)]['T']
                nn,aa,_,_=fitnn(yy,Rt['s'],Tb,Rt['mask']); 
                # compare in integrated flux terms (a*Tsum): true N/S flux ratio vs recovered
                rec.append((aa[0]*Tb['N'].sum())/(aa[1]*Tb['S'].sum()+1e-12))
        truth=(a2[0]*Rt['T']['N'].sum())/(a2[1]*Rt['T']['S'].sum())
        rec=np.array(rec); inj.append(dict(line=lkey,true_NS_flux_ratio=truth,f_amp=f,rec_median=np.median(rec),rec_p16=np.percentile(rec,16),rec_p84=np.percentile(rec,84),n=len(rec)))
I=pd.DataFrame(inj); I.to_csv('/workspace/out/q1/q1_deblend_inject.csv',index=False); print(I.round(3).to_string())

# ---- what the injection implies for the N/S contrast of HCN/HCO+ (ratio of the two species' recovered/true N/S flux ratios at no injected difference) ----
cb=[]
for J,(h,o) in (('2-1',('HCN 2-1','HCO+ 2-1')),('4-3',('HCN 4-3','HCO+ 4-3'))):
    for f in (1.0,1.5,2.0):
        a=I[(I.line==h)&(np.isclose(I.f_amp,f))].iloc[0]; b=I[(I.line==o)&(np.isclose(I.f_amp,f))].iloc[0]
        ra,rb=a.rec_median/a.true_NS_flux_ratio,b.rec_median/b.true_NS_flux_ratio
        cb.append(dict(J=J,injected_f_per_species=f,HCN_rec_over_true=ra,HCO_rec_over_true=rb,contrast_bias_factor=ra/rb))
CB=pd.DataFrame(cb); CB.to_csv('/workspace/out/q1/q1_deblend_contrast_bias.csv',index=False); print(CB.round(3).to_string())
