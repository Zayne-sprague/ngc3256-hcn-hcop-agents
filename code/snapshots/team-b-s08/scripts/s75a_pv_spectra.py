# Q2 step b: PV cuts along the S-outflow axis and region spectra with 1- vs 2-component Gaussian decomposition.
# Axis: PA=185 deg E of N through S (s74a: wing annulus peaks at PA 180-210; s11b: PA_centroid 183-189); the same line passes through the Bridge and N (N is 5.2" N, 0.16" E of S).
# Data: Xb4d HCN/HCO+ 3-2 at native beam (0.94x0.74) and 1.05" circular; Xb4b smoothed to 0.5" and native 0.25".  pbcor, baseline per-pixel median(-900..-450), VSYS 2780 radio.
# Output: out/q2/s75a.txt, s75a.json, s75a_pv.npz, s75a_pv.png, s75a_spec.png
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
PA=185.0; S0=REGIONS['S'][:2]
def prep(ous,fw,half):
    o={}
    for mol,spw in (('HCN','spw27'),('HCO+','spw25')):
        c,v,ps,wc,x0,y0,b0=load_cube(ous,spw,F0[mol][3],half=half)
        if fw: c=conv(c,ps,b0,(fw,fw,0.0)); b=(fw,fw,0.0)
        else: b=b0
        c=c-np.median(c[(v>-900)&(v<-450)],axis=0)
        o[mol]=dict(c=c,v=v,dv=abs(np.median(np.diff(v))),ps=ps,n=c.shape[1],beam=b,ba=1.1331*b[0]*b[1]/ps**2,pxy=pxy_fn(wc,x0,y0))
    return o
SET={'Xb4d_1.05':prep('Xb4d',1.05,24.0),'Xb4d_nat':prep('Xb4d',None,24.0),'Xb4b_0.5':prep('Xb4b',0.5,10.0)}
def line_pos(base,off,pa=PA):
    th=np.deg2rad(pa); dE=off*np.sin(th); dN=off*np.cos(th)       # off>0 along PA (south for PA=185 -> use sign below)
    return offset_pos(base,dE,dN)
def pv(R,offs,width=1.0,pa=PA):
    """offs: arcsec along the cut, positive = NORTH of S (opposite to PA 185 direction). slit width averaged with 5 parallel strips."""
    th=np.deg2rad(pa); ex=np.array([np.sin(th),np.cos(th)])  # unit (E,N) along PA
    perp=np.array([ex[1],-ex[0]]); out=np.zeros((len(R['v']),len(offs)))
    pS=np.array(R['pxy'](*S0))
    for j,o in enumerate(offs):
        acc=0
        for w in np.linspace(-width/2,width/2,5):
            d=(-o)*ex+w*perp                      # (E,N) arcsec; north-positive o
            px=pS[0]-d[0]/R['ps']; py=pS[1]+d[1]/R['ps']       # x increases to the West (RA decreases)
            acc=acc+np.array([map_coordinates(R['c'][i],[[py],[px]],order=1)[0] for i in range(len(R['v']))])
        out[:,j]=acc/5
    return out
offs=np.arange(-7,7.01,0.25)
pvs={}
for k in ('Xb4d_1.05','Xb4b_0.5'):
    for mol in ('HCN','HCO+'): pvs[f'{k}_{mol}']=pv(SET[k][mol],offs,1.0 if k=='Xb4d_1.05' else 0.5)
np.savez(OUT+'/s75a_pv.npz',offs=offs,v_HCN=SET['Xb4d_1.05']['HCN']['v'],**{'v_HCO+':SET['Xb4d_1.05']['HCO+']['v']},**pvs)
# rms per channel in PV (line-free channels of the PV map)
fig,ax=plt.subplots(2,4,figsize=(20,9))
for r,k in enumerate(('Xb4d_1.05','Xb4b_0.5')):
    for j,mol in enumerate(('HCN','HCO+')):
        M=pvs[f'{k}_{mol}']; v=SET[k][mol]['v']; lf=(v>-900)&(v<-450); s=M[lf].std()
        P(f'PV {k} {mol}: rms per channel (line-free, slit-averaged) {s*1e3:.2f} mJy/beam ; peak {M.max()*1e3:.1f} mJy/beam')
        a=ax[r,2*j]; a.imshow(M,origin='lower',aspect='auto',extent=[offs[0],offs[-1],v[0],v[-1]],cmap='viridis',vmin=-3*s,vmax=M.max())
        a.contour(offs,v,M,levels=[3*s,5*s,10*s],colors='w',linewidths=.5); a.set_ylim(-700,600); a.set_xlabel('offset N of S along PA185 (")'); a.set_ylabel('v-VSYS'); a.set_title(f'{k} {mol} PV')
        a=ax[r,2*j+1]; a.imshow(M/np.maximum(M.max(0,keepdims=True),1e-9),origin='lower',aspect='auto',extent=[offs[0],offs[-1],v[0],v[-1]],cmap='viridis',vmin=-0.05,vmax=0.6); a.set_ylim(-700,600); a.set_title('column-normalised')
plt.tight_layout(); plt.savefig(OUT+'/s75a_pv.png',dpi=70)
# ---- region spectra & decomposition
def g1(v,a,c,s): return a*np.exp(-0.5*((v-c)/s)**2)
def g2(v,a1,c1,s1,a2,c2,s2): return g1(v,a1,c1,s1)+g1(v,a2,c2,s2)
def fit(v,y,e):
    m=(v>-800)&(v<650); v,y=v[m],y[m]; res={}
    a0=y.max(); c0=v[np.argmax(y)]
    try:
        p,cv=curve_fit(g1,v,y,p0=[a0,c0,60],sigma=np.full(len(v),e),absolute_sigma=True,maxfev=20000)
        chi=np.sum(((y-g1(v,*p))/e)**2); res['g1']=dict(p=p.tolist(),err=np.sqrt(np.diag(cv)).tolist(),chi2=chi,n=len(v))
    except Exception as ex: res['g1']=None
    best=None
    for c2 in (-300,-200,-100,0,100,200,300):
        for s2 in (60,120,200):
            try:
                p,cv=curve_fit(g2,v,y,p0=[a0,c0,50,0.2*a0,c2,s2],sigma=np.full(len(v),e),absolute_sigma=True,maxfev=20000,bounds=([0,-250,10,0,-700,40],[10*abs(a0)+1,250,200,10*abs(a0)+1,600,400]))
                chi=np.sum(((y-g2(v,*p))/e)**2)
                if best is None or chi<best[0]: best=(chi,p,np.sqrt(np.diag(cv)))
            except Exception: pass
    if best: res['g2']=dict(p=best[1].tolist(),err=best[2].tolist(),chi2=best[0],n=len(v))
    return res
def aperture_spec(R,pos,rad):
    return aspec(R['c'],R['ba'],apmask(R['n'],R['ps'],R['pxy'](*pos),rad))
def peak_spec(R,pos,search=0.6):
    """spectrum at the pixel of max moment-0 (|v|<350) within `search` arcsec of pos, native units Jy/beam"""
    px,py=R['pxy'](*pos); n=R['n']; Y,X=np.mgrid[:n,:n]; m=np.hypot(X-px,Y-py)*R['ps']<=search
    M0=R['c'][np.abs(R['v'])<350].sum(0); i=np.argmax(np.where(m,M0,-1e9)); y,x=np.unravel_index(i,M0.shape)
    return R['c'][:,y,x],(x,y)
REG={'S':0.0,'Slobe':-1.85,'Bridge_mid':+2.6,'Red3.3N':+3.3,'N':+5.22}
out={}
fig,ax=plt.subplots(len(REG),4,figsize=(20,3.2*len(REG)))
for ri,(rn,o) in enumerate(REG.items()):
    pos=offset_pos(S0,0.0 if rn!='N' else 0.156,o) if True else None
    if rn=='Slobe': pos=sout_pos((-0.2+0.15,-1.85)) if False else offset_pos(S0,0.15,-1.85)
    for ci,(k,kind) in enumerate((('Xb4d_1.05','ap'),('Xb4d_nat','peak'),('Xb4b_0.5','ap'),('Xb4b_0.5','peak'))):
        if k=='Xb4b_0.5' and kind=='peak' and False: pass
        for mol in ('HCN','HCO+'):
            R=SET[k][mol]; v=R['v']; lf=(v>-900)&(v<-450)
            if kind=='ap': sp=aperture_spec(R,pos,1.0); unit='Jy (r=1" aperture sum / beam area = beam-avg Jy/beam*N_beams)'
            else: sp,xy=peak_spec(R,pos); unit='Jy/beam'
            e=sp[lf].std(); f=fit(v,sp,e)
            key=f'{rn}|{k}|{kind}|{mol}'; out[key]=dict(rms_chan=float(e),fit=f)
            ax[ri,ci].plot(v,sp*1e3,lw=.8,label=mol)
            if f.get('g2'): ax[ri,ci].plot(v[(v>-800)&(v<650)],g2(v[(v>-800)&(v<650)],*f['g2']['p'])*1e3,'k--',lw=.6)
            s1=f['g1']; s2=f.get('g2')
            P(f'{key}: rms/chan {e*1e3:.2f} mJy ; 1G: A={s1["p"][0]*1e3:.1f} c={s1["p"][1]:.0f} FWHM={2.355*abs(s1["p"][2]):.0f} chi2={s1["chi2"]:.0f}/{s1["n"]-3}'+
              (f' ; 2G: A1={s2["p"][0]*1e3:.1f} c1={s2["p"][1]:.0f} FWHM1={2.355*s2["p"][2]:.0f} | A2={s2["p"][3]*1e3:.1f} c2={s2["p"][4]:.0f}+-{s2["err"][4]:.0f} FWHM2={2.355*s2["p"][5]:.0f} chi2={s2["chi2"]:.0f}/{s2["n"]-6} dchi2={s1["chi2"]-s2["chi2"]:.0f}' if s2 else ''))
        ax[ri,ci].set_title(f'{rn} {k} {kind}',fontsize=8); ax[ri,ci].set_xlim(-800,650)
plt.tight_layout(); plt.savefig(OUT+'/s75a_spec.png',dpi=60)
json.dump(out,open(OUT+'/s75a.json','w'),indent=1,default=float); open(OUT+'/s75a.txt','w').write('\n'.join(txt))
