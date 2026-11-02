# Q17 round 2, step b: independent-velocity-bin evidence and resolution/product dependence of the ridge (reviewer fix 2).
# (i) 40 km/s bins (~9 channels of 4.4 km/s, ~independent given Hanning) from +120..+400 (ridge) and mirrored -400..-120 (lobe) at 1.05": ridge/lobe flux per bin, matched-core control scatter,
#     number of bins >3 sigma_ctrl; maps in the +-(160-360) window in 4 bins -> figure.
# (ii) products: std_Xb4d (0.94x0.74"), tap_Xb4d (0.96x0.84" hand-tapered), std_Xb4b (0.25") convolved to 1.05", HCN & HCO+; ridge/lobe 180-320 flux. Xb4f (5.7x3.5") cannot resolve ridge from N/S: wing flux in a 4" aperture only.
# usage: python3 -I s106b_q17_channels_res.py
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.signal import fftconvolve
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q17r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
RIDGE=offset_pos(S0,0.15,3.3); LOBE_=offset_pos(S0,0.15,-1.85)
BINS=[(120,160),(160,200),(200,240),(240,280),(280,320),(320,360),(360,400)]
fig,ax=plt.subplots(4,5,figsize=(18,14))
for jm,mol in enumerate(('HCN','HCO+')):
    R=load_any('std_Xb4d',mol,1.05,24.0); v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
    h=(np.hypot(*np.mgrid[-int(1.2/ps):int(1.2/ps)+1,-int(1.2/ps):int(1.2/ps)+1])*ps<=1.0).astype(float)/R['ba']
    def Wm(a,b):
        k=np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
    Mc=fftconvolve(Wm(-150,150),h,mode='same'); Y,X=np.mgrid[:n,:n]
    def dist(pos): x,y=R['pxy'](*pos); return np.hypot(X-x,Y-y)*ps
    dR,dL,dN,dS=dist(RIDGE),dist(LOBE_),dist(N0),dist(S0)
    mx,my=R['pxy'](*offset_pos(S0,0.1,2.55)); sel=np.hypot(X-mx,Y-my)*ps<10; step=int(round(1/ps)); grid=np.zeros((n,n),bool); grid[::step,::step]=True
    ir=tuple(int(round(t)) for t in R['pxy'](*RIDGE)[::-1]); il=tuple(int(round(t)) for t in R['pxy'](*LOBE_)[::-1])
    okc=sel&grid&~((dR<1.5)|(dL<1.5)|(dN<1.5)|(dS<1.5))
    mcr=okc&(Mc>0.5*Mc[ir])&(Mc<2*Mc[ir]); mcl=okc&(Mc>0.5*Mc[il])&(Mc<2*Mc[il])
    fm=field_pixels(R,rmin=3.5,rmax=22.0); P(f'### {mol}: matched-core control grid points: ridge-like {mcr.sum()}, lobe-like {mcl.sum()}; all-sky controls {okc.sum()}')
    for tag,sg_,pos_,ix,mc_ in (('ridge',1,RIDGE,ir,mcr),('lobe',-1,LOBE_,il,mcl)):
        out=[]
        for a,b in BINS:
            M=fftconvolve(Wm(a,b) if sg_>0 else Wm(-b,-a),h,mode='same'); f=M[ix]; ctrl=M[okc]; cm=M[mc_]
            sig=np.std(ctrl); sigm=np.std(cm) if mc_.sum()>=5 else sig
            out.append((a,b,f,sig,float(cm.mean()) if mc_.sum() else np.nan,sigm))
        P(f'  {mol} {tag} (|v| bins, mean-over-r=1" aperture Jy km/s; sigma=sd of all-control apertures [matched-core mean/sd]): '+'; '.join(f'{a}-{b}: {f:.3f} ({f/s:.1f}s) [{m:+.3f}/{sm:.3f}]' for a,b,f,s,m,sm in out))
        res[f'{mol}|{tag}|bins']=[dict(a=a,b=b,f=float(f),sig=float(s),ctrl_mean=m,ctrl_sd=float(sm)) for a,b,f,s,m,sm in out]
        nb=sum(1 for a,b,f,s,m,sm in out if f>3*s)
        P(f'     bins above 3 sigma(all-control): {nb}/{len(out)}; consecutive bins 160-320 all >3 sigma: {all(o[2]>3*o[3] for o in out[1:5])}')
    # maps: rows: red HCN/HCO+ ; blue
    rowr=jm*2; 
    for j,(a,b) in enumerate([(160,200),(200,240),(240,280),(280,320),(320,360)]):
        for row,sg_ in ((rowr,1),(rowr+1,-1)):
            M=Wm(a,b) if sg_>0 else Wm(-b,-a); M=fftconvolve(M,h,mode='same')*R['ba']*0+fftconvolve(M,h,mode='same')*0+fftconvolve(M,h,mode='same')
            a_=ax[row,j]; x0,y0=R['pxy'](*S0); ext=7/ps
            a_.imshow(M[int(y0-ext*0.6):int(y0+ext*1.2),int(x0-ext):int(x0+ext)],origin='lower',cmap='magma',vmin=-0.03,vmax=0.15 if mol=='HCN' else 0.12)
            sl=(slice(int(y0-ext*0.6),int(y0+ext*1.2)),slice(int(x0-ext),int(x0+ext)))
            for nm,pos,col in (('S',S0,'c'),('N',N0,'w'),('R',RIDGE,'lime'),('L',LOBE_,'lime')):
                x,y=R['pxy'](*pos); a_.plot(x-sl[1].start,y-sl[0].start,'+' if nm in 'SN' else 'o',color=col,mfc='none',ms=9)
                a_.text(x-sl[1].start+2,y-sl[0].start+2,nm,color=col,fontsize=8)
            a_.set_title(f'{mol} {"+" if sg_>0 else "-"}{a}..{b} km/s (r=1" mean)',fontsize=8); a_.set_xticks([]); a_.set_yticks([])
plt.tight_layout(); plt.savefig(OUT+'/s106b_channel_maps.png',dpi=70); plt.close()
# ---- (ii) products / resolution
for mol in ('HCN','HCO+'):
    rows=[]
    for prod in ('std_Xb4d','tap_Xb4d','std_Xb4b'):
        try: R=load_any(prod,mol,1.05,12.0)
        except Exception as e: P(f'{prod} {mol}: not loadable ({e})'); continue
        v=R['v']; dv=R['dv']; ps=R['ps']; n=R['n']
        if prod=='tap_Xb4d': c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0) if ((v>-900)&(v<-450)).sum()>3 else R['c']
        else: c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0)
        def Wm(a,b):
            k=np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
        fm=field_pixels(R,rmin=3.5,rmax=11.0)
        h=(np.hypot(*np.mgrid[-int(1.2/ps):int(1.2/ps)+1,-int(1.2/ps):int(1.2/ps)+1])*ps<=1.0).astype(float)/R['ba']
        mr=fftconvolve(Wm(180,320),h,mode='same'); mb=fftconvolve(Wm(-320,-180),h,mode='same')
        ir=tuple(int(round(t)) for t in R['pxy'](*RIDGE)[::-1]); il=tuple(int(round(t)) for t in R['pxy'](*LOBE_)[::-1])
        P(f'  {mol} {prod}: dv={dv:.1f} km/s, nchan in window {((v>180)&(v<320)).sum()}; ridge {mr[ir]:.3f}+-{np.std(mr[fm]):.3f}, lobe {mb[il]:.3f}+-{np.std(mb[fm]):.3f} Jy km/s (r=1", 1.05" beam), ridge/lobe {mr[ir]/mb[il]:.2f}')
        res[f'{mol}|{prod}']=dict(ridge=float(mr[ir]),lobe=float(mb[il]),sr=float(np.std(mr[fm])),sl=float(np.std(mb[fm])))
# Xb4f (5.7x3.5") : ridge, N and S blend: red wing in r=4" aperture centred on the ridge, vs blue wing in the mirrored aperture at the lobe (r=4": contains S)
for mol,spw in (('HCN','spw18'),('HCO+','spw16')):
    c,v,ps,wc,x0,y0,b0=load_cube('Xb4f',spw,REST[mol],-1000,700,half=24)
    dv=abs(np.median(np.diff(v))); c=c-np.median(c[(v>-900)&(v<-450)],axis=0); pxy=pxy_fn(wc,x0,y0); n=c.shape[1]
    ba=1.1331*b0[0]*b0[1]/ps**2
    def Wm(a,b):
        k=np.clip(np.minimum(v+dv/2,b)-np.maximum(v-dv/2,a),0,dv)/dv; return np.tensordot(k,c,axes=(0,0))*dv
    r_=Wm(180,320); b_=Wm(-320,-180)
    for nm,pos in (('ridge',RIDGE),('N',N0),('S',S0),('lobe',LOBE_)):
        m=apmask(n,ps,pxy(*pos),2.5); P(f'  {mol} Xb4f (beam {b0[0]:.1f}x{b0[1]:.1f}"): r=2.5" aperture at {nm}: red180-320 {r_[m].sum()/ba:.3f}, blue-320..-180 {b_[m].sum()/ba:.3f} Jy km/s')
        res[f'{mol}|Xb4f|{nm}']=dict(red=float(r_[m].sum()/ba),blue=float(b_[m].sum()/ba))
json.dump(res,open(OUT+'/s106b.json','w'),indent=1); open(OUT+'/s106b.txt','w').write('\n'.join(txt))
