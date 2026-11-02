# Q2 round 2, fix 6: geometry vs literature and a disk-rotation / warp null test.
# (a) positions/PA: lobe peak (s75b: dE=-0.10+-0.04, dN=-1.88+-0.11 from S) and Michiyama+2018 South-OF-Blue (Table 1) vs the CO outflow axis (perpendicular to the S nuclear disk PA~90 deg, Sakamoto+2014 Table 8; curves west).
# (b) narrow-component velocity field from Xb4b (0.25", HCO+ and HCN, 3.9 MHz cubes): per-pixel narrow peak velocity (parabolic peak of 3-channel-smoothed spectrum in |v|<200) and flux; plane-fit gradient at the lobe.
# (c) BEAM-SMEARED ROTATION NULL MODEL: synthetic cube where every 0.25" pixel carries ONLY a narrow Gaussian (sigma_n=25/40/60 km/s) at its own measured narrow-peak velocity with its own |v|<150 flux (no wing at all); smoothed to the 1.05" circular beam exactly like the data;
#     predicted wing flux in -405..-150 and +180..320 at the lobe / S / red feature / N compared with the observed wing from the same Xb4b cube smoothed to 1.05" (the pure velocity-gradient/beam-smearing content of a thin disk with the observed velocity structure).
# (d) E-W PV through the lobe (PA 95): narrow-peak velocity and wing centroid (-420..-100) vs E-W offset (Xb4d 1.05").
# Output out/q2r2/s77g.txt/json, s77g_null.png
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s77_lib import *
from scipy.ndimage import uniform_filter1d, map_coordinates
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q2r2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True)); res={}
def wsum(v,dv,sp,lo,hi):
    ov=np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv; return float((sp*ov).sum()*dv)
# (a)
M18=((10+27/60+51.17/3600)*15, -(43+54/60+21.43/3600)); M18S=((10+27/60+51.22/3600)*15,-(43+54/60+19.20/3600))
dE=(M18[0]-M18S[0])*np.cos(np.deg2rad(M18S[1]))*3600; dN=(M18[1]-M18S[1])*3600
pa=lambda e,n: np.rad2deg(np.arctan2(e,n))%360
P(f'(a) M18 South-OF-Blue relative to M18 S nucleus: dE={dE:+.2f}" dN={dN:+.2f}" -> PA {pa(dE,dN):.0f} deg, r={np.hypot(dE,dN):.2f}" = {np.hypot(dE,dN)*194:.0f} pc')
P(f'    our HCN blue wing peak (Xb4d@1.05"): dE=-0.10+-0.04 dN=-1.88+-0.11 -> PA {pa(-0.10,-1.88):.0f} deg (+-{np.rad2deg(np.hypot(0.04,0.11/1.88*0.1))+0:.0f} from positional errors), r=1.88" = 365 pc; HCO+: dE=-0.31+-0.12 dN=-1.79+-0.45 -> PA {pa(-0.31,-1.79):.0f}; Xb4b@1.05 HCN: dE=-0.37 dN=-2.62 -> PA {pa(-0.37,-2.62):.0f}')
P('    Sakamoto+2014: S nuclear disk PA~90 deg (edge-on, i~80), outflow axis along the rotation axis => PA ~0/180 deg, blue cone south, curving west with distance; 20 deg opening; CO high-velocity emission peaks at 1.8" (310 pc) each side. Pereira-Santaella+2024: outflow perpendicular to the disk, i~80 deg (near the plane of the sky).')
# (b),(c)
def prep(mol,half=7.0):
    c,v,ps,wc,x0,y0,b0=load_cube('Xb4b',SPW[mol],REST[mol],-1000,700,half=half)
    return c,v,ps,wc,x0,y0,b0
out={}
for mol in ('HCO+','HCN'):
    c,v,ps,wc,x0,y0,b0=prep(mol); dv=abs(np.median(np.diff(v))); pxy=pxy_fn(wc,x0,y0); n=c.shape[1]
    c=c-np.median(c[(v>-900)&(v<-450)],axis=0)
    lf=(v>-900)&(v<-450); rms=c[lf].std()
    kc=(v>-200)&(v<200); sm=uniform_filter1d(c[kc],3,axis=0); vv=v[kc]
    ip=np.argmax(sm,axis=0); ip=np.clip(ip,1,len(vv)-2)
    yy,xx=np.mgrid[:n,:n]; a=sm[ip-1,yy,xx]; b=sm[ip,yy,xx]; cc=sm[ip+1,yy,xx]; den=a-2*b+cc; off=np.where(np.abs(den)>0,0.5*(a-cc)/den,0)
    vpk=vv[ip]+np.clip(off,-1,1)*dv
    F=c[(v>=-150)&(v<=150)].sum(0)*dv; sgF=np.median(np.abs(F-np.median(F)))*1.4826
    ok=(b>4*rms/np.sqrt(3))&(F>5*sgF)
    P(f'### {mol} Xb4b: rms/chan {rms*1e3:.2f} mJy/b; pixels used for the null model {ok.sum()} of {n*n}')
    # plane fit around the lobe (r<1.5")
    for nm,pos in (('lobe',LOBE),('S',S0),('N',N0)):
        px,py=pxy(*pos); m=ok&(np.hypot(xx-px,yy-py)*ps<1.5)
        if m.sum()>10:
            A=np.c_[np.ones(m.sum()),(xx[m]-px)*ps,(yy[m]-py)*ps]; sol=np.linalg.lstsq(A,vpk[m],rcond=None)[0]; resid=vpk[m]-A@sol
            P(f'  narrow-peak velocity field within r<1.5" of {nm}: v0={sol[0]:+.0f} km/s, gradient dv/dx(pixel x, W->E negative sign convention: pixel x increases westward)={sol[1]:+.0f}, dv/dy(N)={sol[2]:+.0f} km/s/arcsec, residual rms {resid.std():.0f} km/s, v range {vpk[m].min():.0f}..{vpk[m].max():.0f}; => max |dv| across a 1.05" beam ~ {np.hypot(sol[1],sol[2])*1.05:.0f} km/s; rotational smearing sigma ~ {np.hypot(sol[1],sol[2])*1.05/2.355:.0f} km/s')
            out[f'{mol}|grad|{nm}']=(float(sol[0]),float(sol[1]),float(sol[2]),float(resid.std()),int(m.sum()))
    fig,ax=plt.subplots(1,2,figsize=(12,5)) if mol=='HCO+' else (None,[None,None])
    if mol=='HCO+':
        im=ax[0].imshow(np.where(ok,vpk,np.nan),origin='lower',cmap='RdBu_r',vmin=-150,vmax=150); ax[0].set_title('Xb4b HCO+ narrow-peak velocity'); ax[1].imshow(F,origin='lower',vmax=np.percentile(F,99.5)); plt.savefig(OUT+'/s77g_vfield.png',dpi=60); plt.close()
    for sg in (25.,40.,60.):
        syn=np.zeros_like(c); amp=F/(np.sqrt(2*np.pi)*sg)*ok
        syn=amp[None]*np.exp(-0.5*((v[:,None,None]-vpk[None])/sg)**2)
        # exact same smoothing as the data
        cs=conv(syn.astype(np.float32),ps,b0,(1.05,1.05,0.0)); co=conv(c.astype(np.float32),ps,b0,(1.05,1.05,0.0)); ba=1.1331*1.05**2/ps**2
        for nm,pos,rad in (('lobe',LOBE,1.0),('S',S0,1.0),('red',RED,1.0),('N',N0,1.0)):
            m=apmask(n,ps,pxy(*pos),rad); so=aspec(co,ba,m); ss_=aspec(cs,ba,m)
            r_={f'{w[0]}..{w[1]}':(wsum(v,dv,so,*w),wsum(v,dv,ss_,*w)) for w in ((-405,-150),(150,405),(180,320),(-150,150))}
            out[f'{mol}|null|{sg:.0f}|{nm}']=r_
            P(f'  null sigma_n={sg:.0f}: {nm:4s} r=1": observed wing(-405..-150) {r_['-405..-150'][0]:+.3f} vs predicted {r_['-405..-150'][1]:+.3f}; (+150..405) obs {r_['150..405'][0]:+.3f} pred {r_['150..405'][1]:+.3f}; core(-150..150) obs {r_['-150..150'][0]:.2f} pred {r_['-150..150'][1]:.2f} Jy km/s')
# (d) E-W PV at the lobe and at S in Xb4d 1.05"
for mol in ('HCN','HCO+'):
    R=load_any('std_Xb4d',mol,1.05,12.0); v=R['v']; c=R['c']-np.median(R['c'][(v>-900)&(v<-450)],axis=0); n=R['n']; ps=R['ps']; dv=R['dv']
    for pn,pos in (('lobe',LOBE),('S',S0)):
        row=[]
        for o in np.arange(-3,3.01,0.5):   # E positive
            p2=offset_pos(pos,o,0.0); m=apmask(n,ps,R['pxy'](*p2),0.5); sp=aspec(c,R['ba'],m)*1.0
            k=(np.abs(v)<200); vp=v[k][np.argmax(uniform_filter1d(sp[k],3))]; kw=(v>=-420)&(v<=-100); w=np.clip(sp[kw],0,None); vw=(w*v[kw]).sum()/w.sum() if w.sum()>0 else np.nan
            row.append((o,vp,vw,wsum(v,dv,sp,-405,-150)))
        out[f'{mol}|EW|{pn}']=row
        P(f'(d) {mol} E-W cut through {pn} (r=0.5" apertures, E offset: narrow-peak v / wing centroid(-420..-100) / wing flux(-405..-150)): '+'; '.join(f'{o:+.1f}":{vp:+.0f}/{vw:+.0f}/{fl:+.2f}' for o,vp,vw,fl in row))
P('Analytic bound: to put gas at |v_los| = 180 km/s inside a 1.05" beam by rotation alone the velocity must change by >=180 km/s across <= half the beam, i.e. a gradient >= 180/(0.5*1.05) = 340 km/s/arcsec; the measured narrow-peak gradient at the lobe is 28-45 km/s/arcsec (HCO+, HCN), i.e. 8-12 times too small.')
json.dump(out,open(OUT+'/s77g.json','w'),indent=1,default=float); open(OUT+'/s77g.txt','w').write('\n'.join(txt))
