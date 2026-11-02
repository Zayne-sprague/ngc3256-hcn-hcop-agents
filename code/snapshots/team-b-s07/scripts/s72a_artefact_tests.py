# Q11 round 2, fix 1: artefact tests of the Xb4d-only HCN skirt at Soutflow.
# (a) maps of HCN, HCO+ and E=HCN-k*HCO+ integrated over -405..-205 (blue) and +205..+405 (red null) km/s about VSYS, whole Xb4d field, 1.05" beam
# (b) same quantity in r=1" apertures on a ring (8 azimuths) 1.85" from S and from N, + Soutflow itself;  noise = 60 random line-free apertures
# (c) Xb4d clean-mask cube in the same channels: mask fraction in the apertures, window flux inside vs outside the mask (native resolution)
# (d) channel-by-channel mean (bowl) and rms around S, HCN (spw27) vs HCO+ (spw25), raw and baseline-subtracted
# run: python3 -I scripts/s72a_artefact_tests.py   -> out/q11r2/s72a.txt, s72a.json, s72a_maps.png, s72a_mask.png, s72a_bowl.png
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q11r2'; rng=np.random.default_rng(72)
txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
BL=(-405,-205); RD=(205,405); KS,KO=0.37,0.534     # k: S-core ratio (Jy), Soutflow core ratio (Jy) from s70
def base(c,v): return c-np.median(c[(v>-900)&(v<-450)],axis=0)
RAW={};CV={};MK={}
for mol,spw in (('HCN','spw27'),('HCO+','spw25')):
    c,v,ps,wc,x0,y0,b0=load_cube('Xb4d',spw,F0[mol][3]); m,vm,_,_,_,_,_=load_cube('Xb4d',spw,F0[mol][3],kind='mask')
    RAW[mol]=dict(c=c,v=v,ps=ps,b0=b0,pxy=pxy_fn(wc,x0,y0),mask=m>0.5)
    cc=conv(c,ps,b0,(1.05,1.05,0.0)); cc=base(cc,v); CV[mol]=dict(c=cc,v=v,dv=abs(np.median(np.diff(v))),ps=ps,ba=1.1331*1.05**2/ps**2,pxy=RAW[mol]['pxy'])
H=CV['HCN']; O=CV['HCO+']; v=H['v']; dv=abs(np.median(np.diff(v))); n=H['c'].shape[1]; ps=H['ps']; pxy=H['pxy']
P(f'# Xb4d spw27 (HCN 3-2) / spw25 (HCO+ 3-2), pixel {ps:.3f}", cube {H["c"].shape}, dv {dv:.2f} km/s, native beam HCN {RAW["HCN"]["b0"][0]:.2f}x{RAW["HCN"]["b0"][1]:.2f} HCO+ {RAW["HCO+"]["b0"][0]:.2f}x{RAW["HCO+"]["b0"][1]:.2f}; analysis beam 1.05" circular; baseline per-pixel median -900..-450')
P(f'# windows (about VSYS=2780 km/s, radio): blue {BL}, red mirror {RD}; E=HCN-k*HCO+ with k={KS} (S core) and {KO} (Soutflow core)')
Spos=REGIONS['S'][:2]; Npos=REGIONS['N'][:2]; SO=sout_pos()
# ---------------- (a) maps
def mom(R,a,b): k=(R['v']>=a)&(R['v']<=b); return R['c'][k].sum(0)*R['dv']      # Jy/beam km/s
M={}
for nm,(a,b) in (('blue',BL),('red',RD)):
    M['HCN_'+nm]=mom(H,a,b); M['HCO_'+nm]=mom(O,a,b)
    for kk,kn in ((KS,'kS'),(KO,'kO')): M[f'E_{nm}_{kn}']=M['HCN_'+nm]-kk*M['HCO_'+nm]
m0=H['c'][np.abs(v)<350].sum(0)*dv; sg0=np.median(np.abs(m0-np.median(m0)))*1.4826
bright=binary_dilation(m0>4*sg0,iterations=int(1.5/ps)); Y,X=np.mgrid[:n,:n]
cx,cy=pxy(*CEN); rcen=np.hypot(X-cx,Y-cy)*ps; field=(~bright)&(rcen>4)&(rcen<12)
SIG={k:float(np.std(M[k][field])) for k in M}
P('## (a) field maps: pixel rms in line-free field (4<r<12" from N-S midpoint, bright masked), Jy/beam km/s: '+', '.join(f'{k} {s:.4f}' for k,s in SIG.items()))
sx,sy=pxy(*Spos); nx,ny=pxy(*Npos); ox,oy=pxy(*SO)
near=np.hypot(X-sx,Y-sy)*ps<4.0
RES={'a':{}}
for k in ('E_blue_kS','E_blue_kO','E_red_kS','E_red_kO','HCN_blue','HCN_red','HCO_blue','HCO_red'):
    mp=np.where(near,M[k],-1e9); iy,ix=np.unravel_index(np.argmax(mp),mp.shape)
    # contiguous area above half-max around peak
    from scipy.ndimage import label
    lab,_=label(M[k]>0.5*M[k][iy,ix]); area=(lab==lab[iy,ix]).sum()*ps**2
    # field-wide pixel statistic: number of independent-beam positions in r<12" (not nuclei) with value >= this
    nfield=((rcen<12)&(~bright)).sum(); nex=(M[k][(rcen<12)&(~bright)]>=M[k][iy,ix]).sum()
    dxE=(ix-sx)*ps*-1; dyN=(iy-sy)*ps    # +x pixel = west for standard WCS (RA decreases): E offset = -dx
    RES['a'][k]=dict(peak=float(M[k][iy,ix]),sig=float(M[k][iy,ix]/SIG[k]),dE=float(dxE),dN=float(dyN),area_halfmax_as2=float(area),n_exceed_field=int(nex),n_field=int(nfield))
    P(f'   {k:10s} peak within 4" of S: {M[k][iy,ix]:+.3f} Jy/b km/s = {M[k][iy,ix]/SIG[k]:+.1f} sigma at offset (E {dxE:+.2f}", N {dyN:+.2f}") from S, half-max area {area:.1f} arcsec2 (beam 0.87); field pixels (non-bright, r<12") >= this value: {nex}/{nfield}')
# values at Soutflow
for k in ('HCN_blue','HCO_blue','E_blue_kO','HCN_red','HCO_red','E_red_kO'):
    P(f'   at Soutflow pixel {k}: {M[k][int(round(oy)),int(round(ox))]:+.3f} ({M[k][int(round(oy)),int(round(ox))]/SIG[k]:+.1f} sigma)')
np.savez(OUT+'/s72a_maps.npz',**M,ps=ps,s=(sx,sy),n=(nx,ny),so=(ox,oy))
fig,ax=plt.subplots(2,3,figsize=(15,9)); ext=lambda:None
tit=[('HCN_blue','HCN 3-2 -405..-205'),('HCO_blue','HCO+ 3-2 -405..-205'),('E_blue_kO','E=HCN-0.534 HCO+ (blue)'),('HCN_red','HCN +205..+405'),('HCO_red','HCO+ +205..+405'),('E_red_kO','E red mirror')]
for a_,(k,t) in zip(ax.ravel(),tit):
    lim=4*SIG[k]*2.5; im=a_.imshow(M[k],origin='lower',cmap='RdBu_r',vmin=-0.3,vmax=0.3,extent=[0,n*ps,0,n*ps])
    a_.contour((np.arange(n)+.5)*ps,(np.arange(n)+.5)*ps,m0,levels=[0.1*m0.max(),0.3*m0.max(),0.6*m0.max()],colors='k',linewidths=.6)
    for (xx,yy,c_,l_) in ((sx,sy,'g','S'),(nx,ny,'m','N'),(ox,oy,'y','Sout')): a_.plot(xx*ps,yy*ps,'o',color=c_,ms=5); 
    th=np.linspace(0,2*np.pi,100); a_.plot((sx+1.85/ps*np.sin(th))*ps,(sy+1.85/ps*np.cos(th))*ps,'g--',lw=.7)
    a_.set_xlim((sx-8/ps)*ps,(sx+8/ps)*ps); a_.set_ylim((sy-8/ps)*ps,(sy+8/ps)*ps); a_.set_title(t+f' (pix rms {SIG[k]*1e3:.0f} mJy/b km/s)',fontsize=9); plt.colorbar(im,ax=a_,fraction=.046)
plt.tight_layout(); plt.savefig(OUT+'/s72a_maps.png',dpi=90); plt.close()
# ---------------- (b) ring apertures
RP=random_pos(H['c'],v,ps,pxy,60,1.0,5.0,9.5,rng,excl=[pxy(*Spos),pxy(*Npos),pxy(*SO)])
def ap_sums(xy,rad=1.0):
    out={}
    for mol,R in (('HCN',H),('HCO',O)):
        sp=aspec(R['c'],R['ba'],apmask(n,ps,xy,rad)); out[mol]={'blue':wsum(R['v'],sp,*BL,R['dv']),'red':wsum(R['v'],sp,*RD,R['dv']),'core':wsum(R['v'],sp,-205,95,R['dv']),'sp':sp}
    return out
nz={'HCN':np.array([[ap_sums(xy)['HCN'][w] for w in ('blue','red')] for xy in RP]),'HCO':np.array([[ap_sums(xy)['HCO'][w] for w in ('blue','red')] for xy in RP])}
sgm={m:nz[m].std(0) for m in nz}; ME={'E_blue_kS':None}
P(f'\n## (b) apertures r=1" (Jy km/s), noise from {len(RP)} random apertures: sigma HCN blue {sgm["HCN"][0]:.3f} red {sgm["HCN"][1]:.3f}; HCO+ blue {sgm["HCO"][0]:.3f} red {sgm["HCO"][1]:.3f}; random-aperture means (HCN blue {nz["HCN"][:,0].mean():+.3f}, red {nz["HCN"][:,1].mean():+.3f}) = baseline bias check')
def Esig(h,c,k,w): return h-k*c, np.sqrt(sgm['HCN'][w]**2+(k*sgm['HCO'][w])**2)
RB={}
def report(label,xy):
    s=ap_sums(xy); row={}
    for w,wn in ((0,'blue'),(1,'red')):
        h,c=s['HCN'][wn],s['HCO'][wn]
        e1,s1=Esig(h,c,KS,w); e2,s2=Esig(h,c,KO,w); row[wn]=dict(HCN=h,HCO=c,E_kS=e1,E_kO=e2,sig_kS=s1,sig_kO=s2,zE_kO=e2/s2,zHCN=h/sgm['HCN'][w],zHCO=c/sgm['HCO'][w])
    row['core']=(s['HCN']['core'],s['HCO']['core']); RB[label]=row
    b,r=row['blue'],row['red']
    P(f'  {label:14s} core(-205..95) HCN {row["core"][0]:.2f} HCO+ {row["core"][1]:.2f} | BLUE HCN {b["HCN"]:+.3f} ({b["zHCN"]:+.1f}s) HCO+ {b["HCO"]:+.3f} ({b["zHCO"]:+.1f}s) E(k.534) {b["E_kO"]:+.3f}+-{b["sig_kO"]:.3f} ({b["zE_kO"]:+.1f}s) | RED HCN {r["HCN"]:+.3f} ({r["zHCN"]:+.1f}s) HCO+ {r["HCO"]:+.3f} ({r["zHCO"]:+.1f}s) E {r["E_kO"]:+.3f} ({r["zE_kO"]:+.1f}s)')
report('Soutflow',pxy(*SO))
for cen,cn,cpos in (('S','S',Spos),('N','N',Npos)):
    P(f'  -- ring 1.85" from {cn} (azimuth deg E of N; Soutflow offset is az ~186 deg from S)')
    for az in range(0,360,45):
        pos_=offset_pos(cpos,1.85*np.sin(np.deg2rad(az)),1.85*np.cos(np.deg2rad(az))); report(f'{cn}_az{az:03d}',pxy(*pos_))
    ring=[RB[f'{cn}_az{a:03d}'] for a in range(0,360,45)]
    zb=np.array([r['blue']['zE_kO'] for r in ring]); zr=np.array([r['red']['zE_kO'] for r in ring]); hb=np.array([r['blue']['zHCN'] for r in ring])
    P(f'     ring {cn}: E_blue z {np.round(zb,1)} | E_red z {np.round(zr,1)} | HCN blue z {np.round(hb,1)} ; ring mean E_blue(k.534) {np.mean([r["blue"]["E_kO"] for r in ring]):+.3f}, azimuthal rms {np.std([r["blue"]["E_kO"] for r in ring]):.3f} (stat sigma {ring[0]["blue"]["sig_kO"]:.3f})')
RES['b']={k:{kk:vv for kk,vv in r.items()} for k,r in RB.items()}
# ---------------- (c) mask
P('\n## (c) Xb4d clean masks (native pixels; mask value>0.5): fraction of r=1" aperture pixels inside the mask, per 50 km/s of channels')
mk={}; 
for mol in ('HCN','HCO+'):
    R=RAW[mol]; nn=R['c'].shape[1]; mk[mol]=R
posd={'Soutflow':SO,'S':Spos,'N':Npos}
for az in (0,90,180,270): posd[f'S_az{az}']=offset_pos(Spos,1.85*np.sin(np.deg2rad(az)),1.85*np.cos(np.deg2rad(az)))
posd['S_az135']=offset_pos(Spos,1.85*np.sin(np.deg2rad(135)),1.85*np.cos(np.deg2rad(135)))
posd['S_az225']=offset_pos(Spos,1.85*np.sin(np.deg2rad(225)),1.85*np.cos(np.deg2rad(225)))
RC={}
for mol in ('HCN','HCO+'):
    R=RAW[mol]; vv=R['v']; nn=R['c'].shape[1]; P(f'  {mol}: mask nonzero channels span {vv[R["mask"].any((1,2))].min():.0f}..{vv[R["mask"].any((1,2))].max():.0f} km/s; total masked pixels in channels -405..-205: {R["mask"][(vv>=-405)&(vv<=-205)].sum()} of {R["mask"][0].size*((vv>=-405)&(vv<=-205)).sum()}')
    bins=[(a,a+50) for a in range(-550,200,50)]
    P('     bin(km/s):        '+' '.join(f'{(a+b)//2:+5d}' for a,b in bins))
    for pn,pp in posd.items():
        xy=R['pxy'](*pp); am=apmask(nn,R['ps'],xy,1.0); fr=[R['mask'][(vv>=a)&(vv<b)][:,am].mean() for a,b in bins]
        RC.setdefault(mol,{})[pn]=fr
        P(f'     {pn:9s} mask frac '+' '.join(f'{f:5.2f}' for f in fr))
    # window flux inside vs outside mask at native resolution, aperture r=1"
    for pn in ('Soutflow','S_az0','S_az90','S_az270'):
        xy=R['pxy'](*posd[pn]); am=apmask(nn,R['ps'],xy,1.0); k=(vv>=-405)&(vv<=-205)
        cub=R['c']-np.median(R['c'][(vv>-900)&(vv<-450)],axis=0); ba=1.1331*R['b0'][0]*R['b0'][1]/R['ps']**2; dvm=abs(np.median(np.diff(vv)))
        ins=(cub[k]*(R['mask'][k]&am)).sum()*dvm/ba; out_=(cub[k]*((~R['mask'][k])&am)).sum()*dvm/ba
        P(f'     {mol} {pn}: native window (-405..-205) flux in aperture: inside mask {ins:+.3f}, outside mask {out_:+.3f} Jy km/s (total {ins+out_:+.3f})')
        RC[mol][pn+'_inout']=(float(ins),float(out_))
RES['c']=RC
fig,ax=plt.subplots(2,5,figsize=(17,7))
for i,mol in enumerate(('HCN','HCO+')):
    R=RAW[mol]; vv=R['v']
    for j,vc in enumerate((-405,-355,-305,-255,-205)):
        ch=np.argmin(abs(vv-vc)); a_=ax[i,j]; a_.imshow(R['c'][ch],origin='lower',cmap='gray_r',vmin=-0.02,vmax=0.05); a_.contour(R['mask'][ch],levels=[0.5],colors='r',linewidths=.8)
        for pn,c_ in (('Soutflow','y'),('S','g'),('N','m')):
            xy=R['pxy'](*posd[pn]); a_.add_patch(plt.Circle(xy,1.0/R['ps'],fill=False,color=c_))
        sxy=R['pxy'](*Spos); a_.set_xlim(sxy[0]-60,sxy[0]+60); a_.set_ylim(sxy[1]-60,sxy[1]+60); a_.set_title(f'{mol} {vv[ch]:+.0f} km/s',fontsize=9)
plt.tight_layout(); plt.savefig(OUT+'/s72a_mask.png',dpi=80); plt.close()
# ---------------- (d) channel statistics around S
P('\n## (d) per-channel statistics around S (native-resolution raw pbcor cubes unless stated); annulus 1.5-4.0" from S excluding r<2" of N and Bridge and bright-mask pixels; off-source: r>8" from midpoint, 60-pixel blocks')
RD_={}
for mol in ('HCN','HCO+'):
    R=RAW[mol]; vv=R['v']; nn=R['c'].shape[1]; ps_=R['ps']; sxy=R['pxy'](*Spos); Y_,X_=np.mgrid[:nn,:nn]
    ann=(np.hypot(X_-sxy[0],Y_-sxy[1])*ps_>=1.5)&(np.hypot(X_-sxy[0],Y_-sxy[1])*ps_<=4.0)
    for q in ('N','Bridge'):
        qx=R['pxy'](*REGIONS[q][:2]); ann&=np.hypot(X_-qx[0],Y_-qx[1])*ps_>2.0
    cxy=R['pxy'](*CEN); off=np.hypot(X_-cxy[0],Y_-cxy[1])*ps_>8.0
    ba=1.1331*R['b0'][0]*R['b0'][1]/ps_**2
    mean_ann=R['c'][:,ann].mean(1); rms_off=R['c'][:,off].std(1); mean_off=R['c'][:,off].mean(1)
    ll=(vv>-900)&(vv<-500)
    P(f'  {mol}: off-source per-channel rms (median over channels -900..-500) {np.median(rms_off[ll])*1e3:.2f} mJy/beam; annulus mean in same channels {mean_ann[ll].mean()*1e3:+.2f} mJy/beam')
    for (a,b) in ((-405,-305),(-305,-205),(-205,-105),(-105,0),(150,250),(250,350),(350,450)):
        k=(vv>=a)&(vv<b)
        P(f'     v {a:+4d}..{b:+4d}: annulus mean {mean_ann[k].mean()*1e3:+.2f} +- {np.std(mean_ann[ll])/np.sqrt(1)*1e3:.2f}(chan scatter of annulus mean at line-free v) mJy/b ; offsource mean {mean_off[k].mean()*1e3:+.2f}, rms {rms_off[k].mean()*1e3:.2f} mJy/b; ratio rms/(line-free rms) {rms_off[k].mean()/np.median(rms_off[ll]):.2f}')
    RD_[mol]=dict(v=vv.tolist(),mean_ann=mean_ann.tolist(),rms_off=rms_off.tolist(),mean_off=mean_off.tolist())
RES['d']=RD_
fig,ax=plt.subplots(2,1,figsize=(9,7),sharex=True)
for mol,c_ in (('HCN','C0'),('HCO+','C1')):
    d=RD_[mol]; ax[0].plot(d['v'],np.array(d['mean_ann'])*1e3,c_,label=mol+' annulus mean'); ax[1].plot(d['v'],np.array(d['rms_off'])*1e3,c_,label=mol+' off-source rms')
ax[0].axvspan(-405,-205,color='y',alpha=.2); ax[1].axvspan(-405,-205,color='y',alpha=.2); ax[0].set_ylabel('mJy/beam'); ax[1].set_ylabel('mJy/beam'); ax[0].legend(); ax[1].legend(); ax[1].set_xlabel('v-VSYS km/s')
plt.tight_layout(); plt.savefig(OUT+'/s72a_bowl.png',dpi=80); plt.close()
json.dump(RES,open(OUT+'/s72a.json','w'),indent=1,default=float); open(OUT+'/s72a.txt','w').write('\n'.join(txt))
