"""Q29: Band 3 (2016.1.00965.S) 1-0 line integrals at N, S, Bridge, S blue lobe, ridge; matched circular beams (1.0", 1.29", 1.7"), plus Xb4d 3-2 at the same beams.
Method (as s170): pbcor cube -> RJ T_b (native beam) -> unit-area-kernel convolution -> per-pixel median baseline from line-free channels -> moment 0 in an aperture
(radius = beam FWHM) centred on the HCO+ 1-0 moment peak within 0.4" of nominal. Windows: N,S,Bridge -350..+350 km/s about 2780 (and +-150 as Harada); Slobe -320..-180, Ridge +180..+320.
Weak lines (H13CN, H13CO+, SiO, HN13C, CCH, HNC, HC3N): box +-150 km/s moment and a joint template fit (template = HCO+ 1-0 aperture spectrum of same region, free amplitude per line).
usage: python3 -I s182_q29_measure.py  -> out/q29/s182.json"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from regions import REGIONS
OUT='/workspace/out/q29/'
def off(pos,dE,dN): return (pos[0]+dE/3600/np.cos(np.deg2rad(pos[1])),pos[1]+dN/3600)
S0=REGIONS['S'][:2]
REG={k:REGIONS[k][:2] for k in ('N','S','Bridge')}; REG['Slobe']=off(S0,0.15,-1.85); REG['Ridge']=off(S0,0.15,3.3)
SPECIAL={'Slobe':(-320,-180),'Ridge':(180,320)}
THETAS=[1.0,1.29,1.7]
# rest freqs GHz
F0=dict(H13CN=86.3402,H13COp=86.7543,SiO=86.8470,HN13C=87.0908,CCHa=87.3169,CCHb=87.4020,HCN=88.6318,HCOp=89.1885,HNC=90.6636,HC3N=90.9790,HCN32=265.8864,HCOp32=267.5576)
CUBE_LINES={'X94a_spw25':['H13CN','H13COp','SiO','HN13C','CCHa','CCHb'],'X94a_spw27':['HCN','HCOp'],'X94e_spw25':['HNC','HC3N']}
def vel(nu,f0): return C*(1-nu/f0)-V0
def kern_w(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
def prep(key,theta,lines):
    cut,nu,pix,bm,pxy=load(key,half=16)
    fn=np.mean(nu)
    Tb=to_tb_conv(cut,nu,pix,bm,theta,fn)*1.0   # K
    # line-free mask: channels farther than 600 km/s from every known line in this cube (plus a 3rd line list: vib HCN, HOC+, HNCO)
    allf=[F0[l] for l in lines]+[89.487,87.925,88.8,87.402,87.317]
    free=np.ones(len(nu),bool)
    for f in allf: free&=np.abs(vel(nu,f))>600
    cont=np.median(Tb[free],axis=0); Tb=Tb-cont
    return Tb,nu,pix,bm,pxy,free,cont
def band3(theta):
    R={}; D={}
    for key,lines in CUBE_LINES.items():
        D[key]=prep(key,theta,lines)
    Tb,nu,pix,bm,pxy,free,cont=D['X94a_spw27']; vO=vel(nu,F0['HCOp']); vH=vel(nu,F0['HCN']); dv=abs(vO[1]-vO[0])
    ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    for rn,pos in REG.items():
        px,py=pxy(pos)
        if rn in SPECIAL: lo,hi=SPECIAL[rn]; rad=max(theta,0.6)
        else: lo,hi=-350,350; rad=theta
        kO=kern_w(vO,dv,lo,hi); m0=np.tensordot(kO,Tb,axes=(0,0))*dv
        r0=np.hypot(X-px,Y-py)*pix
        if rn in SPECIAL: mpos=(px,py)
        else:
            mm=np.where(r0<0.4,m0,-np.inf); j=np.argmax(mm); mpos=(float(X.flat[j]),float(Y.flat[j]))
        out=dict(pos_pix=mpos,rad=rad)
        # template from HCO+ spectrum in the aperture
        ap=lambda ny_,nx_,pix_,mp: (np.hypot(*(np.mgrid[:ny_,:nx_]-np.array(mp[::-1])[:,None,None])*pix_)<=rad)
        for key,lines in CUBE_LINES.items():
            Tb_,nu_,pix_,bm_,pxy_,free_,cont_=D[key]
            # same sky position in this cube's pixel grid
            a,b=pxy_(pos if rn in SPECIAL else None) if rn in SPECIAL else (None,None)
        # sky position of the peak (convert pixel -> ra,dec via the spw27 wcs: use the offset in arcsec from nominal)
        dx=(mpos[0]-px)*pix; dy=(mpos[1]-py)*pix     # arcsec offsets in pixel axes (x=-RA)
        for mode in ('ap','pix'):
            specs={}
            for key,lines in CUBE_LINES.items():
                Tb_,nu_,pix_,bm_,pxy_,free_,cont_=D[key]
                x0,y0=pxy_(pos); x0+=dx/pix_; y0+=dy/pix_
                ny_,nx_=Tb_.shape[1:]; Y_,X_=np.mgrid[:ny_,:nx_]
                if mode=='ap': apm=(np.hypot(X_-x0,Y_-y0)*pix_<=rad); sp=Tb_[:,apm].mean(axis=1); sp_c=cont_[apm].mean()
                else: ix_,iy_=int(round(x0)),int(round(y0)); sp=Tb_[:,iy_,ix_]; sp_c=cont_[iy_,ix_]
                specs[key]=(nu_,sp,free_,sp_c,float(np.std(sp[free_])))
            o2={}
            for ln,key in (('HCN','X94a_spw27'),('HCOp','X94a_spw27'),('HNC','X94e_spw25'),('HC3N','X94e_spw25'),('H13CN','X94a_spw25'),('H13COp','X94a_spw25'),('SiO','X94a_spw25'),('HN13C','X94a_spw25'),('CCHa','X94a_spw25')):
                nu_,sp,free_,sp_c,sig=specs[key]; v=vel(nu_,F0[ln]); dvc=abs(v[1]-v[0])
                res={}
                for wn,(l2,h2) in {'win':(lo,hi),'h150':(-150,150)}.items():
                    if rn in SPECIAL and wn=='h150': continue
                    k=kern_w(v,dvc,l2,h2)
                    if (k>0).sum()==0: res[wn]=[float('nan')]*4; continue
                    m=float((k*sp).sum()*dvc); e=sig*np.sqrt((k**2).sum())*dvc
                    res[wn]=[m,e,float(v[k>0][0]),float(v[k>0][-1])]
                res['Tpk']=float(np.max(np.convolve(sp,np.ones(3)/3,'same')[np.abs(v)<350])); res['sig_ch']=sig
                o2[ln]=res
            o2['cont_K']={k:specs[k][3] for k in specs}
            nu_,sp,free_,sp_c,sig=specs['X94a_spw25']; nuH,spH,_,_,sigH=specs['X94a_spw27']
            vtO=vel(nuH,F0['HCOp']); order=np.argsort(vtO); tmpl_v=vtO[order]; tmpl=spH[order]
            sel=np.abs(tmpl_v)<(500 if rn not in SPECIAL else 400)
            tm=lambda vv: np.interp(vv,tmpl_v[sel],tmpl[sel],left=0,right=0)/(np.sum(tmpl[sel])*dv)
            lines=['H13CN','H13COp','SiO','HN13C','CCHa','CCHb']
            A=np.array([tm(vel(nu_,F0[ln])) for ln in lines]).T
            use=np.zeros(len(nu_),bool)
            for ln in lines: use|=np.abs(vel(nu_,F0[ln]))<700
            cf=np.linalg.lstsq(A[use],sp[use],rcond=None)[0]
            cov_=np.linalg.inv(A[use].T@A[use])*sig**2
            o2['fit']={ln:[float(cf[i]),float(np.sqrt(cov_[i,i]))] for i,ln in enumerate(lines)}
            o2['fit_resid_rms']=float(np.std((sp-A@cf)[use])); o2['fit_sig_free']=sig
            out[mode]=o2
        out['ap_dummy']=0
        out['spec_HCOp_pix']=None
        out.update({k:v for k,v in out['ap'].items()}) if False else None
        # HCO+ and HCN with same template for consistency
        R[rn]=out
        out=out['ap']
        h=out['HCN']['win'][0]; o=out['HCOp']['win'][0]
        print(theta,rn,'pk(%.0f,%.0f) HCN %.1f+-%.1f HCO+ %.1f+-%.1f R10=%.3f | h150 HCN %.1f HCO+ %.1f HNC %.1f HC3N %.1f | tmpl H13CO+ %.1f+-%.1f SiO %.1f+-%.1f H13CN %.1f+-%.1f HN13C %.1f+-%.1f'%(
            mpos[0],mpos[1],h,out['HCN']['win'][1],o,out['HCOp']['win'][1],h/o,*(out[l].get('h150',out[l]['win'])[0] for l in('HCN','HCOp','HNC','HC3N')),*out['fit']['H13COp'],*out['fit']['SiO'],*out['fit']['H13CN'],*out['fit']['HN13C']),flush=True)
    return R
def xb4d(theta):
    """Xb4d 3-2 at the same circular beam; same apertures and windows; positions: HCO+ 3-2 moment peak"""
    Dd={}
    for ln,(path,f0) in XB4D.items():
        cut,nu,pix,bm,pxy=load(path,half=16); Tb=to_tb_conv(cut,nu,pix,bm,theta,f0*(1-V0/C))*1e3  # path loaded Jy -> mK? tbk uses 1.222e6 -> K
        Tb=Tb/1e3
        v=vel(nu,f0); lf=(v<-450)|(v>600); Tb=Tb-np.median(Tb[lf],axis=0); Dd[ln]=(Tb,v,pix,pxy,lf)
    R={}
    Tb,v,pix,pxy,lf=Dd['HCO+32']; ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    # HCN grid offset ~117 km/s: interpolate HCN onto HCO+ velocity grid
    TH,vH,_,_,_=Dd['HCN32']; fi=np.interp(v,vH[::-1],np.arange(len(vH))[::-1]); i0=np.clip(np.floor(fi).astype(int),0,len(vH)-2); fr=(fi-i0).astype(np.float32)[:,None,None]
    TH=TH[i0]*(1-fr)+TH[i0+1]*fr; dv=abs(v[1]-v[0]); ok=(v>-1400)
    for rn,pos in REG.items():
        px,py=pxy(pos)
        if rn in SPECIAL: lo,hi=SPECIAL[rn]; rad=max(theta,0.6)
        else: lo,hi=-350,350; rad=theta
        out={}
        for wn,(l2,h2) in {'win':(lo,hi),'h150':(-150,150)}.items():
            if rn in SPECIAL and wn=='h150': continue
            k=kern_w(v,dv,l2,h2); mo=np.tensordot(k,Tb,axes=(0,0))*dv; mh=np.tensordot(k,TH,axes=(0,0))*dv
            if wn=='win':
                if rn in SPECIAL: mp=(px,py)
                else:
                    r0=np.hypot(X-px,Y-py)*pix; mm=np.where(r0<0.4,mo,-np.inf); j=np.argmax(mm); mp=(float(X.flat[j]),float(Y.flat[j]))
            ap=(np.hypot(X-mp[0],Y-mp[1])*pix<=rad)
            sO=Tb[:,ap].mean(axis=1); sH=TH[:,ap].mean(axis=1); lfk=(v<-450)|(v>600)
            eO=np.std(sO[lfk])*np.sqrt((k**2).sum())*dv; eH=np.std(sH[lfk])*np.sqrt((k**2).sum())*dv
            out[wn]=dict(HCN32=[float(mh[ap].mean()),float(eH)],HCOp32=[float(mo[ap].mean()),float(eO)])
        R[rn]=out
        o=out['win']; print('Xb4d',theta,rn,'HCN32 %.1f+-%.1f HCO+32 %.1f+-%.1f R32=%.3f'%(*o['HCN32'],*o['HCOp32'],o['HCN32'][0]/o['HCOp32'][0]),flush=True)
    return R
if __name__=='__main__':
    res={}
    for th in THETAS:
        res['%.2f'%th]={'b3':band3(th),'x4d':xb4d(th)}
    json.dump(res,open(OUT+'s182.json','w'),indent=1,default=float)
