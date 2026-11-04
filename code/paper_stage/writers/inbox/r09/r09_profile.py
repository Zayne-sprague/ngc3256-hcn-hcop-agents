"""r09: HCN/HCO+ 3-2 T_b ratio vs radius in annuli at N and S, same pipeline as s170/s171 (pbcor cube -> RJ T_b -> flux-conserving
convolution to a circular beam -> per-pixel median baseline of line-free channels (v<-450,>600) -> +-350 km/s moment, HCN put on HCO+ velocity grid).
Runs Xb4b at 0.26" (matched-native) and 1.0", Xb4d at 1.0" (same 1.0" beam => uv/flux cross-check). Annuli centred on the Xb4b 0.26" HCO+ moment peak
(within 0.4" of nominal). Error per annulus = per-pixel sigma(line-free-channel std x sqrt(sum k^2) x dv, includes spatial correlation) x sqrt(1/N_ind)... see code.
usage: python3 -I r09_profile.py -> r09_profile.json"""
import sys, json, numpy as np
from scipy.signal import fftconvolve
def mk(obs):
    sys.argv=['x',obs]; ns={'__name__':'s170'}
    src=open('/workspace/scripts/s170_q26_ratios.py').read(); src=src[:src.index('def main():')]
    exec(src,ns); return ns
ANN=[(0,0.25),(0.25,0.5),(0.5,1.0),(1.0,2.0),(2.0,3.0),(0.5,1.5)]
NUR=(267.5576/265.8864)**2
def maps(ns,theta,peak_world=None):
    DAT={ln:ns['load'](*ns['LINES'][ln]) for ln in ns['LINES']}; sm={}
    for ln,(Tb,v,pix,bm,pxy) in DAT.items():
        k=ns['kernel'](ns['cov'](theta,theta,0),ns['cov'](*bm),pix); c=np.empty_like(Tb)
        for i in range(Tb.shape[0]): c[i]=fftconvolve(np.nan_to_num(Tb[i]),k,'same')
        c-=np.median(c[ns['LF'](v)],axis=0); sm[ln]=(c,v)
    vh=DAT['HCO+'][1]; keep=vh>-1400; cH,vH=sm['HCN']
    fi=np.interp(vh[keep],vH[::-1],np.arange(len(vH))[::-1]); i0=np.clip(np.floor(fi).astype(int),0,len(vH)-2); fr=(fi-i0).astype(np.float32)[:,None,None]
    H=cH[i0]*(1-fr)+cH[i0+1]*fr; O=sm['HCO+'][0][keep]; v=vh[keep]; dv=abs(v[1]-v[0]); lf=ns['LF'](v)
    kw=ns['kern_w'](v,dv,-350,350); sk=np.sqrt((kw**2).sum())*dv
    M={'HCN':np.tensordot(kw,H,axes=(0,0))*dv,'HCO+':np.tensordot(kw,O,axes=(0,0))*dv}
    S={'HCN':np.std(H[lf],axis=0)*sk,'HCO+':np.std(O[lf],axis=0)*sk}
    return M,S,DAT['HCO+']
def run(obs,theta,centres=None):
    ns=mk(obs); M,S,(Tb,v,pix,bm,pxy)=maps(ns,theta); ny,nx=M['HCO+'].shape; Y,X=np.mgrid[:ny,:nx]
    from astropy.io import fits; from astropy.wcs import WCS
    hd=fits.getheader(ns['P'](ns['LINES']['HCO+'][0])); w=WCS(hd).celestial
    n=int(ns['HALF']/pix); cx,cy=[float(t) for t in w.world_to_pixel_values(*ns['CEN'])]; x0=max(int(round(cx))-n,0); y0=max(int(round(cy))-n,0)
    out={}; ba=1.1331*theta**2/pix**2   # circular-Gaussian beam area in pixels
    for rn in ('N','S'):
        if centres is None:
            px,py=pxy(ns['REG'][rn]); r0=np.hypot(X-px,Y-py)*pix; j=np.argmax(np.where(r0<0.4,M['HCO+'],-np.inf)); cpx=(float(X.flat[j]),float(Y.flat[j]))
            cw=[float(t) for t in w.pixel_to_world_values(cpx[0]+x0,cpx[1]+y0)]
        else:
            cw=centres[rn]; a,b=pxy(cw); cpx=(float(a),float(b))
        rr=np.hypot(X-cpx[0],Y-cpx[1])*pix; out[rn]={'centre_world':cw,'rows':{}}
        for lo,hi in ANN:
            m=(rr>=lo)&(rr<hi); nind=max(m.sum()/ba,1.0); d={}
            for ln in ('HCN','HCO+'):
                d[ln]=float(M[ln][m].mean()); d['e'+ln]=float(np.sqrt((S[ln][m]**2).mean())/np.sqrt(nind)) if True else 0
            R=d['HCN']/d['HCO+']; eR=abs(R)*np.hypot(d['eHCN']/d['HCN'],d['eHCO+']/d['HCO+']) if d['HCN']>0 else float('nan')
            d.update(R=R,eR=eR,npix=int(m.sum()),nind=nind); out[rn]['rows'][f'{lo}-{hi}']=d
    return out
if __name__=='__main__':
    res={}
    res['Xb4b_0.26']=run('Xb4b',0.26); cen={k:res['Xb4b_0.26'][k]['centre_world'] for k in 'NS'}
    res['Xb4b_1.0']=run('Xb4b',1.0,cen); res['Xb4d_1.0']=run('Xb4d',1.0,cen)
    json.dump(res,open('/workspace/paper_answers/r09/r09_profile.json','w'),indent=1)
    for k,r in res.items():
        for rn in 'NS':
            print(k,rn,'centre',r[rn]['centre_world'])
            for a,d in r[rn]['rows'].items(): print('  %-8s HCN %6.2f+-%.2f HCO+ %6.2f+-%.2f K km/s  R=%.3f+-%.3f'%(a,d['HCN'],d['eHCN'],d['HCO+'],d['eHCO+'],d['R'],d['eR']))
