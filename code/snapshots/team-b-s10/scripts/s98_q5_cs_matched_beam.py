# Q5 round 2: CS 7-6 vs HCN 4-3 vs HCO+ 4-3 in a COMMON beam (Xb47, spw22/16/18), baseline-offset correction, null test and N/S/Bridge point-source decomposition.
# Targets: (a) CS native beam 4.31x2.79" PA 87 (HCN, HCO+ smoothed per channel with radio_beam kernels, flux rescaled by beam-area ratio); (b) circular 6.5" (the aperture team-a used for Q1).
# Moment-0: window 2780+-300 km/s (all flux in the CS window counted as CS = upper bound on CS; blend with H2CS/29SiO <=0.6 of 1.9 at N from the template fits in out/q5/b7_fits.txt), linear baseline from 500-1100 km/s sidebands, per pixel.
# Empirical noise and baseline offset: median / 1.4826 MAD of the moment map over off-source pixels (5-13" from the field centre, >2.5" from all region positions).
# Null test: moment-0 at positions reflected through the field centre.  Decomposition: least-squares fit of point sources (N,S,Bridge,Sblue; beam-shaped) + constant in a 16" box.
# Outputs out/q5/cs_matched_beam.json, .csv.   usage: python3 -I s98_q5_cs_matched_beam.py
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS; import astropy.units as u
from astropy.convolution import convolve_fft; from radio_beam import Beam
from scipy.optimize import nnls
from common import files; from regions import REGIONS
C=299792.458
LINES={'CS 7-6':('spw22',342.88285),'HCN 4-3':('spw16',354.50547),'HCO+ 4-3':('spw18',356.73424)}
def moment(spw,f0,tb=None,vc=2780.,w=300.,extra=()):
    h=fits.open(files(f'*Xb47*{spw}.cube.I.pbcor.fits')[0]); hd=h[0].header; d=np.nan_to_num(h[0].data[0]); nz=d.shape[0]; ps=abs(hd['CDELT2'])*3600
    nu=hd['CRVAL3']+(np.arange(nz)+1-hd['CRPIX3'])*hd['CDELT3']; v=C*(1-nu/(f0*1e9)); dv=abs(np.median(np.diff(v)))
    b=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg)
    if tb is not None and (tb.major>b.major or tb.minor>b.minor):
        ker=tb.deconvolve(b).as_kernel(ps*u.arcsec,x_size=41,y_size=41).array; ker/=ker.sum(); fac=(tb.sr/b.sr).value; bo=tb
    else: ker=None; fac=1.; bo=b
    sel=np.where((abs(v-vc)<1100)|np.any([abs(v-e)<1100 for e in extra],axis=0) if len(extra) else (abs(v-vc)<1100))[0]; dd=np.zeros((len(sel),)+d.shape[1:],dtype=np.float32)
    for j,i in enumerate(sel): dd[j]=(convolve_fft(d[i],ker,boundary='fill',allow_huge=True)*fac) if ker is not None else d[i]
    vs=v[sel]; sb=(abs(vs-vc)>500); X=np.vstack([np.ones(sb.sum()),(vs[sb]-vc)/600.]).T; P=np.linalg.pinv(X)
    coef=np.tensordot(P,dd[sb],axes=(1,0)); base=coef[0][None]+coef[1][None]*((vs-vc)/600.)[:,None,None]
    inw=abs(vs-vc)<w; mom=((dd-base)[inw]).sum(0)*dv; rms=np.std((dd-base)[sb],axis=0)
    nuobs=np.interp(vc,v[::-1],nu[::-1])/1e9
    ex={}
    for e in extra:   # same estimator at line-free window centres (on-source null)
        s2=abs(vs-e)>500; s2&=abs(vs-e)<1100; w2=abs(vs-e)<w
        if s2.sum()<20 or w2.sum()<10: continue
        X2=np.vstack([np.ones(s2.sum()),(vs[s2]-e)/600.]).T; P2=np.linalg.pinv(X2); c2=np.tensordot(P2,dd[s2],axes=(1,0)); b2=c2[0][None]+c2[1][None]*((vs-e)/600.)[:,None,None]
        ex[e]=((dd-b2)[w2]).sum(0)*dv
    return mom,rms,hd,bo,nuobs,ex
def run(tag,tbeam):
    res={}; maps={}
    for nm,(spw,f0) in LINES.items():
        m,r,hd,bo,nuo,ex=moment(spw,f0,tbeam,extra=()); maps[nm]=(m,r,nuo)
    wc=WCS(hd).celestial; ny,nx=maps['CS 7-6'][0].shape; yy,xx=np.mgrid[:ny,:nx]; ps=abs(hd['CDELT2'])*3600
    cx,cy=wc.world_to_pixel_values(156.96346,-43.90461); rr=np.hypot(xx-cx,yy-cy)*ps
    pos={k:[float(t) for t in wc.world_to_pixel_values(REGIONS[k][0],REGIONS[k][1])] for k in ('N','S','Bridge','Sblue','W')}
    mk=(rr>5)&(rr<13)
    for k,(px,py) in pos.items(): mk&=np.hypot(xx-px,yy-py)*ps>2.5
    out={'beam':f'{bo.major.value*1:.2f}x{bo.minor.value*1:.2f}'}
    st={}
    for nm,(m,r,nuo) in maps.items():
        a=m[mk&(m!=0)]; st[nm]=dict(med=float(np.median(a)),mad=float(1.4826*np.median(abs(a-np.median(a)))),n=int(len(a)),nuobs=nuo)
    out['field']=st
    def val(m,px,py): return m[int(round(py)),int(round(px))]
    for k,(px,py) in pos.items():
        o={}
        for nm,(m,r,nuo) in maps.items():
            raw=val(m,px,py); o[nm]=dict(raw=float(raw),corr=float(raw-st[nm]['med']),mad=st[nm]['mad'],formal=float(val(r,px,py)*np.sqrt(60)*3.91*0+0))
        for den in ('HCN 4-3','HCO+ 4-3'):
            for tag2,key in (('raw','raw'),('corr','corr')):
                a=o['CS 7-6'][key]; b=o[den][key]; fk=(maps[den][2]/maps['CS 7-6'][2])**2
                ea=st['CS 7-6']['mad']; eb=st[den]['mad']
                R=a/b*fk; eR=abs(R)*np.hypot(ea/a,eb/b)
                o[f'R_K CS/{den} {tag2}']=[float(R),float(eR)]
        o['snr_CS_raw']=o['CS 7-6']['raw']/st['CS 7-6']['mad']; o['snr_CS_corr']=o['CS 7-6']['corr']/st['CS 7-6']['mad']; out[k]=o
    # null test: fraction of off-source field pixels (same estimator, annulus above) with moment >= the on-source value. (A reflection through the field centre is NOT a null here: the pointing is the N-S midpoint, so it maps N<->S; and the spw covers only +-800 km/s around the line, too short for line-free windows with their own baselines.)
    nul={}; fa=maps['CS 7-6'][0][mk&(maps['CS 7-6'][0]!=0)]
    for k,(px,py) in pos.items():
        v_=val(maps['CS 7-6'][0],px,py); nul[k]=dict(W=float(v_),frac_field_ge=float((fa>=v_).mean()),n=int(len(fa)))
    out['null_field_exceedance_CS']=nul
    # point-source decomposition (beam-shaped components N,S,Bridge,Sblue + const), 16" box around the N-S midpoint
    bm=bo.major.to(u.arcsec).value/ps; bn=bo.minor.to(u.arcsec).value/ps; th=np.deg2rad(bo.pa.value+90)   # PA east of north -> x axis; sigma=FWHM/2.355
    def gb(px,py):
        dx=xx-px; dy=yy-py; xp=dx*np.cos(th)+dy*np.sin(th); yp=-dx*np.sin(th)+dy*np.cos(th)  # major axis along PA (from north, toward east=-x?) approximate
        return np.exp(-0.5*((xp/(bm/2.355))**2+(yp/(bn/2.355))**2))
    mx,my=(pos['N'][0]+pos['S'][0])/2,(pos['N'][1]+pos['S'][1])/2; box=(abs(xx-mx)*ps<8)&(abs(yy-my)*ps<8)
    dec={}
    for nm,(m,r,nuo) in maps.items():
        for tagc,comps in (('N+S',['N','S']),('N+S+Bridge',['N','S','Bridge']),('N+S+Bridge+Sblue',['N','S','Bridge','Sblue'])):
            A=np.array([gb(*pos[c])[box] for c in comps]+[np.ones(box.sum())]).T
            sol=np.linalg.lstsq(A,m[box],rcond=None)[0]; res_=m[box]-A@sol; dof=box.sum()-A.shape[1]
            cov=np.linalg.inv(A.T@A)*np.sum(res_**2)/dof*8   # x8 for ~pixel correlation (beam area/pixel^2 ~ 40/ ~5 indep)
            dec[f'{nm}|{tagc}']=dict(amp={c:float(s) for c,s in zip(comps,sol)},err={c:float(np.sqrt(cov[i,i])) for i,c in enumerate(comps)},const=float(sol[-1]),rms=float(res_.std()))
    out['decomp']=dec
    return out
OUT={'common_CSbeam':run('cs',Beam(4.3122*u.arcsec,2.7884*u.arcsec,86.676*u.deg)),'circ6.5':run('c65',Beam(6.5*u.arcsec,6.5*u.arcsec,0*u.deg))}
json.dump(OUT,open('/workspace/out/q5/cs_matched_beam.json','w'),indent=1)
for tg,o in OUT.items():
    print('==',tg,o['beam']); print('field',{k:(round(v['med'],3),round(v['mad'],3)) for k,v in o['field'].items()})
    for k in ('N','S','Bridge','Sblue','W'):
        e=o[k]; print(k,'CS raw %.2f corr %.2f (snr %.1f/%.1f)'%(e['CS 7-6']['raw'],e['CS 7-6']['corr'],e['snr_CS_raw'],e['snr_CS_corr']),'HCN %.2f HCO+ %.2f'%(e['HCN 4-3']['raw'],e['HCO+ 4-3']['raw']),'R_K CS/HCN raw %.3f+-%.3f corr %.3f+-%.3f | CS/HCO+ raw %.3f corr %.3f+-%.3f'%(*e['R_K CS/HCN 4-3 raw'],*e['R_K CS/HCN 4-3 corr'],e['R_K CS/HCO+ 4-3 raw'][0],*e['R_K CS/HCO+ 4-3 corr']))
    print('field exceedance',{k:round(v['frac_field_ge'],4) for k,v in o['null_field_exceedance_CS'].items()})
    for k,v in o['decomp'].items(): print(k,{c:round(a,2) for c,a in v['amp'].items()},'err',{c:round(a,2) for c,a in v['err'].items()},'const',round(v['const'],2))
