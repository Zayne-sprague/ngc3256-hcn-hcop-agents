# Q35: abundance-normalised comparison HCN vs HCO+ against CS 7-6 and the 346 GHz continuum at N, S, Bridge (+disk limits), one common beam (4.31x2.79" PA 87, Xb47).
# Inputs: out/q5/cs_matched_beam.json (s98: line moments, common beam, beam-shaped point-source decomposition N,S,Bridge(+Sblue) + const in a 16" box; errors x sqrt(8)),
#         Xb47 continuum (spw16_18_20_22.cont, 4.16x2.70", convolved to the CS beam), same decomposition. Line amplitudes: Jy/beam km/s; continuum Jy/beam.
# Outputs: out/q35/s200.json, s200.txt.  usage: python3 -I s200_q35_cs_norm.py
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, os, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS; import astropy.units as u
from astropy.convolution import convolve_fft; from radio_beam import Beam
from regions import REGIONS
os.makedirs('/workspace/out/q35',exist_ok=True)
J=json.load(open('/workspace/out/q5/cs_matched_beam.json'))['common_CSbeam']
L=('CS 7-6','HCN 4-3','HCO+ 4-3'); tb=Beam(4.3122*u.arcsec,2.7884*u.arcsec,86.676*u.deg)
h=fits.open('/data/member.uid___A001_X1273_Xb47.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'); hd=h[0].header; d=np.nan_to_num(h[0].data[0,0]); ps=abs(hd['CDELT2'])*3600
b=Beam(hd['BMAJ']*u.deg,hd['BMIN']*u.deg,hd['BPA']*u.deg)
ker=tb.deconvolve(b).as_kernel(ps*u.arcsec,x_size=41,y_size=41).array; ker/=ker.sum(); m=convolve_fft(d,ker,boundary='fill',allow_huge=True)*(tb.sr/b.sr).value
wc=WCS(hd).celestial; ny,nx=m.shape; yy,xx=np.mgrid[:ny,:nx]
cx,cy=wc.world_to_pixel_values(156.96346,-43.90461); rr=np.hypot(xx-cx,yy-cy)*ps
pos={k:[float(t) for t in wc.world_to_pixel_values(REGIONS[k][0],REGIONS[k][1])] for k in ('N','S','Bridge','Sblue','E1','E2','W')}
mk=(rr>5)&(rr<13)
for k,(px,py) in pos.items(): mk&=np.hypot(xx-px,yy-py)*ps>2.5
a=m[mk&(m!=0)]; med=np.median(a); mad=1.4826*np.median(abs(a-med))
bm=tb.major.value/ps; bn=tb.minor.value/ps; th=np.deg2rad(tb.pa.value+90)
def gb(px,py):
    dx=xx-px; dy=yy-py; xp=dx*np.cos(th)+dy*np.sin(th); yp=-dx*np.sin(th)+dy*np.cos(th)
    return np.exp(-0.5*((xp/(bm/2.355))**2+(yp/(bn/2.355))**2))
mx,my=(pos['N'][0]+pos['S'][0])/2,(pos['N'][1]+pos['S'][1])/2; box=(abs(xx-mx)*ps<8)&(abs(yy-my)*ps<8)
dec={}
for tag,comps in (('N+S',['N','S']),('N+S+Bridge',['N','S','Bridge'])):
    A=np.array([gb(*pos[c])[box] for c in comps]+[np.ones(box.sum())]).T; sol=np.linalg.lstsq(A,m[box],rcond=None)[0]; r_=m[box]-A@sol
    cov=np.linalg.inv(A.T@A)*np.sum(r_**2)/(box.sum()-A.shape[1])*8; dec[tag]=({c:float(s) for c,s in zip(comps,sol)},{c:float(np.sqrt(cov[i,i])) for i,c in enumerate(comps)})
cont={'field_med':float(med),'field_mad':float(mad),'decomp':dec,'peak':{k:float(m[int(round(py)),int(round(px))]) for k,(px,py) in pos.items()}}
out={'cont':cont}; T=[]
def P(s): print(s); T.append(s)
P('Continuum 346 GHz in the 4.31x2.79" beam (Jy/beam): decomposition N+S+Bridge '+str({k:round(v,4) for k,v in dec['N+S+Bridge'][0].items()})+' +- '+str({k:round(v,4) for k,v in dec['N+S+Bridge'][1].items()})+f' ; field med {med:.4f} MAD {mad:.4f}')
P('Peak pixels: '+str({k:round(v,4) for k,v in cont['peak'].items()}))
def ratio(a,ea,b,eb): r=a/b; return r,abs(r)*np.hypot(ea/a,eb/b)
tab={}
for tag in ('N+S','N+S+Bridge'):
    cd,ce=dec[tag]; comps=list(cd); P(f'\n== decomposition {tag} (deblended beam-shaped components)')
    amp={ln:(J['decomp'][f'{ln}|{tag}']['amp'],J['decomp'][f'{ln}|{tag}']['err']) for ln in L}
    amp['cont']=(cd,ce); rows={}
    for c in comps:
        r={}
        for ln in L: r[ln]=(amp[ln][0][c],amp[ln][1][c])
        r['cont']=(cd[c],ce[c])
        fk43={'CS 7-6':339.7,'HCN 4-3':351.2,'HCO+ 4-3':353.4}
        def K(num,den): # RJ brightness ratio of integrated lines (nu^-2 conversion)
            x,ex=ratio(*r[num],*r[den]); f=(fk43[den]/fk43[num])**2 if num in fk43 else 1.; return x*f,ex*f
        rows[c]={'CS/HCN':K('CS 7-6','HCN 4-3'),'CS/HCO+':K('CS 7-6','HCO+ 4-3'),'HCN/CS':tuple(1/np.array(K('CS 7-6','HCN 4-3')[:1]))+(0,),
                 'HCN/HCO+':K('HCN 4-3','HCO+ 4-3'),'CS/cont':ratio(*r['CS 7-6'],*r['cont']),'HCN/cont':ratio(*r['HCN 4-3'],*r['cont']),'HCO+/cont':ratio(*r['HCO+ 4-3'],*r['cont'])}
        rows[c]['HCN/CS']=(1/rows[c]['CS/HCN'][0],rows[c]['CS/HCN'][1]/rows[c]['CS/HCN'][0]**2)
        rows[c]['HCO+/CS']=(1/rows[c]['CS/HCO+'][0],rows[c]['CS/HCO+'][1]/rows[c]['CS/HCO+'][0]**2)
        P(f'{c}: '+'; '.join(f'{k} {v[0]:.3g}+-{v[1]:.2g}' for k,v in rows[c].items()))
    tab[tag]=rows
    # double ratios relative to S, with +-0.05 fractional systematic added in quadrature
    for ref in ('S',):
        for c in comps:
            if c==ref: continue
            s=[]
            for k in ('HCN/CS','HCO+/CS','HCN/HCO+','CS/cont','HCN/cont','HCO+/cont'):
                x=rows[c][k][0]/rows[ref][k][0]; e=x*np.sqrt((rows[c][k][1]/rows[c][k][0])**2+(rows[ref][k][1]/rows[ref][k][0])**2+0.05**2)
                s.append(f'{k} {c}/{ref} = {x:.2f}+-{e:.2f}'); tab[tag].setdefault('double',{})[f'{k} {c}/{ref}']=(float(x),float(e))
            P('  '+'; '.join(s))
# disk limits (peak pixels, offset-corrected amplitudes, noise = field MAD) E1, E2, W
P('\n== Disk/off-nuclear clumps (peak pixel, 4.31x2.79", offset-corrected; 3-sigma limits when <3 sigma)')
dk={}
for k in ('W','Sblue'):
    e=J[k]; cs,sc=e['CS 7-6']['corr'],J['field']['CS 7-6']['mad']; hcn=e['HCN 4-3']['corr']; hco=e['HCO+ 4-3']['corr']; lim=3*sc
    P(f"{k}: CS {cs:.2f}+-{sc:.2f} (3sig limit {lim:.2f}); HCN {hcn:.2f}+-{J['field']['HCN 4-3']['mad']:.2f}; HCO+ {hco:.2f}+-{J['field']['HCO+ 4-3']['mad']:.2f} Jy/beam km/s; CS/HCO+ (flux) < {lim/hco:.2f} (3 sigma) vs N {J['N']['CS 7-6']['corr']/J['N']['HCO+ 4-3']['corr']:.3f}")
    dk[k]=dict(cs=cs,sig=sc,hcn=hcn,hco=hco)
out['tab']=tab; out['disk']=dk
json.dump(out,open('/workspace/out/q35/s200.json','w'),indent=1,default=float); open('/workspace/out/q35/s200.txt','w').write('\n'.join(T))
