# Q35 fix 3: Monte Carlo injection into the beam-shaped deblend of s98/s200 (replaces the ad hoc sqrt(8) error inflation).
# Truth: measured N,S,Bridge amplitudes (+const) of CS 7-6, HCN 4-3, HCO+ 4-3, 346 GHz continuum (out/q5/cs_matched_beam.json, out/q35/s200.json), in the 4.31x2.79" PA87 beam, 16" box,
# noise = white noise smoothed by the beam kernel and scaled to the measured off-source MAD (same estimator as the data). Fit = the s98 linear least squares (N,S,Bridge + const).
# Scenarios: (a) ideal (point sources at the nominal positions); (b) sources resolved (intrinsic 1.0" FWHM, circular) ; (c) positions of the true sources shifted by 0.3" from the fit positions;
# (d) a 4th (Sblue-like) unmodelled component at 25% of N amplitude; each with and without noise-correlation. Output: bias and scatter of N/S ratios, empirical error vs analytic sqrt(8)-inflated error.
# usage: python3 -I s210b_q35_injection.py -> out/q35/s210b.txt/json
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits; from astropy.wcs import WCS; import astropy.units as u
from scipy.ndimage import gaussian_filter
from regions import REGIONS
J=json.load(open('/workspace/out/q5/cs_matched_beam.json'))['common_CSbeam']; S0=json.load(open('/workspace/out/q35/s200.json'))
hd=fits.getheader('/data/member.uid___A001_X1273_Xb47.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'); ps=abs(hd['CDELT2'])*3600
wc=WCS(hd).celestial; ny,nx=hd['NAXIS2'],hd['NAXIS1']; yy,xx=np.mgrid[:ny,:nx]
pos={k:np.array([float(t) for t in wc.world_to_pixel_values(REGIONS[k][0],REGIONS[k][1])]) for k in ('N','S','Bridge','Sblue')}
bmaj,bmin,bpa=4.3122,2.7884,86.676; bm=bmaj/ps; bn=bmin/ps; th=np.deg2rad(bpa+90)
def gb(px,py,extra=0.):   # extra: intrinsic FWHM (arcsec) added in quadrature to the beam axes (resolved source, brightness peak lowered by the area ratio)
    sa=np.sqrt(bm**2+(extra/ps)**2)/2.355; sb=np.sqrt(bn**2+(extra/ps)**2)/2.355; fac=(bm*bn)/(sa*sb*2.355**2)
    dx=xx-px; dy=yy-py; xp=dx*np.cos(th)+dy*np.sin(th); yp=-dx*np.sin(th)+dy*np.cos(th)
    return fac*np.exp(-0.5*((xp/sa)**2+(yp/sb)**2))
mx,my=(pos['N']+pos['S'])/2; box=(abs(xx-mx)*ps<8)&(abs(yy-my)*ps<8)
comps=['N','S','Bridge']
A0=np.array([gb(*pos[c])[box] for c in comps]+[np.ones(box.sum())]).T; Ai=np.linalg.inv(A0.T@A0)@A0.T
# noise smoothing: white noise convolved with an elliptical gaussian of the beam, normalised to unit std
rng=np.random.default_rng(7); sig=np.array([bm,bn])/2.355
def noise(mad):
    w=rng.standard_normal((ny,nx)); n=gaussian_filter(w,(sig[1]/1.0*0+sig[1],sig[0]))  # axes: beam major ~E-W (PA 87): x gets bmaj; small rotation ignored
    return n/ n[box].std()*mad
TR={}   # truth amplitudes
cd=S0['cont']['decomp']['N+S+Bridge']; 
NOISE=sys.argv[1] if len(sys.argv)>1 else 'box'   # 'box': std of the real-fit residual in the 16" box (fit-region noise incl. misfit); 'field': off-source annulus MAD (5-13" from centre, higher primary-beam noise)
LN=['CS 7-6','HCN 4-3','HCO+ 4-3']
for ln in LN: d=J['decomp'][f'{ln}|N+S+Bridge']; TR[ln]=(d['amp'],d['const'],d['rms'] if NOISE=='box' else J['field'][ln]['mad'])
TR['cont']=(cd[0],S0['cont']['field_mad']*0+0.0,S0['cont']['field_mad'])
fk={'CS 7-6':339.7,'HCN 4-3':351.2,'HCO+ 4-3':353.4}
def doubles(a):   # a: dict line-> (N,S,Bridge) fitted amps
    r={}
    f=lambda x,y,ln1,ln2:(x/y)*((fk[ln2]/fk[ln1])**2 if ln1 in fk and ln2 in fk else 1.)
    for nm,(num,den) in {'HCN/CS':('HCN 4-3','CS 7-6'),'HCO+/CS':('HCO+ 4-3','CS 7-6'),'HCN/HCO+':('HCN 4-3','HCO+ 4-3'),'HCN/cont':('HCN 4-3','cont'),'HCO+/cont':('HCO+ 4-3','cont'),'CS/cont':('CS 7-6','cont')}.items():
        rr=[(a[num][i]/a[den][i]) for i in (0,1)]; r[nm]=rr[0]/rr[1]     # N/S double ratio (frequency factors cancel in N/S)
    return r
def truth_map(scn):
    out={}
    for ln,(amp,c0,mad) in TR.items():
        m=np.zeros((ny,nx))
        ext=1.0 if scn=='resolved' else 0.
        for k in comps:
            p=pos[k].copy()
            if scn=='shifted': p=p+np.array([0.3/ps,0.0]) if k=='N' else p-np.array([0.0,0.3/ps]) if k=='S' else p
            m+=amp[k]*gb(*p,extra=ext)
        if scn=='fourth': m+=0.25*amp['N']*gb(*pos['Sblue'])
        out[ln]=m+c0
    return out
A2=np.array([gb(*pos[c])[box] for c in ['N','S']]+[np.ones(box.sum())]).T; Ai2=np.linalg.inv(A2.T@A2)@A2.T
NSIM=600; res={}; T=[]
def P(s): print(s,flush=True); T.append(s)
truth_ratios=None
for scn in ('ideal','resolved','shifted','fourth'):
    tm=truth_map(scn); fits=[]; fits2=[]; amps=[]; singles=[]
    for i in range(NSIM):
        a={}; a2={}
        for ln,(amp,c0,mad) in TR.items():
            y=(tm[ln]+noise(mad))[box]; sol=Ai@y; a[ln]=sol[:3]; a2[ln]=Ai2@y
        singles.append({c:{'CS/HCN':a['CS 7-6'][i]/a['HCN 4-3'][i]*(fk['HCN 4-3']/fk['CS 7-6'])**2,'CS/HCO+':a['CS 7-6'][i]/a['HCO+ 4-3'][i]*(fk['HCO+ 4-3']/fk['CS 7-6'])**2,'HCN/HCO+':a['HCN 4-3'][i]/a['HCO+ 4-3'][i]*(fk['HCO+ 4-3']/fk['HCN 4-3'])**2,'cont':a['cont'][i]} for i,c in enumerate(comps)}); fits.append(doubles(a)); fits2.append(doubles(a2)); amps.append({ln:a[ln] for ln in a})
    if scn=='ideal':
        res['singles']={}
        for c in comps:
            for k in ('CS/HCN','CS/HCO+','HCN/HCO+'):
                v=np.array([x[c][k] for x in singles]); q=np.percentile(v,[16,50,84]); res['singles'][f'{c}|{k}']=q.tolist()
                P(f'  [ideal] single ratio {c} {k} (K km/s): median {q[1]:.3f} [{q[0]:.3f},{q[2]:.3f}] (frac {(q[2]-q[0])/2/q[1]:.0%})')
        for ln in TR:
            v=np.array([x[ln] for x in amps]); P(f'  [ideal, N+S+Bridge fit] {ln}: amplitude scatter N {v[:,0].std():.3g} S {v[:,1].std():.3g} Bridge {v[:,2].std():.3g}; analytic (s98) err '+(str({k:round(x,3) for k,x in J["decomp"][ln+"|N+S+Bridge"]["err"].items()}) if ln!='cont' else str({k:round(x,4) for k,x in cd[1].items()})))
        P('  N+S-only fit to a truth that has a Bridge component: '+'; '.join(f"{k} {np.median([f[k] for f in fits2]):.2f} [{np.percentile([f[k] for f in fits2],16):.2f},{np.percentile([f[k] for f in fits2],84):.2f}]" for k in fits2[0]))
    # input (unbiased) double ratios from the nominal amplitudes
    a0={ln:np.array([TR[ln][0][c] for c in comps]) for ln in TR}; t=doubles(a0)
    P(f'\n== scenario {scn}: recovered N/S double ratios (median, 16-84%, bias=median/input-1) vs input; N={NSIM}')
    for k in t:
        v=np.array([f[k] for f in fits]); q=np.percentile(v,[16,50,84]); res[f'{scn}|{k}']=dict(inp=t[k],q=q.tolist())
        P(f'  {k:10s} input {t[k]:.2f}  recovered {q[1]:.2f} [{q[0]:.2f},{q[2]:.2f}]  bias {q[1]/t[k]-1:+.0%}  frac.scatter {(q[2]-q[0])/2/q[1]:.0%}')
# analytic (s200) double-ratio errors for comparison
dd=S0['tab']['N+S+Bridge']['double']
P('\nAnalytic sqrt(8)-inflated errors (s200, include +-5% systematic): '+'; '.join(f"{k} {v[0]:.2f}+-{v[1]:.2f} ({v[1]/v[0]:.0%})" for k,v in dd.items() if k.endswith('N/S')))
open('/workspace/out/q35/s210b_'+NOISE+'.txt','w').write('\n'.join(T)); json.dump(res,open('/workspace/out/q35/s210b_'+NOISE+'.json','w'),indent=1)
