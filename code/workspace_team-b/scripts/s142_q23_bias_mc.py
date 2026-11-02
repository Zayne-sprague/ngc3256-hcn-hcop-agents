# Q23: bias of ladder-line integrals and ratios from the low-edge ramp, marginalised over the ramp's unknown extent under the line.
# Ramp model y(ch)=f*cont_loc*max(0,1-ch/c0) (cont_loc from s131, per OUS beam; regions with no s131 continuum (E1,E2) get no ramp -> bias is a lower bound there).
# f fitted per (OUS,spw) at each c0 in GRID on channels v<-450 km/s (line-free) of N,S,Bridge,W; MC over f (Gaussian, fit error), c0 (uniform over GRID).
# Integrals as s30: baseline=median(v<-600 or 450..600), window +-350 km/s, per-channel rms from ref*1.3*sqrt(Nch)*dv; peak-pixel (pix) and aperture (ap) spectra at native beam.
import sys, json, numpy as np; sys.path.insert(0,'/workspace/scripts')
from common import vrad, VSYS
rng=np.random.default_rng(1)
S131=json.load(open('/workspace/out/q22/s131.json'))
F0={'HCN':{2:177.26111e9,3:265.88618e9,4:354.50547e9},'HCO+':{2:178.37506e9,3:267.55763e9,4:356.73424e9}}
LINES={('Xb55',2):((16,'HCN'),(18,'HCO+')),('Xb4f',3):((18,'HCN'),(16,'HCO+')),('Xb47',4):((16,'HCN'),(18,'HCO+'))}
REG=['N','S','Bridge','W']; FITREG=['N','S','Bridge','W']; GRIDS={'narrow':[60,90,120],'wide':[160,200,300,400]}; GRID=sorted(sum(GRIDS.values(),[])); W=350; CAL=0.07
def integ(y,v):
    ref=(v<-600)|((v>450)&(v<600)); dv=abs(np.median(np.diff(v))); b=np.median(y[ref]); m=abs(v)<W
    return (y[m]-b).sum()*dv/1e3, np.std(y[ref]-b)*1.3*np.sqrt(m.sum())*dv/1e3
res={}
for kind in ('pix','ap'):
 for scen,GR in GRIDS.items():
  for (ous,J),lst in LINES.items():
    D={}
    for spw,mol in lst:
        z=np.load(f'/workspace/out/regspec/{ous}.spw{spw}.cube.pbcor.npz'); nu=z['nu']; v=vrad(nu,F0[mol][J])-VSYS; ch=np.arange(len(nu)); dv=abs(np.median(np.diff(v)))
        Y={r:np.nan_to_num(z[f'{r}_{kind}'])*1e3 for r in REG}
        cont={r:S131[f'{ous}.spw{spw}.{r}']['cont'] for r in REG}
        chl=int(np.argmin(abs(v))); fm=ch<chl-int(450/dv)
        if kind=='ap': # ap spectra carry aperture-mean cont; scale from pix s131 cont is approximate -> reuse pix-based f (per cont) 
            pass
        D[mol]=dict(v=v,Y=Y,cont=cont,ch=ch,fm=fm,chl=chl,spw=spw)
    # fits
    for mol,d in D.items():
        d['fit']={}
        for c0 in GRID:
            sh=np.clip(1-d['ch']/c0,0,None)
            num=0;den=0
            for r in FITREG:
                x=d['cont'][r]*sh[d['fm']]; y=d['Y'][r][d['fm']]; num+=np.sum(x*y); den+=np.sum(x*x)
            f=num/den; 
            resid=np.concatenate([d['Y'][r][d['fm']]-f*d['cont'][r]*sh[d['fm']] for r in FITREG]); s=resid.std()*1.3  # 1.3 correlation inflation
            d['fit'][c0]=(f,s/np.sqrt(den),float(np.sum(resid**2)))
    # MC
    N=3000; ci=rng.choice(GR,N); out={}
    for r in REG:
        I={};E={};B={}
        for mol,d in D.items():
            raw,err=integ(d['Y'][r],d['v']); bs=np.empty(N)
            for k in range(N):
                f=rng.normal(d['fit'][ci[k]][0],d['fit'][ci[k]][1]); ramp=f*d['cont'][r]*np.clip(1-d['ch']/ci[k],0,None)
                bs[k]=raw-integ(d['Y'][r]-ramp,d['v'])[0]      # raw - corrected (negative: ramp lowers raw)
            I[mol]=raw;E[mol]=err;B[mol]=bs
        Rraw=I['HCN']/I['HCO+']; Rcorr=(I['HCN']-B['HCN'])/(I['HCO+']-B['HCO+'])
        dR=Rcorr/Rraw-1
        # total 1-sigma on the ratio: stat + 7% cal per line (as in our ratio errors)
        sR=np.sqrt((E['HCN']/I['HCN'])**2+(E['HCO+']/I['HCO+'])**2+2*CAL**2)
        out[r]=dict(I_HCN=I['HCN'],I_HCO=I['HCO+'],stat_HCN=E['HCN'],stat_HCO=E['HCO+'],
          bias_HCN=[float(np.median(B['HCN'])),float(np.percentile(B['HCN'],16)),float(np.percentile(B['HCN'],84))],
          bias_HCO=[float(np.median(B['HCO+'])),float(np.percentile(B['HCO+'],16)),float(np.percentile(B['HCO+'],84))],
          frac_HCN=float(np.median(-B['HCN'])/I['HCN']), frac_HCO=float(np.median(-B['HCO+'])/I['HCO+']),
          R_raw=Rraw, dR_med=float(np.median(dR)), dR_16=float(np.percentile(dR,16)), dR_84=float(np.percentile(dR,84)), sigR_frac=float(sR),
          dR_over_sig=float(np.median(dR)/sR))
    res[f'{kind}.{scen}.{ous}.J{J}']=dict(out=out,fits={m:{str(c):v for c,v in d['fit'].items()} for m,d in D.items()},line_ch={m:d['chl'] for m,d in D.items()})
    print(f'\n=== {kind} {scen} {ous} J={J}-{J-1}  (line ch: {[ (m,d["chl"]) for m,d in D.items()]})')
    for m,d in D.items(): print('  ramp fits',m,'spw',d['spw'],{c:(round(v[0],3),round(v[1],3)) for c,v in d['fit'].items()})
    for r in REG:
        o=out[r]; print(f'  {r:6s} I_HCN={o["I_HCN"]:7.3f}±{o["stat_HCN"]:.3f} bias {o["bias_HCN"][0]:+.3f} [{o["bias_HCN"][1]:+.3f},{o["bias_HCN"][2]:+.3f}] ({100*o["frac_HCN"]:+.1f}%)'
              f' | I_HCO+={o["I_HCO"]:7.3f}±{o["stat_HCO"]:.3f} bias {o["bias_HCO"][0]:+.3f} [{o["bias_HCO"][1]:+.3f},{o["bias_HCO"][2]:+.3f}] ({100*o["frac_HCO"]:+.1f}%)'
              f' | R_raw={o["R_raw"]:.3f} dR={100*o["dR_med"]:+.1f}% [{100*o["dR_16"]:+.1f},{100*o["dR_84"]:+.1f}] vs sigR={100*o["sigR_frac"]:.1f}% -> {o["dR_over_sig"]:+.2f} sigma')
json.dump(res,open('/workspace/out/q23/s142.json','w'),default=float)
