# Q23: bias of ladder-line integrals from the continuum-proportional low-edge ramp.
# 1) fit per (OUS,spw) a ramp model  y(ch)=f*cont_loc*max(0,1-ch/c0)  to line-free low-edge channels of N,S,Bridge,E1,E2,W (peak-pixel spectra)
# 2) integrate each ladder line exactly as s30 (baseline = median of v<-600 and 450..600 km/s about VSYS, window +-W) on raw and ramp-subtracted spectra
# 3) direct test: mean levels on blue and red sides of each line.
import sys, json, numpy as np; sys.path.insert(0,'/workspace/scripts')
from common import vrad, VSYS, C
from scipy.optimize import least_squares
S131=json.load(open('/workspace/out/q22/s131.json')); print(list(S131['Xb47.spw16.N'].keys()) if isinstance(S131['Xb47.spw16.N'],dict) else S131['Xb47.spw16.N'])
F0={'HCN':{2:177.26111e9,3:265.88618e9,4:354.50547e9},'HCO+':{2:178.37506e9,3:267.55763e9,4:356.73424e9}}
LINES=[('Xb55',16,'HCN',2),('Xb55',18,'HCO+',2),('Xb4f',18,'HCN',3),('Xb4f',16,'HCO+',3),('Xb47',16,'HCN',4),('Xb47',18,'HCO+',4),
       ('Xb4d',27,'HCN',3),('Xb4d',25,'HCO+',3),('Xb4b',27,'HCN',3),('Xb4b',25,'HCO+',3)]
REG=['N','S','Bridge','E1','E2','W']
out={}
def spec(ous,spw,kind='pix'):
    z=np.load(f'/workspace/out/regspec/{ous}.spw{spw}.cube.pbcor.npz'); return z['nu'],{r:np.nan_to_num(z[f'{r}_{kind}'])*1e3 for r in REG}
def integ(y,v,W):
    ref=(v<-600)|((v>450)&(v<600)); dv=abs(np.median(np.diff(v)))
    b=np.median(y[ref]); m=abs(v)<W
    sig=np.std(y[ref]-b); return ((y[m]-b).sum()*dv/1e3, sig*1.3*np.sqrt(m.sum())*dv/1e3, b, sig)
for ous,spw,mol,J in LINES:
    nu,Y=spec(ous,spw); v=vrad(nu,F0[mol][J])-VSYS
    n=len(nu); ch=np.arange(n); chl=int(np.argmin(abs(v)))
    dv=abs(np.median(np.diff(v)))
    fit_mask=ch<chl-int(450/dv)    # up to v<-450: beyond 450 km/s blue of the line (line-free)
    cont={r:S131.get(f'{ous}.spw{spw}.{r}',{}).get('cont',np.nan) if isinstance(S131.get(f'{ous}.spw{spw}.{r}'),dict) else np.nan for r in REG}
    # cont: use s131 cont (matched to each OUS beam)
    cvals=np.array([cont[r] for r in REG])
    ok=np.isfinite(cvals)
    def model(p,c,chs): return p[0]*c*np.clip(1-chs/p[1],0,None)
    def res(p):
        return np.concatenate([(Y[r][fit_mask]-model(p,cont[r],ch[fit_mask]))/1.0 for r,o in zip(REG,ok) if o])
    best=None
    for c0 in (40,80,120,160): 
        s=least_squares(res,[-0.1,c0],bounds=([-1,10],[1,400]))
        if best is None or s.cost<best.cost: best=s
    p=best.x; dof=len(best.fun)-2
    J_=best.jac; rms=np.sqrt(2*best.cost/dof); cov=np.linalg.inv(J_.T@J_)*rms**2; pe=np.sqrt(np.diag(cov))
    # bootstrap not needed; also f at fixed c0 grid
    o={'f':p[0],'f_err':pe[0],'c0':p[1],'c0_err':pe[1],'line_ch':chl,'dv':dv}
    print(f'\n{ous} spw{spw} {mol} {J}-{J-1}: line at ch {chl}, ramp fit f={p[0]:.3f}+-{pe[0]:.3f} x cont, c0={p[1]:.0f}+-{pe[1]:.0f} ch (rms resid {rms:.2f})')
    # alt ramp: fit only f with c0 free is above; alternative extrapolation: ramp zero beyond measured range = same model (linear to c0)
    for W in (250,350,500):
        for r in REG:
            y=Y[r]; ramp=model(p,cont[r],ch) if np.isfinite(cont[r]) else 0*y
            a=integ(y,v,W); b=integ(y-ramp,v,W)
            # side levels
            blue=np.mean(y[(v>-600)&(v<-420)]); red=np.mean(y[(v>420)&(v<600)])
            o[f'{r}_W{W}']=dict(raw=a[0],corr=b[0],bias=a[0]-b[0],err=a[1],base_raw=a[2],base_corr=b[2],blue=blue,red=red,rms=a[3])
    out[f'{ous}.spw{spw}.{mol}{J}']=o
    for W in (350,):
        for r in REG:
            q=o[f'{r}_W{W}']; print(f'  W{W} {r:6s} I_raw={q["raw"]:8.3f} I_corr={q["corr"]:8.3f} bias={q["bias"]:+.3f} ({100*q["bias"]/q["corr"] if q["corr"]!=0 else 0:+.1f}%) stat={q["err"]:.3f} bias/stat={q["bias"]/q["err"]:+.2f}  baseline raw {q["base_raw"]:+.2f} corr {q["base_corr"]:+.2f}; side levels blue(-600..-420) {q["blue"]:+.2f} red(420..600) {q["red"]:+.2f} mJy rms/ch {q["rms"]:.2f}')
json.dump(out,open('/workspace/out/q23/s141.json','w'),default=float)
