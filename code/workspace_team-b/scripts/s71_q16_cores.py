# Q16 (1): independent test of team-a's p0012 N/S deblended core ladders (q1_deblend_ladders_xb4f.csv: S HCN 4-3/3-2 Tb 0.268, N 0.395; HCO+ S 0.442, N 0.482).
# Forward model with OUR aperture machinery: Xb4d native 3-2 map (0.94x0.74") = morphology template, split into core N, core S (r<rc of each nucleus) and the rest;
# each component smoothed to the 6.5" and 10" circular beams (s30 maps/ratios use the same beams); fractions of the 3-2 flux of each component at the six s30 positions
# (N,S,Bridge,E1,E2,W). Observed T43/T32 (K km/s ratio, s30_ratios.json, same beam, same position) = sum_c f_c r_c  -> solve r_N, r_S, r_rest (and a 4-comp variant with a Bridge core).
# v2: fractions are referenced to the OBSERVED Xb4f 3-2 flux at the position (rest = Xb4f - cores, clipped at 0) because Xb4d loses extended flux; template clipped at 0; only N,S,Bridge (+E1,E2 at 10") are used.
# MC over rc (1.0-2.5"), position shift (0.3"), template floor 3-5% (+ stat), position-dependent rest ladder gradient.   usage: python3 -I s71_q16_cores.py  -> out/q16/s71.txt/json
import json,numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import nnls
import sys; sys.path.insert(0,'/workspace/scripts')
from regions import REGIONS
OUT='/workspace/out/q16/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
d=np.load('/workspace/out/s30_maps_Xb4d_J3_native.npz'); ps=float(d['ps']); cx,cy=float(d['cx']),float(d['cy']); bn=d['beam']
from common import NUC
RA0=0.5*(NUC['N'][0]+NUC['S'][0]); DE0=0.5*(NUC['N'][1]+NUC['S'][1]); cd=np.cos(np.radians(DE0))
def off(r): return (r[0]-RA0)*cd*3600,(r[1]-DE0)*3600
POS={k:off(REGIONS[k]) for k in ('N','S','Bridge','E1','E2','W')}
R={r['tag']:r for r in json.load(open('/workspace/out/s30_ratios.json'))}
F={'HCN':{2:177.26111,3:265.88618,4:354.50547},'HCO+':{2:178.37506,3:267.55763,4:356.73424}}
def pix(dx,dy,sgn): return cx+sgn*dx/ps,cy+dy/ps
# orientation check: HCN peak near N position
M=np.nan_to_num(d['HCN']); Y,X=np.mgrid[:M.shape[0],:M.shape[1]]
for sgn in (-1,1):
    x,y=pix(*POS['N'],sgn); w=np.hypot(X-x,Y-y)<4/ps; i=np.argmax(np.where(w,M,0)); P(f'# orientation sgn={sgn}: max HCN within 4" of N at pixel dist {np.hypot(X.ravel()[i]-x,Y.ravel()[i]-y)*ps:.2f}" ; value {M.ravel()[i]:.2f}')
SG=-1
def obs_ratio(mol,fw,Jup):
    out={}
    key='HCN' if mol=='HCN' else 'HCOp'
    a=R[f'Xb47_J4_{fw}as']['rows']; b=R[f'Xb4f_J3_{fw}as']['rows']
    for p in POS:
        S4,S3=a[p+'_pos'][key],b[p+'_pos'][key]; e4,e3=a[p+'_pos']['e'+key],b[p+'_pos']['e'+key]
        # Jy/beam km/s at same circular beam -> K ratio = (S4/S3) * (nu3/nu4)^2 ; line-free noise as stat
        k=(F[mol][3]/F[mol][4])**2; r=S4/S3*k; out[p]=(r,r*np.hypot(e4/S4,e3/S3))
    return out
def frac(mol,fw,rc,dshift=(0,0),comps=('N','S','rest'),sgn=SG):
    key='HCN' if mol=='HCN' else 'HCOp'
    M=np.clip(np.nan_to_num(d[key]),0,None); Y,X=np.mgrid[:M.shape[0],:M.shape[1]]
    msk={}
    for c in ('N','S','Bridge'):
        x,y=pix(*POS[c],sgn); x+=dshift[0]/ps*sgn; y+=dshift[1]/ps; msk[c]=np.hypot(X-x,Y-y)*ps<rc[c]
    C={'N':np.where(msk['N'],M,0),'S':np.where(msk['S'],M,0)}
    if 'Bridge' in comps: C['Bridge']=np.where(msk['Bridge']&~msk['N']&~msk['S'],M,0)
    sig=np.sqrt(max(fw**2-(bn[0]*bn[1]),1e-3))/2.3548/ps; fac=fw*fw/(bn[0]*bn[1])
    Sm={k:gaussian_filter(v,sig)*fac for k,v in C.items()}
    ref=R['Xb4f_J3_%sas'%('6.5' if fw==6.5 else '10')]['rows']
    out={}
    for p in POS:
        x,y=pix(*POS[p],sgn); v={k:float(m[int(round(y)),int(round(x))]) for k,m in Sm.items()}
        tot=ref[p+'_pos'][key]; f={k:v[k]/tot for k in v}; f['rest']=max(1-sum(f.values()),0.0); out[p]=f
    return out
def solve(mol,rc,dshift,floor,comps,rng=None,use=('6.5','10'),pos_use=None):
    rows=[];y=[];s=[]
    for fw in (6.5,10):
        if str(fw).rstrip('0').rstrip('.') not in use and str(fw) not in use: continue
        ob=obs_ratio(mol,'6.5' if fw==6.5 else '10',4); fr=frac(mol,fw,rc,dshift,comps)
        for p in POS:
            if p in ('E1','E2','W') and (fw==6.5 or p=='W'): continue
            rows.append([fr[p][c] for c in comps]); y.append(ob[p][0]); s.append(np.hypot(ob[p][0]*floor,ob[p][1]))
    A=np.array(rows);y=np.array(y);s=np.array(s)
    if rng is not None: y=y*(1+floor*rng.standard_normal(len(y)))
    r,rn=nnls(A/s[:,None],y/s); return r,float(np.sum(((A@r-y)/s)**2)),len(y)-len(comps),A,y
def obs10(mol):
    return obs_ratio(mol,'10',4)
P('# observed T43/T32 (K ratio, same beam) at the six positions')
for mol in ('HCN','HCO+'):
    for fw in ('6.5','10'): P(f'  {mol:5s} {fw:>4s}": '+'  '.join(f'{p} {v[0]:.3f}' for p,v in obs_ratio(mol,fw,4).items()))
res={}
TEAMA={'HCN':{'N':0.395,'S':0.268},'HCO+':{'N':0.482,'S':0.442}}
for comps in (('N','S','rest'),('N','S','Bridge','rest')):
    P(f'\n# model components {comps}')
    for mol in ('HCN','HCO+'):
        rng=np.random.default_rng(1); sol=[];chis=[]
        for it in range(600):
            rcv=rng.uniform(1.0,2.5); rc={'N':rcv,'S':rcv,'Bridge':rng.uniform(0.8,1.5)}; sh=tuple(rng.normal(0,0.3,2)); fl=rng.uniform(0.03,0.05)
            r,c2,dof,A,y=solve(mol,rc,sh,fl,comps,rng)
            sol.append(r);chis.append(c2/max(dof,1))
        sol=np.array(sol); q=np.percentile(sol,[16,50,84],axis=0)
        r0,c0,dof,A,y=solve(mol,{'N':1.5,'S':1.5,'Bridge':1.0},(0,0),0.04,comps)
        P(f'  {mol:5s} nominal(rc=1.5") '+' '.join(f'r_{c}={v:.3f}' for c,v in zip(comps,r0))+f'  chi2/dof {c0:.1f}/{dof}')
        P(f'        MC 16/50/84: '+'  '.join(f'r_{c} {q[1,i]:.3f} [{q[0,i]:.3f},{q[2,i]:.3f}]' for i,c in enumerate(comps))+f'   median chi2/dof {np.median(chis):.2f}')
        P(f'        team-a (p0012) core ladders: N {TEAMA[mol]["N"]}  S {TEAMA[mol]["S"]};  S core: ours {q[1,1]:.3f} vs team-a {TEAMA[mol]["S"]} -> ratio {q[1,1]/TEAMA[mol]["S"]:.2f}')
        res[f'{"+".join(comps)}|{mol}']=dict(nominal=dict(zip(comps,map(float,r0))),mc=q.tolist(),chi2dof=float(np.median(chis)))
# consistency of team-a's S core with the S 6.5" aperture: what rest ladder is implied?
P('\n# implied ladder of the non-core gas if the cores have the team-a values (fix r_N,r_S; solve r_rest), per position set')
for mol in ('HCN','HCO+'):
    for rcv in (1.0,1.5,2.5):
        rc={'N':rcv,'S':rcv,'Bridge':1.0}; out=[]
        for fw in ('6.5','10'):
            ob=obs_ratio(mol,fw,4); fr=frac(mol,float(fw),rc)
            for p in ('N','S','Bridge','E1','E2'):
                if fw=='6.5' and p in ('E1','E2'): continue
                f=fr[p]; imp=(ob[p][0]-f['N']*TEAMA[mol]['N']-f['S']*TEAMA[mol]['S'])/f['rest']; out.append(f'{p}@{fw}: f_S={f["S"]:.2f} f_N={f["N"]:.2f} rest_impl={imp:.3f}')
        P(f'  {mol:5s} rc={rcv}": '+'; '.join(out))
json.dump(res,open(OUT+'s71.json','w'),indent=1); open(OUT+'s71.txt','w').write('\n'.join(log))
