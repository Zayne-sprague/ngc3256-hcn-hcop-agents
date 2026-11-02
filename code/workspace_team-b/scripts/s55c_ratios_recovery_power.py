# Q3 round 2: (1) empirical flux-recovery calibration from the nested Xb4b/Xb4d/Xb4f set; (2) same-band HCN/HCO+ ratios and their J-dependence per position with full noise,
# centre-vs-disk significance (stat only / +7% per line independent / calibration common), differential-recovery needed to erase the disk trend;
# (3) dilution power test for N vs S vs Bridge at 10" (two-zone injection, overlapping-aperture noise covariance).
# usage: python3 -I s55c_ratios_recovery_power.py -> out/q3b/s55c_*.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s55_lib as S
L=S.L; log=[]
def P(s=''): print(s); log.append(s)
POS=['N','S','Bridge','E1','E2','W']
# ---------------------------------------------------------------- (1) recovery calibration
R26=json.load(open('/workspace/out/s26_recovery.json'))
BEAM={'Xb4b':0.24,'Xb4d':0.83,'Xb4f':4.5,'Xb47':3.4,'Xb55':7.3}   # geometric-mean beam FWHM ["] (native)
P('# (1) Recovery relative to Xb4f (same line, convolved to the Xb4f beam 5.7x3.5", |v|<350 km/s, s26): Xb4d (0.83") and Xb4b (0.24")')
P('pos      mol   Xb4d/Xb4f  Xb4b/Xb4f   slope s=-d ln R/d log2(beam) (from the 2 points)')
rec={}
for pos in POS:
    for mol in ('HCN','HCO+'):
        k=f'pk_{pos}_core'; f=R26[mol]['Xb4f'][k]; d=R26[mol]['Xb4d'][k]/f; b=R26[mol]['Xb4b'][k]/f
        x=np.log2(np.array([BEAM['Xb4b'],BEAM['Xb4d']])/BEAM['Xb4f']); y=np.log(np.array([b,d])); s=float(np.sum(x*y)/np.sum(x*x))   # ln R = s log2(b/b_f)
        rec[(pos,mol)]=(d,b,s); P(f'{pos:7s}  {mol:5s}  {d:6.2f}     {b:6.2f}       s={s:.3f}')
P('\n# Expected relative recovery vs Xb4f if MRS scales with beam (ln R = s log2(b/b_Xb4f)): Xb47 (3.4"): R47/R4f = exp(s log2(3.4/4.5)); Xb55 (7.3"): exp(s log2(7.3/4.5)). (Band 5/7 baseline sets differ; treated as k in 0.5-1.5.)')
for pos in POS:
    sH=rec[(pos,'HCN')][2]; sO=rec[(pos,'HCO+')][2]
    P(f'{pos:7s} HCN: R47/R4f={np.exp(sH*np.log2(3.4/4.5)):.2f} R55/R4f={np.exp(sH*np.log2(7.3/4.5)):.2f} | HCO+: R47/R4f={np.exp(sO*np.log2(3.4/4.5)):.2f} R55/R4f={np.exp(sO*np.log2(7.3/4.5)):.2f} | HCN/HCO+ relative recovery at Xb4d beam = {rec[(pos,"HCN")][0]/rec[(pos,"HCO+")][0]:.2f}, Xb4b {rec[(pos,"HCN")][1]/rec[(pos,"HCO+")][1]:.2f}')
# ---------------------------------------------------------------- (2) ratios
P('\n# (2) Same-band HCN/HCO+ (K km/s, RJ) at 10" per position, with errors: stat | +7% per line independent per position (conservative) | calibration common to all positions (cancels in centre/disk), + PB term 0.3(1-PB) per band')
PB=json.load(open(S.OUT+'s55a_pb_snr.json'))['pb']
def Rj(pos,J):
    I,e=L.DATA[pos]; j=J-2; r=I[j]/I[3+j]; return r,r*np.hypot(e[j],e[3+j])   # ratio, absolute stat error
out={}
for pos in POS:
    r=[Rj(pos,J) for J in (2,3,4)]; out[pos]=r
    P(f'{pos:7s} R2={r[0][0]:.3f}+-{r[0][1]:.3f}  R3={r[1][0]:.3f}+-{r[1][1]:.3f}  R4={r[2][0]:.3f}+-{r[2][1]:.3f} (abs stat)  R4/R2={r[2][0]/r[0][0]:.2f}  R3/R2={r[1][0]/r[0][0]:.2f}  R4/R3={r[2][0]/r[1][0]:.2f}')
res={'R':{p:[(a,b) for a,b in v] for p,v in out.items()}}
def sig_tot(pos,J,mode):
    R=out[pos][J-2][0]; s_stat=out[pos][J-2][1]; pb=0.3*(1-PB[pos][J-2])*np.sqrt(2)*R if pos in ('E1','E2','W') else 0.0   # absolute errors on R
    if mode=='stat': return s_stat
    if mode=='indep7': return np.sqrt(s_stat**2+2*(0.07*R)**2+pb**2)
    if mode=='calcommon': return np.sqrt(s_stat**2+pb**2)
P('\n# Centre (N) vs disk position: ln[R(disk)/R(N)] with significance, per J; and trend d=ln[R(4-3)/R(2-1)] disk vs N')
for mode in ('stat','indep7','calcommon'):
    P(f'## error model: {mode}')
    for pos in ('S','Bridge','E1','E2','W'):
        zs=[]
        for J in (2,3,4):
            if pos=='W' and J==4: zs.append('  (HCN 4-3 2.7 sigma: limit)'); continue
            dl=np.log(out[pos][J-2][0]/out['N'][J-2][0]); e=np.hypot(sig_tot(pos,J,mode)/out[pos][J-2][0],sig_tot('N',J,mode)/out['N'][J-2][0])
            zs.append(f'J={J}: ratio {np.exp(dl):.2f}  ({dl/e:+.1f} sigma)')
        dd=np.log(out[pos][2][0]/out[pos][0][0])-np.log(out['N'][2][0]/out['N'][0][0]); ed=np.sqrt(sum((sig_tot(p,J,mode)/out[p][J-2][0])**2 for p in (pos,'N') for J in (2,4)))
        P(f'  {pos:6s} '+' | '.join(zs)+f' | trend ln[(R4/R2)_{pos}/(R4/R2)_N] = {dd:+.2f}+-{ed:.2f} ({dd/ed:+.1f} sigma)'); res[f'{mode}|{pos}']=dict(trend=dd,e=ed)
P('\n## Bound allowing HCN-vs-HCO+ differential recovery as large as the Xb4f->Xb4d gap measured at that position (disk R can only be raised: R_true=R_obs/(1-delta), delta=1-rec(HCN)/rec(HCO+) at Xb4d, same for each band; calibration common; centre uncorrected)')
for pos in ('E1','E2','W'):
    dl=max(0.0,1-rec[(pos,'HCN')][0]/rec[(pos,'HCO+')][0]); zs=[]
    for J in (2,3,4):
        if pos=='W' and J==4: continue
        Rc=out[pos][J-2][0]/(1-dl); dd=np.log(Rc/out['N'][J-2][0]); e=np.hypot(sig_tot(pos,J,'calcommon')/out[pos][J-2][0],sig_tot('N',J,'calcommon')/out['N'][J-2][0]); zs.append(f'J={J}: R_disk/R_N={np.exp(dd):.2f} ({dd/e:+.1f} sigma)')
    P(f'  {pos}: delta_max={dl:.2f}; '+' | '.join(zs)); res[f'bound|{pos}']=dl
# differential recovery needed to erase the trend: HCN loses an extra fraction x at Xb47 relative to Xb55 (HCO+ unchanged)
P('\n# Differential HCN-vs-HCO+ recovery (Xb47 vs Xb55) needed to explain the disk R4/R2 deficit vs N: (1-x47)/(1-x55) = (R4/R2)_disk/(R4/R2)_N')
for pos in ('E1','E2'):
    rr=(out[pos][2][0]/out[pos][0][0])/(out['N'][2][0]/out['N'][0][0]); P(f'  {pos}: needs (1-x47)/(1-x55)={rr:.2f}, e.g. HCN extra loss of {100*(1-rr):.0f}% at Xb47 with none at Xb55; measured HCN/HCO+ relative recovery at the 0.83" beam (Xb4d/Xb4f) is {rec[(pos,"HCN")][0]/rec[(pos,"HCO+")][0]:.2f} (a 4.5 -> 0.83" change; Xb47 is 3.4", only {np.log2(4.5/3.4)/np.log2(4.5/0.83):.2f} of that log-baseline step) -> predicted differential at Xb47 ~{100*(1-(rec[(pos,"HCN")][0]/rec[(pos,"HCO+")][0])**(np.log2(4.5/3.4)/np.log2(4.5/0.83))):.0f}%, relative to Xb55 (larger beam) a further ~{100*(1-(rec[(pos,"HCN")][0]/rec[(pos,"HCO+")][0])**(np.log2(4.5/3.4)/np.log2(4.5/0.83)+np.log2(7.3/4.5)/np.log2(4.5/0.83))):.0f}%')
    res[f'needed|{pos}']=rr
# ---------------------------------------------------------------- (3) dilution power test
P('\n# (3) Dilution power test at 10": two-zone injection. Core fluxes (Xb4b r<0.5", |v|<350, Jy km/s) from s43; 10" beam flux = core*1 + envelope (Jy/beam km/s at the position).')
import pandas as pd
b=pd.read_csv('/workspace/out/s43_profiles.csv'); bm=1.1331*0.2456*0.2373/0.047**2
core={z:{m:float((b[b.region==z].reset_index(drop=True)[c]*b[b.region==z].reset_index(drop=True).npix/bm)[:2].sum()) for m,c in (('HCN','I_HCN'),('HCO+','I_HCO'))} for z in ('N','S')}
R30={r['tag']:r for r in json.load(open('/workspace/out/s30_ratios.json'))}
rows=R30['Xb4f_J3_10as']['rows']
P(f'core flux (Jy km/s, 3-2): N HCN {core["N"]["HCN"]:.2f} HCO+ {core["N"]["HCO+"]:.2f} (R={core["N"]["HCN"]/core["N"]["HCO+"]:.2f}); S HCN {core["S"]["HCN"]:.2f} HCO+ {core["S"]["HCO+"]:.2f} (R={core["S"]["HCN"]/core["S"]["HCO+"]:.2f})')
for z in ('N','S','Bridge'):
    r=rows[z+'_pos']; P(f'10" Xb4f {z}: HCN {r["HCN"]:.2f} HCO+ {r["HCOp"]:.2f} Jy/beam km/s -> core fractions {core.get(z,{}).get("HCN",0)/r["HCN"]:.2f} (HCN), {core.get(z,{}).get("HCO+",0)/r["HCOp"]:.2f} (HCO+); R_Jy={r["HCN"]/r["HCOp"]:.3f}')
# null hypothesis: same envelope ratio R_env at N and S (no difference in the gas other than the Q1 cores). Predicted 10" R_N/R_S.
oN,oS=rows['N_pos'],rows['S_pos']; eH=lambda z,r: r['HCN']-core[z]['HCN']; eO=lambda z,r: r['HCOp']-core[z]['HCO+']
RenvN=eH('N',oN)/eO('N',oN); RenvS=eH('S',oS)/eO('S',oS)
P(f'envelope (10" flux minus compact core) HCN/HCO+ (Jy): N {RenvN:.3f}, S {RenvS:.3f}; env ratio N/S = {RenvN/RenvS:.3f}; observed 10" N/S = {(oN["HCN"]/oN["HCOp"])/(oS["HCN"]/oS["HCOp"]):.3f}; observed N/S (total) given next')
Renv=0.5*(RenvN+RenvS); pred=[]
for z,o in (('N',oN),('S',oS)):
    envO=o['HCOp']-core[z]['HCO+']; pred.append((core[z]['HCN']+Renv*envO)/(core[z]['HCO+']+envO))
P(f'   -> injected common-envelope model (R_env={Renv:.3f}) predicts 10" ratios N {pred[0]:.3f}, S {pred[1]:.3f}: N/S = {pred[0]/pred[1]:.3f} (observed {(oN["HCN"]/oN["HCOp"])/(oS["HCN"]/oS["HCOp"]):.3f}+-0.14 incl. 7% per line)')
res['dilution']=dict(core=core,pred_NS=pred[0]/pred[1],obs_NS=(oN["HCN"]/oN["HCOp"])/(oS["HCN"]/oS["HCOp"]),RenvN=RenvN,RenvS=RenvS)
# full ladder injection into the RADEX fit: envelope = best single zone for the mean N/S envelope; add cores with ladders thick/mid/thin; compare dchi2 N vs S ... (pair of ln-flux vectors, error 7%/line + band scales common to N,S)
P('   Power: separation of N and S in the single-zone ladder is quantified by the minimum detectable ratio-of-ratios: 2-sigma threshold = exp(2*0.14)=1.32 (1.15 if the 7% per-line calibration cancels between N and S, since both use the same lines).')
P('   Overlap: aperture-flux correlation between N and S at 10" (Gaussian beam FWHM 10", separation 5.0"): beam overlap = exp(-4 ln2 d^2/(2*FWHM^2))... pixel-noise correlation rho = %.2f; Bridge-N (1.4"): %.2f; Bridge-S (3.6"): %.2f' % tuple(np.exp(-4*np.log(2)*d**2/(2*10.**2)) for d in (5.0,1.4,3.6)))
json.dump(res,open(S.OUT+'s55c_ratios_recovery_power.json','w'),indent=1,default=float); open(S.OUT+'s55c_ratios_recovery_power.txt','w').write('\n'.join(log))
