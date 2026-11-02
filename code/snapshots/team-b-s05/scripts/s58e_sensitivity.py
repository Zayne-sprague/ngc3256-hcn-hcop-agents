# Q12 round 2, step e: the ACA (Xb55, Xb47) vs Xb4f relative-recovery extrapolation, restated as CONDITIONAL on its priors, with posterior sensitivity.
# Model (unchanged from s56e): retained flux R ~ MRS^alpha_loc (alpha from the Xb4d/Xb4f ratio, MRS_d/MRS_f=0.30 [0.20,0.45] ln-U, 7% scale error); for another ACA set
#   R_a/R_4f = (nu_4f/nu_a)^(e*alpha*kJ*sat)  (loss-only: exponent>=0).  Priors varied: sat (U[0,1] | 1 | 0.5), kJ (J-dependent structure: U[1,1.5] for 2-1, U[0.5,1] for 4-3 | 1), e (U[0.5,1.5] | 1+-0.1).
# Inputs: centre positions N,S,Bridge: 10" point values (measure.csv, s30 ladder fluxes); disk E1,E2,W: de-contaminated clump total fluxes (s58d, V0 +- stat, sys) so
#   the alpha and ladder ratios at the clumps no longer contain the nuclear bowl.
# Ladder: Tb ratio = (S/nu^2) ratio; corrected = observed / (R_a/R_4f); 10% calibration per band added.
# usage: python3 -I s58e_sensitivity.py -> out/q12b/sensitivity.txt/.json
import sys,json,csv; sys.path.insert(0,'/workspace/scripts')
import numpy as np
rng=np.random.default_rng(58); N=6000; log=[]
def P(s=''): print(s,flush=True); log.append(s)
q=lambda x:np.nanpercentile(x,[16,50,84])
M={}
for x in csv.DictReader(open('/workspace/out/q12/measure.csv')): M[(x['line'],x['ous'],x['beam'],x['pos'])]=(float(x['val']),float(x['err']))
CF=json.load(open('/workspace/out/q12b/clumpfit.json'))
S30={x['tag']:x['rows'] for x in json.load(open('/workspace/out/s30_ratios.json'))}
tags={2:'Xb55_J2_10as',3:'Xb4f_J3_10as',4:'Xb47_J4_10as'}
F0={'HCN':{2:177.26111,3:265.88618,4:354.50547},'HCO+':{2:178.37506,3:267.55763,4:356.73424}}; key={'HCN':'HCN','HCO+':'HCOp'}
NU={'Xb55':178.0,'Xb4f':266.7,'Xb47':355.6}   # mean of HCN/HCO+ line frequencies [GHz]; only ratios enter
def cf(line,ous,c):
    d=CF[f'{line}|{ous}']; v=d['V0'][c]; e=np.hypot(d['V0']['err'][c],np.std([d[k][c] for k in d])); return v,e
def ratio_draws(line,pos):
    """Xb4d/Xb4f draws: centre -> point values at 10"; clumps -> fitted total clump flux"""
    if pos in ('N','S','Bridge'):
        (d,ed),(f,ef)=M[(line,'Xb4d','10',pos)],M[(line,'Xb4f','10',pos)]
    else:
        (d,ed),(f,ef)=cf(line,'Xb4d',pos),cf(line,'Xb4f',pos)
    return (d+ed*rng.standard_normal(N))/(f+ef*rng.standard_normal(N))*(1+0.07*rng.standard_normal(N)),d/f if f>0 else np.nan
def alpha(r):
    md=1/np.exp(rng.uniform(np.log(0.20),np.log(0.45),N)); ok=r>0.02; return np.where(ok,np.log(np.clip(r,1e-3,None))/-np.log(md),np.nan)
SC={'prior (sat~U[0,1], kJ~U, e~U[0.5,1.5])':dict(sat='U',kJ='U',e='U'),'sat=1':dict(sat=1,kJ='U',e='U'),'sat=0.5':dict(sat=0.5,kJ='U',e='U'),'kJ=1':dict(sat='U',kJ=1,e='U'),
    'e=1+-0.1':dict(sat='U',kJ='U',e='N'),'e=1+-0.1,sat=1,kJ=1 (max)':dict(sat=1,kJ=1,e='N'),'e=1+-0.1,sat=0.5,kJ=1':dict(sat=0.5,kJ=1,e='N')}
def mult(a,J,sc):
    n=len(a); e={'U':lambda:rng.uniform(0.5,1.5,n),'N':lambda:1+0.1*rng.standard_normal(n)}[sc['e']]()
    sat=rng.uniform(0,1,n) if sc['sat']=='U' else sc['sat']*np.ones(n)
    k=(sc['kJ']*np.ones(n)) if sc['kJ']!='U' else {2:rng.uniform(1,1.5,n),4:rng.uniform(0.5,1,n)}[J]
    nu=NU['Xb55'] if J==2 else NU['Xb47']; return (NU['Xb4f']/nu)**(e*np.clip(a,0,None)*k*sat)
POS=['N','S','Bridge','E1','E2','W']; RES={}
P('## alpha_df = d ln(retained flux)/d ln MRS from Xb4d/Xb4f (10" point values at N,S,Bridge; clump total fluxes from the s58d forward model at E1,E2,W)')
AL={}
for line in ('HCO+','HCN'):
    for pos in POS:
        r,rc=ratio_draws(line,pos); a=alpha(r); AL[(line,pos)]=a; v=q(a)
        P(f'{line:5s}{pos:7s} Xb4d/Xb4f = {np.nanmedian(r):.2f} [{q(r)[0]:.2f},{q(r)[2]:.2f}] -> alpha = {v[1]:.2f} [{v[0]:.2f},{v[2]:.2f}]  (fraction of draws with defined alpha {np.mean(np.isfinite(a)):.2f})')
P('\n## R_a/R_4f by scenario (conditional on the priors); median [16,84]')
for sn,sc in SC.items():
    P(f'\n### {sn}')
    for line in ('HCO+','HCN'):
        for pos in POS:
            a=AL[(line,pos)]; a=np.where(np.isfinite(a),a,np.nanmedian(a))
            if np.all(~np.isfinite(a)): continue
            m55=mult(a,2,sc); m47=mult(a,4,sc); RES[f'{sn}|{line}|{pos}']=dict(m55=q(m55).tolist(),m47=q(m47).tolist())
            P(f'{line:5s}{pos:7s} R55/R4f={q(m55)[1]:.2f}[{q(m55)[0]:.2f},{q(m55)[2]:.2f}]  R47/R4f={q(m47)[1]:.2f}[{q(m47)[0]:.2f},{q(m47)[2]:.2f}]')
# ---- ladder ratios: centre from the 10" values, disk from fitted clump fluxes
def ladder_obs(line,pos):
    if pos in ('N','S','Bridge'):
        S={J:S30[tags[J]][pos+'_pos'][key[line]] for J in (2,3,4)}; E={J:S30[tags[J]][pos+'_pos']['e'+key[line]] for J in (2,3,4)}
    else:
        ous={2:'Xb55',3:'Xb4f',4:'Xb47'}; S={};E={}
        for J in (2,3,4): S[J],E[J]=cf(line,ous[J],pos)
    return S,E
P('\n## Tb ladder ratios (S/nu^2): observed vs corrected for the relative recovery, per scenario; 10% calibration per band included. Disk positions use clump total fluxes (nuclear contamination removed).')
LAD={}
for line in ('HCN','HCO+'):
    for pos in POS:
        S,E=ladder_obs(line,pos); T={J:S[J]/F0[line][J]**2 for J in S}; cal=np.exp(0.10*rng.standard_normal((3,N)))
        if not (S[3]>3*E[3] and S[2]>3*E[2]): P(f'{line:5s}{pos:7s} 3-2 or 2-1 flux below 3 sigma in the clump fit (2-1 {S[2]:.2f}+-{E[2]:.2f}, 3-2 {S[3]:.2f}+-{E[3]:.2f}): no ladder ratio'); continue
        d=lambda J:1+E[J]/S[J]*rng.standard_normal(N)
        t21=T[2]*d(2)/(T[3]*d(3))*cal[0]/cal[1]; o21=T[2]/T[3]; det4=S[4]>3*E[4]
        t43=(T[4]*d(4)/(T[3]*d(3))*cal[2]/cal[1]) if det4 else np.full(N,3*E[4]/F0[line][4]**2/T[3]); o43=T[4]/T[3]
        row=f'{line:5s}{pos:7s} obs T21/T32={o21:5.2f} [{q(t21)[0]:.2f},{q(t21)[2]:.2f}] '+(f'T43/T32={o43:5.2f} [{q(t43)[0]:.2f},{q(t43)[2]:.2f}]' if det4 else f'T43/T32 < {t43[0]:.2f} (3 sigma; 4-3 {S[4]:.2f}+-{E[4]:.2f})')+' ->'
        a=AL[(line,pos)]; a=np.where(np.isfinite(a),a,np.nanmedian(a)); P(row)
        for sn in ('prior (sat~U[0,1], kJ~U, e~U[0.5,1.5])','sat=0.5','e=1+-0.1,sat=1,kJ=1 (max)'):
            sc=SC[sn]; m55=mult(a,2,sc); m47=mult(a,4,sc); c21=t21/m55; c43=t43/m47; LAD[f'{line}|{pos}|{sn}']=(q(c21).tolist(),q(c43).tolist())
            P(f'      {sn:42s} T21/T32 -> {q(c21)[1]:.2f}[{q(c21)[0]:.2f},{q(c21)[2]:.2f}]   T43/T32 -> '+(f'{q(c43)[1]:.2f}[{q(c43)[0]:.2f},{q(c43)[2]:.2f}]' if det4 else f'< {q(c43)[1]:.2f}'))
P('\n## disk/centre (E1,E2,W over the centre point N and over the nuclear-complex total) of T21/T32 and T43/T32, observed (calibration errors independent), HCO+ and HCN')
for line in ('HCN','HCO+'):
    Sn,En=ladder_obs(line,'N'); Tn={J:Sn[J]/F0[line][J]**2 for J in Sn}
    for pos in ('E1','E2','W'):
        S,E=ladder_obs(line,pos); T={J:S[J]/F0[line][J]**2 for J in S}
        if not (S[3]>3*E[3] and S[2]>3*E[2]): continue
        P(f'{line:5s}{pos}: T21/T32 disk {T[2]/T[3]:.2f}  N {Tn[2]/Tn[3]:.2f}  ratio {T[2]/T[3]/(Tn[2]/Tn[3]):.2f}; T43/T32 disk {T[4]/T[3]:.3f} N {Tn[4]/Tn[3]:.2f} ratio {T[4]/T[3]/(Tn[4]/Tn[3]):.2f}'+('' if S[4]>3*E[4] else ' (4-3 not detected)'))
# nuclear-complex totals from the s58d fits (all nuclear components, beam-matched by the forward model)
P('\n## nuclear-complex (N+S+Bridge+haloA+haloB) total fluxes from the s58d fits and the same ladder ratios; clumps vs complex')
for line in ('HCN','HCO+'):
    ous={2:'Xb55',3:'Xb4f',4:'Xb47'}; S={J:CF[f'{line}|{ous[J]}']['V0']['nuc_total'] for J in (2,3,4)}; T={J:S[J]/F0[line][J]**2 for J in S}
    P(f'{line:5s}nuclear complex fluxes {S[2]:.1f}/{S[3]:.1f}/{S[4]:.1f} Jy km/s (2-1/3-2/4-3): T21/T32 {T[2]/T[3]:.2f}, T43/T32 {T[4]/T[3]:.2f}')
json.dump(dict(res=RES,lad=LAD),open('/workspace/out/q12b/sensitivity.json','w')); open('/workspace/out/q12b/sensitivity.txt','w').write('\n'.join(log))
