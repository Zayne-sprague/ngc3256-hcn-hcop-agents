# Q12 step e: loss-only extrapolation of the nested-set recovery to the Band 7 (Xb47) and Band 5 (Xb55) ACA data, PB table, ladder correction.
# Model (stated assumptions, see s56d header for the MRS ratios): within a data-calibrated range the retained flux is a power law in the maximum recoverable scale,
#   F_a ~ MRS_a^alpha  (alpha from Xb4d/Xb4f at 10", MRS_d/MRS_f~0.30 [0.20,0.45]; centre and disk positions separately). Slopes of Xb4b/Xb4d and Xb4d/Xb4f agree at N,S,Bridge (s56d 2a),
#   so a single power law describes MRS = 3-25" there; at the clumps it steepens below MRS~7" (Xb4b keeps ~0-30%), so alpha_df is the local slope nearest the ACA scale.
#   Real recovery curves saturate (alpha falls to 0 for MRS >> structure), hence a saturation factor sat~U[0,1] multiplies alpha (loss-only: R<=1; R_55>=R_4f>=R_47 follows from MRS ~ wavelength).
#   MRS_a/MRS_4f = (nu_4f/nu_a)^e, e~U[0.5,1.5] (same shortest baseline in metres; different configurations/dates give e<>1).  J-dependence of the sky structure: alpha_J = alpha*kJ,
#   kJ~U[1,1.5] for 2-1 (more extended), U[0.5,1] for 4-3 (more compact).  Absolute Xb4f recovery (NOT measurable from the nested set): R_4f=(25"/MRS_out)^alpha, MRS_out ln-U[40",100"].
# usage: python3 -I s56e_extrapolate.py -> out/q12/extrapolation.json/.txt, out/q12/pb_table.txt
import sys,json,csv; sys.path.insert(0,'/workspace/scripts')
import numpy as np
from scipy.optimize import least_squares
from regions import REGIONS
import s56_lib as L
rng=np.random.default_rng(777); N=4000; log=[]
def P(s=''): print(s,flush=True); log.append(s)
M={}
for x in csv.DictReader(open('/workspace/out/q12/measure.csv')): M[(x['line'],x['ous'],x['beam'],x['pos'])]=(float(x['val']),float(x['err']))
POS=['N','S','Bridge','E1','E2','W']
NU={'Xb55':177.26e9*0+175.5e9,'Xb4f':265.7e9,'Xb47':351.5e9}     # centre of the observed bands, only ratios enter
# ---- alpha_df per line/position (10" and 6.5")
def alpha_draws(line,pos,beam):
    (d,ed),(f,ef)=M[(line,'Xb4d',beam,pos)],M[(line,'Xb4d' if False else 'Xb4f',beam,pos)]
    md=1/np.exp(rng.uniform(np.log(0.20),np.log(0.45),N)); rd=(d+ed*rng.standard_normal(N))/(f+ef*rng.standard_normal(N))*(1+0.07*rng.standard_normal(N))
    ok=rd>0.02; a=np.where(ok,np.log(np.clip(rd,1e-3,None))/-np.log(md),np.nan); return a
P('## 1. alpha_df (retained flux ~ MRS^alpha), Xb4d/Xb4f at 10" and 6.5"; median [16,84] incl. 7% scale error and MRS-ratio range')
AL={}
for line in ('HCO+','HCN'):
    for pos in POS:
        for B in ('10','6.5'):
            a=alpha_draws(line,pos,B); AL[(line,pos,B)]=a; v=np.nanpercentile(a,[16,50,84]); P(f'{line:5s}{pos:7s}B={B:4s} alpha_df={v[1]:5.2f} [{v[0]:.2f},{v[2]:.2f}]  (Xb4d/Xb4f={M[(line,"Xb4d",B,pos)][0]/M[(line,"Xb4f",B,pos)][0]:.2f})')
# ---- multipliers
def mult(a,J):
    n=len(a); e=rng.uniform(0.5,1.5,n); sat=rng.uniform(0,1,n); k={2:rng.uniform(1,1.5,n),3:np.ones(n),4:rng.uniform(0.5,1,n)}[J]
    nu_a={2:NU['Xb55'],4:NU['Xb47']}[J]; ratio=(NU['Xb4f']/nu_a)**e
    return ratio**(np.clip(a,0,None)*k*sat)          # R_a / R_4f
P('\n## 2. relative multipliers R_a/R_4f (loss-only, a=Xb55 for 2-1, Xb47 for 4-3) at 10", with the Xb4f absolute recovery R_4f (range from outer scale)')
EX={}
for line in ('HCO+','HCN'):
    for pos in POS:
        a=AL[(line,pos,'10')]; a=np.where(np.isfinite(a),a,np.nanmedian(a))
        m55=mult(a,2); m47=mult(a,4); Rf=np.minimum(1,(25/np.exp(rng.uniform(np.log(40),np.log(100),N)))**np.clip(a,0,None))
        EX[(line,pos)]=dict(m55=m55,m47=m47,Rf=Rf,a=a)
        q=lambda x:np.percentile(x,[16,50,84])
        P(f'{line:5s}{pos:7s} R55/R4f={q(m55)[1]:.2f}[{q(m55)[0]:.2f},{q(m55)[2]:.2f}]  R47/R4f={q(m47)[1]:.2f}[{q(m47)[0]:.2f},{q(m47)[2]:.2f}]  R4f(abs)={q(Rf)[1]:.2f}[{q(Rf)[0]:.2f},{q(Rf)[2]:.2f}]  R55={q(Rf*m55)[1]:.2f}  R47={q(Rf*m47)[1]:.2f}')
# ---- ladder correction using the s30 10" fluxes
S30={x['tag']:x['rows'] for x in json.load(open('/workspace/out/s30_ratios.json'))}
tags={2:'Xb55_J2_10as',3:'Xb4f_J3_10as',4:'Xb47_J4_10as'}
F0={'HCN':{2:177.26111,3:265.88618,4:354.50547},'HCO+':{2:178.37506,3:267.55763,4:356.73424}}
key={'HCN':'HCN','HCO+':'HCOp'}
P('\n## 3. brightness-temperature ladder ratios (T ~ S/nu^2 at the common 10" beam): observed vs recovery-corrected (median [16,84]); calibration 10%/band added')
LAD={}; DR={}
for line in ('HCN','HCO+'):
    for pos in POS:
        S={J:S30[tags[J]][pos+'_pos'][key[line]] for J in (2,3,4)}; E={J:S30[tags[J]][pos+'_pos']['e'+key[line]] for J in (2,3,4)}
        T={J:S[J]/F0[line][J]**2 for J in S}
        r21=T[2]/T[3]; r43=T[4]/T[3]; ex=EX[(line,pos)]
        cal=np.exp(0.10*rng.standard_normal((3,N)))
        c21=r21*cal[0]/cal[1]/ex['m55']; c43=r43*cal[2]/cal[1]/ex['m47']       # true = S/R: T2/T3 -> *(R_4f/R_55)
        o21=r21*cal[0]/cal[1]; o43=r43*cal[2]/cal[1]; q=lambda x:np.percentile(x,[16,50,84])
        LAD[(line,pos)]=(q(o21),q(c21),q(o43),q(c43)); DR[(line,pos)]=(o21,c21,o43,c43)
        P(f'{line:5s}{pos:7s} T21/T32 obs {r21:5.2f} (+cal {q(o21)[0]:.2f}-{q(o21)[2]:.2f}) -> corr {q(c21)[1]:.2f}[{q(c21)[0]:.2f},{q(c21)[2]:.2f}]   T43/T32 obs {r43:5.2f} -> corr {q(c43)[1]:.2f}[{q(c43)[0]:.2f},{q(c43)[2]:.2f}]   (S/N 4-3 {S[4]/E[4]:.0f})')

# ---- disk-vs-centre contrast of the ladder ratios and a stress test
P('\n## 6. disk/centre contrast of T21/T32 and T43/T32 (disk position / N), observed vs corrected (same calibration draws NOT shared between positions; recovery draws independent)')
for line in ('HCN','HCO+'):
    for pos in ('E1','E2','W'):
        o21,c21,o43,c43=DR[(line,pos)]; n21,nc21,n43,nc43=DR[(line,'N')]
        q=lambda x:np.percentile(x,[16,50,84])
        P(f'{line:5s}{pos:3s}/N  T21/T32: obs {q(o21/n21)[1]:.2f}[{q(o21/n21)[0]:.2f},{q(o21/n21)[2]:.2f}] corr {q(c21/nc21)[1]:.2f}[{q(c21/nc21)[0]:.2f},{q(c21/nc21)[2]:.2f}]   T43/T32: obs {q(o43/n43)[1]:.2f}[{q(o43/n43)[0]:.2f},{q(o43/n43)[2]:.2f}] corr {q(c43/nc43)[1]:.2f}[{q(c43/nc43)[0]:.2f},{q(c43/nc43)[2]:.2f}]')
P('\n## 7. stress test (extreme corner: alpha at its 84th percentile, saturation 1, e=1.5, kJ=3 for 2-1 and 1.5 for 4-3): minimum corrected T21/T32 and maximum corrected T43/T32 (central observed ratio, no calibration error)')
for line in ('HCN','HCO+'):
    for pos in POS:
        a=np.nanpercentile(AL[(line,pos,'10')],84); S={J:S30[tags[J]][pos+'_pos'][key[line]] for J in (2,3,4)}; T={J:S[J]/F0[line][J]**2 for J in S}
        m55=(NU['Xb4f']/NU['Xb55'])**(1.5*a*3); m47=(NU['Xb4f']/NU['Xb47'])**(1.5*a*1.5)
        P(f'{line:5s}{pos:7s} alpha84={a:.2f}  R55/R4f<={m55:.2f}  R47/R4f>={m47:.2f}   T21/T32 {T[2]/T[3]:.2f} -> >= {T[2]/T[3]/m55:.2f}   T43/T32 {T[4]/T[3]:.2f} -> <= {T[4]/T[3]/m47:.2f}')
# ---- differential HCN vs HCO+ loss for same-J ratios
P('\n## 4. same-J HCN/HCO+ ratio: factor by which the ACA-observed ratio could be low/high if HCN and HCO+ have different slopes (true/obs = (MRS_out/MRS_a)^(alpha_HCN-alpha_HCO+); MRS_out ln-U[40,100]", MRS_4f=25"):')
for pos in POS:
    d=EX[('HCN',pos)]['a']-EX[('HCO+',pos)]['a']; mo=np.exp(rng.uniform(np.log(40),np.log(100),N)); f=(mo/25)**d
    q=np.percentile(f,[16,50,84]); P(f'{pos:7s} delta_alpha={np.median(d):5.2f} [{np.percentile(d,16):.2f},{np.percentile(d,84):.2f}]  true/obs HCN/HCO+ (Xb4f 3-2) = {q[1]:.2f}[{q[0]:.2f},{q[2]:.2f}]')
# ---- PB
P('\n## 5. Primary beam: which dish? noise-map PB fits (pbcor noise ~ 1/PB), 12-m check with Xb4d/Xb4b and ACA candidates')
def fit_pb(ous,line,beam=None):
    d=np.load(f'/workspace/out/q12/maps_{ous}_{line}.npz'); e=d['noise'].std(0); ps=float(d['ps']); n=e.shape[0]; Y,X=np.mgrid[:n,:n]
    ok=np.isfinite(e)&(e>0)
    def f(p): a,x0,y0,b=p; return (a+b*((X-x0)**2+(Y-y0)**2)*ps**2)[ok]-np.log(e[ok])
    s=least_squares(f,[np.log(np.median(e[ok])),n/2,n/2,1e-3]); b=s.x[3]; return np.sqrt(4*np.log(2)/b) if b>0 else np.nan,np.sqrt(np.mean(s.fun**2))
lam={'Xb4f':2.998e8/265.7e9,'Xb47':2.998e8/351.5e9,'Xb55':2.998e8/175.5e9}
for ous in ('Xb4b','Xb4d','Xb4f'):
    fw,rms=fit_pb(ous,'HCO+'); P(f'{ous} (Band 6) crop noise-map FWHM_PB = {fw:.1f}" (rms resid {rms:.2f}); expected 1.13 lambda/D: 12 m {1.13*lam["Xb4f"]/12*206265:.1f}", 7 m {1.13*lam["Xb4f"]/7*206265:.1f}"')
for ous in ('Xb55','Xb4f','Xb47'):
    P(f'{ous}: 12-m expectation {1.13*lam[ous]/12*206265:.1f}", 7-m expectation {1.13*lam[ous]/7*206265:.1f}"')
pb=json.load(open('/workspace/out/q3b/s55a_pb_snr.json'))
FW={'Xb55':57.5,'Xb4f':37.2,'Xb47':30.2}   # fitted PB FWHM from s55a (mean of HCN/HCO+: 58.4/65.9, 37.2/37.4, 29.6/30.8)
FIELD=(156.96346,-43.90461)
P('\npos      r_field["]  PB(Xb55)  PB(Xb4f)  PB(Xb47)  | PB(Xb47) if 12-m (17") | pbcor correction factor 1/PB: Xb47 fitted / 12-m')
for pos in POS:
    ra,dec=REGIONS[pos][:2]; r=np.hypot((ra-FIELD[0])*np.cos(np.deg2rad(FIELD[1])),dec-FIELD[1])*3600
    g=lambda fw:np.exp(-4*np.log(2)*r**2/fw**2); P(f'{pos:7s} {r:6.1f}      {g(FW["Xb55"]):.2f}      {g(FW["Xb4f"]):.2f}      {g(FW["Xb47"]):.2f}     |   {g(16.6):.2f}   |   {1/g(FW["Xb47"]):.2f} / {1/g(16.6):.2f}')
P('PB-model sensitivity: a 10% error in the FWHM changes the pbcor factor at E1 (r~8") of Xb47 by '+f'{abs(np.exp(4*np.log(2)*(7.8**2)*(1/30.2**2-1/(1.1*30.2)**2))-1)*100:.0f}% (ACA) vs 12-m case would be a {abs(np.exp(4*np.log(2)*(7.8**2)*(1/16.6**2-1/(1.1*16.6)**2))-1)*100:.0f}% effect')
out={f'{k[0]}|{k[1]}':[[float(x) for x in v] for v in vals] for k,vals in LAD.items()}
json.dump(out,open('/workspace/out/q12/extrapolation.json','w'),indent=1); open('/workspace/out/q12/extrapolation.txt','w').write('\n'.join(log))
