# Q12 step d: recovery table, saturating loss model R(MRS)=1/(1+(s/MRS)^beta) calibrated per position on the nested Xb4b/Xb4d/Xb4f set, extrapolation to the
# Band 5/7 ACA data sets (Xb55, Xb47), ladder correction, PB table.
# Assumptions (all stated in the report): (i) Xb4f, Xb47, Xb55 are ACA (7-m) data: noise-map PB FWHM 37/30/58" (s55a; re-checked here with Xb4d/Xb4b = 12-m);
#   same shortest baseline in metres => MRS_ACA scales with wavelength: MRS_47/MRS_4f = nu4f/nu47 = 0.75, MRS_55/MRS_4f = 1.52 (exponent drawn U[0.5,1.5] for configuration differences).
# (ii) nominal MRS (Band 6, 265 GHz) from Cycle-5 handbook values scaled in frequency: MRS_f~25", MRS_d/MRS_f = 0.30 (ln-uniform 0.20-0.45), MRS_b/MRS_d = 0.40 (0.30-0.55);
#   only RATIOS enter, and the data test them (alpha_bd vs alpha_df consistency, printed).
# (iii) true flux F = F_obs/R, R<=1 for every array (loss-only). Per position the three measured ratios r_d=F_d/F_f, r_b=F_b/F_f (10" or 6.5" common beam) fix (s,beta,R_f).
# (iv) 5% independent flux-scale error per array; line-free-channel noise from s56b. J=2-1/4-3 sky scale s_J = s*kJ, kJ~U[1,1.6] (2-1), U[0.6,1] (4-3) (more/less extended than 3-2).
# usage: python3 -I s56d_recovery_model.py -> out/q12/recovery_model.json/.txt, recovery_table.csv, pb_table.txt
import sys,json,csv; sys.path.insert(0,'/workspace/scripts')
import numpy as np, s56_lib as L
from scipy.optimize import brentq
from regions import REGIONS
rng=np.random.default_rng(12345); NMC=2000
log=[]
def P(s=''): print(s,flush=True); log.append(s)
M={}; 
for x in csv.DictReader(open('/workspace/out/q12/measure.csv')): M[(x['line'],x['ous'],x['beam'],x['pos'])]=(float(x['val']),float(x['err']))
G={}
for x in csv.DictReader(open('/workspace/out/q12/growth.csv')): G[(x['line'],x['ous'],x['beam'],'ap%s'%x['R'])]=(float(x['val']),float(x['err']))
M.update(G)
POS=['N','S','Bridge','E1','E2','W','Sblue']
# ---------- 1. empirical table
rows=[]
P('## 1. Empirical recovery (value ratio) at common circular beams; err = noise (16 sign-flip realisations) quad. 5% per-array scale')
P('line pos      beam | Xb4b/Xb4d | Xb4d/Xb4f | Xb4b/Xb4f')
def ratio(line,pos,beam,a,b):
    if (line,a,beam,pos) not in M or (line,b,beam,pos) not in M: return None
    (va,ea),(vb,eb)=M[(line,a,beam,pos)],M[(line,b,beam,pos)]
    if vb<=0: return None
    r=va/vb; e=abs(r)*np.hypot(np.hypot(ea/max(abs(va),1e-9),eb/vb),0.07); return r,e
for line in ('HCO+','HCN'):
    for pos in POS+['ap4','ap8','ap12']:
        for beam in ('1.05','2.0','3.0','4.5','6.5','10','xb4f'):
            out=[ratio(line,pos,beam,'Xb4b','Xb4d'),ratio(line,pos,beam,'Xb4d','Xb4f'),ratio(line,pos,beam,'Xb4b','Xb4f')]
            if all(o is None for o in out): continue
            f=lambda o:'   --      ' if o is None else f'{o[0]:6.2f}+-{o[1]:.2f}'
            P(f'{line:5s}{pos:8s}{beam:6s}| {f(out[0])} | {f(out[1])} | {f(out[2])}')
            rows.append(dict(line=line,pos=pos,beam=beam,**{k:(None if o is None else round(o[0],4)) for k,o in zip(('bd','df','bf'),out)},**{'e'+k:(None if o is None else round(o[1],4)) for k,o in zip(('bd','df','bf'),out)}))
with open('/workspace/out/q12/recovery_table.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
# ---------- 2. saturating model
def solve(rd,rb,md,mb):
    """(1+x)/(1+x md^b)=rd, (1+x)/(1+x mb^b)=rb  ->  (x_f,beta) or None"""
    def x_of(b):
        den=rd*md**b-1
        return (1-rd)/den if den>0 else np.inf
    def g(b):
        x=x_of(b)
        return np.nan if not np.isfinite(x) else (1+x)/(1+x*mb**b)-rb
    bs=np.exp(np.linspace(np.log(0.05),np.log(10),120)); gs=[g(b) for b in bs]
    for i in range(len(bs)-1):
        if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and gs[i]*gs[i+1]<0:
            b=brentq(g,bs[i],bs[i+1]); return x_of(b),b
    return None
def Rmod(x,beta,ratio): return 1/(1+x*ratio**beta)       # ratio = MRS_f/MRS_a
def model(line,key,beam,cls_note=''):
    (vf,ef)=M[(line,'Xb4f',beam,key)]; (vd,ed)=M[(line,'Xb4d',beam,key)]; (vb,eb)=M[(line,'Xb4b',beam,key)]
    out=[]; nfail=0
    for _ in range(NMC):
        gf,gd,gb=1+0.05*rng.standard_normal(3)
        f=vf*gf+ef*rng.standard_normal(); d=vd*gd+ed*rng.standard_normal(); b=vb*gb+eb*rng.standard_normal()
        if f<=0 or d<=0 or b<=0: nfail+=1; continue
        md=1/np.exp(rng.uniform(np.log(0.20),np.log(0.45))); mb=md/np.exp(rng.uniform(np.log(0.30),np.log(0.55)))
        s=solve(d/f,b/f,md,mb)
        if s is None: nfail+=1; continue
        x,be=s; out.append((x,be,md,mb))
    return np.array(out),nfail
# data-only slope consistency (central MRS ratios) -------------
P('\n## 2a. effective power-law slopes alpha (F ~ MRS^alpha) with central MRS ratios, B=10" and 6.5" (consistency test of the nominal MRS ratios)')
for line in ('HCO+','HCN'):
    for key in ('N','S','Bridge','Sblue','ap4','ap8'):
        s=f'{line:5s}{key:7s}'
        for B in ('6.5','10'):
            rb=ratio(line,key,B,'Xb4b','Xb4d'); rd=ratio(line,key,B,'Xb4d','Xb4f')
            if rb and rd and rb[0]>0 and rd[0]>0: s+=f' B={B}: alpha_bd={np.log(rb[0])/np.log(0.40):5.2f} alpha_df={np.log(rd[0])/np.log(0.30):5.2f} |'
        P(s)
RES={}
P('\n## 2b. saturating-loss model R(MRS)=1/(1+(s/MRS)^beta): per-position solution at the common beam (median [16,84] over MC; frac = fraction of MC draws with a solution)')
CLASS={'N':'centre','S':'centre','Bridge':'centre','Sblue':'centre','E1':'disk','E2':'disk','W':'disk','ap4':'centre','ap8':'mix','ap12':'mix'}
for line in ('HCO+','HCN'):
    for key in ('N','S','Bridge','Sblue','ap4','ap8','E1','E2','W'):
        for B in ('6.5','10'):
            try: out,nf=model(line,key,B)
            except KeyError: continue
            frac=len(out)/NMC
            if len(out)<100: P(f'{line:5s}{key:7s}B={B:4s}: no solution in {1-frac:.0%} of draws (Xb4b flux <=0 or non-monotonic -> loss not a smooth power-law-like function of MRS; flux at this position is dominated by the neighbours\' bowl)'); RES[f'{line}|{key}|{B}']=dict(frac=frac); continue
            x,be,md,mb=out.T; Rf=1/(1+x); Rd=1/(1+x*md**be); Rb=1/(1+x*mb**be)
            RES[f'{line}|{key}|{B}']=dict(frac=frac,Rf=np.percentile(Rf,[16,50,84]).tolist(),Rd=np.percentile(Rd,[16,50,84]).tolist(),Rb=np.percentile(Rb,[16,50,84]).tolist(),beta=np.percentile(be,[16,50,84]).tolist(),
                                          s_arcsec=np.percentile(25*x**(1/be),[16,50,84]).tolist())
            r=RES[f'{line}|{key}|{B}']
            P(f'{line:5s}{key:7s}B={B:4s} frac={frac:.2f}  R(Xb4f)={r["Rf"][1]:.2f}[{r["Rf"][0]:.2f},{r["Rf"][2]:.2f}]  R(Xb4d)={r["Rd"][1]:.2f}[{r["Rd"][0]:.2f},{r["Rd"][2]:.2f}]  R(Xb4b)={r["Rb"][1]:.2f}[{r["Rb"][0]:.2f},{r["Rb"][2]:.2f}]  beta={r["beta"][1]:.2f}  s={r["s_arcsec"][1]:.1f}"')
json.dump(RES,open('/workspace/out/q12/recovery_model.json','w'),indent=1)
open('/workspace/out/q12/recovery_model.txt','w').write('\n'.join(log))
