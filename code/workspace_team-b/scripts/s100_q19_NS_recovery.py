# Q19: is the N>S HCN/HCO+ T43/T32 contrast (p0015; team-a p0022) consistent with 1 after one-sided differential Xb47/Xb4f recovery,
# the Q23 edge-ramp bias and a continuum anchor?  usage: python3 -I s100_q19_NS_recovery.py -> out/q19/s100.txt/json
# Model: observed N/S = L_true * D, D = (R47/R4f)_N/(R47/R4f)_S (D>1: S loses more of the 4-3 flux relative to 3-2).
import json,numpy as np
rng=np.random.default_rng(19); n=200000; OUT='/workspace/out/q19/'; log=[]
def P(s=''): print(s); log.append(s)
# observed (s82, Xb4f native beam, 5.72x3.54): stat+centring (recovery handled separately)
obs={'HCN':(1.262,np.hypot(0.042,0.061)),'HCO+':(1.164,np.hypot(0.016,0.062))}
# other beams (s82): circular 6.5" 1.16/1.10, Xb55 native 1.24/1.13 (HCN/HCO+); team-a p0022 at Xb55 beam 1.257/1.112
ramp={'HCN':(0.98,0.025),'HCO+':(1.00,0.03)}   # Q23: N/S factor on T43/T32 after raising HCN 4-3 by +4+-4% (N), +6+-5% (S); HCO+ 4-3 S blue baseline +-7% -> +-3% on N/S
# continuum Xb47/Xb4f ratio N/S at the matched B5 beam 8.95x5.98 (out/s14_common_beam.csv); rms 1.77/0.87 mJy
c=json.load(open('/dev/stdin')) if False else None
import csv
rows={(r['beam'],r['ous']):r for r in csv.DictReader(open('/workspace/out/s14_common_beam.csv'))}
B='B5beam_8.95x5.98'
def cont(key):
    out={}
    for ap in ('peak','r3','r6','r10'):
        k=lambda s,o: float(rows[(B,o)][f'{s}_peak_mJy' if ap=='peak' else f'{s}_S_{ap}_mJy'])
        r=lambda s: k(s,'Xb47')/k(s,'Xb4f')
        out[ap]=r('N')/r('S')
    return out
C=cont(0); P('# continuum (Xb47/Xb4f)_N/(Xb47/Xb4f)_S at the B5 beam: '+', '.join(f'{k} {v:.3f}' for k,v in C.items()))
for ap in ('peak','r3'):
    rN=[float(rows[(B,o)]['N_peak_mJy' if ap=='peak' else 'N_S_r3_mJy']) for o in ('Xb47','Xb4f')]; rS=[float(rows[(B,o)]['S_peak_mJy' if ap=='peak' else 'S_S_r3_mJy']) for o in ('Xb47','Xb4f')]
# continuum alpha(Xb55-Xb4f) differential, same beam
g=lambda s,o,k='peak': float(rows[(B,o)][f'{s}_peak_mJy'])
dN=g('N','Xb4f')/g('N','Xb55'); dS=g('S','Xb4f')/g('S','Xb55')
P(f'# continuum 4f/55 N/S ratio of ratios {dN/dS:.3f} (peak) -> d_alpha(170-258) = {np.log(dN/dS)/np.log(257.7/170.7):.2f}; 47/4f: d_alpha(258-347)= {np.log(C["peak"])/np.log(346.5/257.7):.2f} (if recovery were identical)')
sig_c=0.05   # stat 3% peaks + pbcor/registration
res={}
def summ(name,x):
    q=np.percentile(x,[2.5,16,50,84,97.5]); pg=(x>1).mean(); res[name]=dict(p=list(q),P_gt1=pg)
    P(f'  {name:62s} median {q[2]:.3f} 68% [{q[1]:.3f},{q[3]:.3f}] 95% [{q[0]:.3f},{q[4]:.3f}]  P(>1)={pg:.2f}  P(1 inside 68%)={int(q[1]<=1<=q[3])}')
for mol in ('HCN','HCO+'):
    L=rng.normal(*obs[mol],n)*rng.normal(*ramp[mol],n)
    P(f'\n## {mol}: observed N/S (Xb4f beam) {obs[mol][0]:.3f}+-{obs[mol][1]:.3f} (stat+centring), x ramp {ramp[mol][0]}+-{ramp[mol][1]}')
    summ(f'{mol} no recovery correction, ramp corrected',L)
    # A: Q12 one-sided priors. R47/R4f at N,S (HCN prior .97/.96, max .92/.89; HCO+ .97/.97, max .93/.91) -> D=R_N/R_S; one-sided D in [1,Dmax]
    Dmax={'HCN':0.92/0.89,'HCO+':0.93/0.91}[mol]
    DA=rng.uniform(1.0,Dmax,n); A=L/DA; summ(f'A Q12 priors one-sided D~U[1,{Dmax:.3f}]',A)
    # B: continuum anchor, equal alpha (D = C/1.345^dalpha), dalpha ~ N(0,sd); C peak 1.137+-0.05 ; also r3
    for ap in ('peak','r3'):
        for sd in (0.0,0.3,0.5):
            Cv=C[ap]*np.exp(rng.normal(0,sig_c,n)); da=rng.normal(0,sd,n); D=Cv/(346.5/257.7)**da
            summ(f'B continuum anchor ({ap}, C={C[ap]:.2f}) d_alpha sd {sd}',L/D)
    # C: team-a physical-Bmin uv model (NS_correction factor_median 0.78 HCN / 0.84 HCO+; p0022 pooled 50/50 values quoted in text)
    Dteam={'HCN':(1/0.782,0.12),'HCO+':(1/0.842,0.08)}[mol]
    D=rng.normal(*Dteam,n); summ(f'C team-a physical uv model D={Dteam[0]:.2f}+-{Dteam[1]}',L/D)
    # pools
    k=n//3; pool2=np.concatenate([A[:n//2], (L/(C['peak']*np.exp(rng.normal(0,sig_c,n))/(346.5/257.7)**rng.normal(0,0.3,n)))[:n//2]])
    summ('pool A+B(peak, sd 0.3) 50/50 (our data only)',pool2)
    pool3=np.concatenate([A[:k],(L/(C['peak']*np.exp(rng.normal(0,sig_c,n))/(346.5/257.7)**rng.normal(0,0.3,n)))[:k],(L/rng.normal(*Dteam,n))[:k]])
    summ('pool A+B+C equal weights',pool3)
json.dump(res,open(OUT+'s100.json','w'),indent=1,default=float); open(OUT+'s100.txt','w').write('\n'.join(log))
