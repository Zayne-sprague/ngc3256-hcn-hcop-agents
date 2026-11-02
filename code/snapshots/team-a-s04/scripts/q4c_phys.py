# Q4c-5: physical consequences with a single-cloud line width. (1) RADEX depends on N/dv: verify with direct dv=50 runs against dv=150 at 3x column. (2) Posterior (physT30, 'all' data, classes A, B, F2) of N(HCN) per cloud
# (dv=50 km/s), number of clouds along the line of sight to build the observed ~150 km/s line, N(H2) and Sigma for X(HCN)=1e-8 (1e-9..1e-7), RADEX n and Tkin -> P/k, path length n -> volume filling.
# Comparison with CO-based Sigma 8-5500 Msun/pc2 (Brunetti+2021) and cloud turbulent pressure median 1e7 K cm-3 (Brunetti+2022). Output out/q4c/phys_*.csv
import sys,subprocess,os; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
import q4c_joint as J, q4c_grid as GR
O=J.O
# (1) spot checks
rng=np.random.default_rng(3); rows=[]
def run_dv(mol,T,ln,lN,dv):
    GR.DV=dv; r=GR.run((mol,T,ln,lN)); return r
for i in range(16):
    mol=['hcn','hco+'][i%2]; T=float(rng.choice([30,40,65,100])); ln=float(rng.choice(np.arange(4,6.3,0.25))); lN=float(rng.choice(np.arange(15,17.6,0.25)))
    a=run_dv(mol,T,ln,lN,150.); b=run_dv(mol,T,ln,lN-np.log10(3),50.)
    a=np.array(a,float); b=np.array(b,float)   # rows: lines (1-0..4-3), cols TR,tau,Tex
    rows.append(dict(mol=mol,T=T,logn=ln,logN150=lN,max_rel_dTR=np.nanmax(np.abs(b[:,0]/a[:,0]-1)),max_rel_dtau=np.nanmax(np.abs(b[:,1]/a[:,1]-1)),max_dTex_K=np.nanmax(np.abs(b[:,2]-a[:,2]))))
chk=pd.DataFrame(rows); chk.to_csv(O+'phys_dv_check.csv',index=False); print(chk.round(4).to_string()); GR.DV=150.
# (2) posteriors
PC=3.0857e18; MSUN_PC2_PER_CM2=1.36*2*1.6726e-24*(3.0857e18**2)/1.989e33   # N(H2) -> Sigma incl. He
def q(x,w,p):
    i=np.argsort(x); c=np.cumsum(w[i]); c/=c[-1]; return np.interp(p,c,x[i])
out=[]
for c in ('S','N'):
    for cl in ('A both thick, kappa=1','B both thick, same Tex (0.8-1.25), kappa free','F2 separate n, kappa=1, both thick'):
        P=J.posterior(c,'all',[cl]); w=P['w']; w=w/w.sum()
        nH2=10**P['logn_H']; T=P['T']
        d=dict(region=c,cls=cl)
        def add(name,x,lg=False):
            d[name+'_16'],d[name+'_50'],d[name+'_84']=[q(x,w,p) for p in (.16,.5,.84)]; d[name+'_lo95'],d[name+'_hi95']=q(x,w,.025),q(x,w,.975)
        add('logNHCN_dv150',P['logNH']); add('logNHCN_percloud_dv50',P['logNH']-np.log10(3)); add('logNHCO+_dv150',P['logNO'])
        add('logn',P['logn_H']); add('T',T); add('logPk',np.log10(nH2*T)); add('logPturb23_k',np.log10(nH2*2.8*1.6726e-24*(2.3e6)**2/1.380649e-16)); add('tau32H',P['tau32H'])
        for X in (1e-8,1e-9,1e-7):
            NH2=10**P['logNH']/X; add(f'logNH2_X{X:.0e}',np.log10(NH2)); add(f'Sigma_X{X:.0e}',NH2*MSUN_PC2_PER_CM2)
            add(f"path_pc_X{X:.0e}",NH2/nH2/PC)       # single zone along the line of sight: N(H2)/n
        d['frac_Sigma_gt_5500_X1e-8']=float(w[10**P['logNH']/1e-8*MSUN_PC2_PER_CM2>5500].sum()); d['frac_Sigma_gt_5500_X1e-7']=float(w[10**P['logNH']/1e-7*MSUN_PC2_PER_CM2>5500].sum())
        d['frac_P_ge_1e7']=float(w[nH2*T>=1e7].sum()); d['frac_tau32H_gt3']=float(w[P['tau32H']>3].sum()); d['frac_tau32H_gt10']=float(w[P['tau32H']>10].sum()); d['tau32H_lo95_onesided_5pct']=q(P['tau32H'],w,.05)
        out.append(d)
R=pd.DataFrame(out); R.round(3).to_csv(O+'phys_posteriors.csv',index=False)
pd.set_option('display.width',250,'display.max_columns',60)
print(R[['region','cls','logNHCN_dv150_50','logNHCN_dv150_lo95','logNHCN_dv150_hi95','logNHCN_percloud_dv50_50','logn_50','logn_lo95','logn_hi95','T_50','T_lo95','T_hi95','logPk_50','logPk_lo95','logPk_hi95','tau32H_50','tau32H_lo95','tau32H_hi95','tau32H_lo95_onesided_5pct','frac_tau32H_gt3','frac_tau32H_gt10']].round(2).to_string())
print(R[['region','cls','Sigma_X1e-08_50','Sigma_X1e-08_lo95','Sigma_X1e-08_hi95','Sigma_X1e-07_50','path_pc_X1e-08_50','path_pc_X1e-08_lo95','path_pc_X1e-08_hi95','frac_Sigma_gt_5500_X1e-8','frac_Sigma_gt_5500_X1e-7','frac_P_ge_1e7']].round(2).to_string())
