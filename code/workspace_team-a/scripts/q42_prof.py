# Q42 sampler-free check: profile chi2_min(Tc_N) over the other 16 parameters (differential evolution, 3 restarts) on a grid of Tc_N, joint N+S shared envelope.
# usage: q42_prof.py <mode ''|drop1|inf30> <scen>  -> out/q42/prof_<mode>_<scen>.csv
import sys,os,json
mode,scen=sys.argv[1],sys.argv[2]; os.environ['Q39_MODE']=mode; src=open('/workspace/scripts/q39_q22.py').read().split('def smc(')[0]
sys.argv=['q39_q22.py','sh','q36'+scen,'1']; exec(compile(src,'q39_q22.py','exec'))
from scipy.optimize import differential_evolution
iN=ix['lTc_N']; rows=[]
for T in [10,11,12,14,16,20,25,30,40,50,70,100,150]:
    lT=np.log10(T); idx_=[i for i in range(ND) if i!=iN]
    def f(z):
        x=np.empty(ND); x[iN]=lT; x[idx_]=z; c=chi2(x[None,:])[0]; return c
    best=1e9
    for s in range(3):
        r=differential_evolution(f,[(LO[i],HI[i]) for i in idx_],popsize=20,maxiter=400,tol=1e-8,seed=s,polish=True,init='sobol'); best=min(best,r.fun)
    rows.append((T,best)); print(mode,scen,T,round(best,2),flush=True)
pd.DataFrame(rows,columns=['TcN','chi2prof']).to_csv(f'/workspace/out/q42/prof_{mode or "base"}_{scen}.csv',index=False)
