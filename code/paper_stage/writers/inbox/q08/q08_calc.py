# q08: N/S contrast of H13CN/H13CO+ (1-0, matched 2.2" beam, 2015.1.00993 + X961) from out/q47/ratios10.json (scripts/q47_ratios10.py)
import json,numpy as np
from scipy.stats import chi2,norm
d=json.load(open('/workspace/out/q47/ratios10.json')); rng=np.random.default_rng(1); n=400000
def g(k,m): v=d[f'{k}|{m}']; return v['R'],v['err'],v['I12_Jybkms']
for m in 'TAP':
    s={}
    for r in 'NS':
        (a,ea,ia),(b,eb,ib)=g(f'HCN|{r}',m),g(f'HCO+|{r}',m)  # 13C/12C flux ratios; I12 integrals
        # H13CN/H13CO+ (Jy ratio) = (a*ia)/(b*ib); Tb ratio = that*(nu13CO+/nu13CN)^2 (cancels in N/S)
        s[r]=(a,ea,b,eb,ia,ib)
    out=[]
    for r in 'NS':
        a,ea,b,eb,ia,ib=s[r]; x=rng.normal(a,ea,n); y=rng.normal(b,eb,n)
        out.append((x,y,ia/ib,a*ia/(b*ib)))
    # 12C HCN/HCO+ (same beam, same estimator): ia/ib
    rN=out[0][0]/out[0][1]*out[0][2]; rS=out[1][0]/out[1][1]*out[1][2]
    ok=(out[0][1]>0)&(out[1][1]>0)&(out[1][0]>0)
    c=(rN/rS)[ok]; q=np.percentile(c,[16,50,84,2.5])
    print(m,'H13CN/H13CO+ N=%.2f S=%.2f (point, Jy ratio)'%(out[0][3],out[1][3]),'12C HCN/HCO+ (Jy-ratio) N=%.2f S=%.2f contrast %.2f'%(out[0][2],out[1][2],out[0][2]/out[1][2]),
      'thin contrast point %.2f; MC (valid draws %.2f) median %.2f [16,84]=%.2f-%.2f, 2.5pct %.2f; P(c>1)=%.2f P(c>2)=%.2f'%(out[0][3]/out[1][3],ok.mean(),q[1],q[0],q[2],q[3],(c>1).mean(),(c>2).mean()))
    # Gaussian error propagation of point value
    def rel(r): a,ea,b,eb,*_=s[r]; return np.hypot(ea/a,eb/b)
    p=out[0][3]/out[1][3]; print('   linear propagation: %.2f +- %.2f'%(p,p*np.hypot(rel('N'),rel('S'))))
for dof in (1,2,3,4):
    print('dof',dof,[ '%.2f sigma'%norm.isf(chi2.sf(x,dof)/2) for x in (5,8)])
