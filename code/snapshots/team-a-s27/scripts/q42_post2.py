# Q42 post-processing of out/q42/<dir>/*.npz -> csv + modes txt. Mode weights: 2-component Gaussian mixture in log10 T on weight-resampled samples
# (accepted if BIC(2)<BIC(1)-10 and the components are separated by >2 sigma), cross-checked against a KDE valley split. Validated on synthetic mixtures (--validate).
import glob,json,sys,numpy as np,pandas as pd
D=sys.argv[1] if len(sys.argv)>1 and sys.argv[1]!='--validate' else 'runs2'
rs=np.random.default_rng(0)
def wq(x,w,q): o=np.argsort(x); c=np.cumsum(w[o]); c/=c[-1]; return np.interp(q,c,x[o])
def em(x,k,n_init=5,it=300):
    """1D Gaussian mixture by EM; returns (weights,means,sds,loglik)"""
    best=None
    for _ in range(n_init):
        m=np.sort(rs.choice(x,k)) if k>1 else np.array([x.mean()]); s=np.full(k,x.std()/k); p=np.full(k,1/k)
        for _ in range(it):
            L=p*np.exp(-0.5*((x[:,None]-m)/s)**2)/(s*np.sqrt(2*np.pi))+1e-300; R=L/L.sum(1,keepdims=True)
            nk=R.sum(0); p=nk/len(x); m=(R*x[:,None]).sum(0)/nk; s=np.sqrt((R*(x[:,None]-m)**2).sum(0)/nk)+1e-3
        ll=np.log(L.sum(1)).sum()
        if best is None or ll>best[3]: best=(p,m,s,ll)
    return best
def gmm_modes(x,w,n=6000):
    xs=x[rs.choice(len(x),n,p=w/w.sum())]
    p1,m1,s1,l1=em(xs,1); p2,m2,s2,l2=em(xs,2); dB=-2*l1+3*np.log(n)-(-2*l2+5*np.log(n)); sep=abs(m2[0]-m2[1])/(0.5*s2.sum())
    if dB>10 and sep>2:
        o=np.argsort(m2); return [(10**m2[i],float(p2[i])) for i in o],dB
    return [(10**m1[0],1.0)],dB
if '--validate' in sys.argv:
    for w1,m1,s1,m2,s2 in [(0.3,1.15,.08,1.6,.1),(0.5,1.2,.1,1.6,.1),(0.1,1.15,.08,1.6,.1),(1,1.4,.2,0,1),(0.3,1.3,.1,1.45,.1)]:
        k=rs.random(20000)<w1; x=np.where(k,rs.normal(m1,s1,20000),rs.normal(m2,s2,20000)) if w1<1 else rs.normal(m1,s1,20000)
        print('truth w=%.2f at %.0f K / %.0f K ->'%(w1,10**m1,10**m2),gmm_modes(x,np.ones(len(x))))
    sys.exit()
rows=[];txt=[]
for f in sorted(glob.glob(f'/workspace/out/q42/{D}/*.npz')):
    d=np.load(f,allow_pickle=True); S,w,nm=d['S'],d['w'],list(d['names']); J=json.load(open(f.replace('.npz','.json')))
    P=lambda n:S[:,nm.index(n)]; TN,TS=10**P('lTc_N'),10**P('lTc_S'); r=dict(J)
    for n,T in (('TcN',TN),('TcS',TS)): r[n+'_16'],r[n+'_50'],r[n+'_84']=wq(T,w,[.16,.5,.84])
    if 'lTe' in nm: Te=10**P('lTe'); r['Te_16'],r['Te_50'],r['Te_84']=wq(Te,w,[.16,.5,.84]); r['P50']=float(w[Te>=50].sum()); r['P_TeGtBoth']=float(w[(Te>TN)&(Te>TS)].sum())
    else:
        for z in 'NS': r['P50_'+z]=float(w[10**P('lTe_'+z)>=50].sum())
    r['P_TcN_gt_TcS']=float(w[TN>TS].sum()); dl=P('lTc_N')-P('lTc_S'); r['dlogT_16'],r['dlogT_50'],r['dlogT_84']=wq(dl,w,[.16,.5,.84])
    r['P_TcN_le30']=float(w[TN<=30].sum()); r['P_TcN_lt20']=float(w[TN<20].sum())
    key=f.split('/')[-1][:-4]; txt.append(key)
    for n in ('lTc_N','lTc_S')+(('lTe',) if 'lTe' in nm else ()):
        g,dB=gmm_modes(P(n),w); txt.append(f'  {n}: '+'; '.join(f'T~{t:.0f} K weight {wt:.2f}' for t,wt in g)+f' (dBIC {dB:.0f})')
        if n=='lTc_N': r['TcN_modes']=' | '.join(f'{t:.0f}K:{wt:.2f}' for t,wt in g)
    rows.append(r)
df=pd.DataFrame(rows); df.to_csv(f'/workspace/out/q42/q42_summary_{D}.csv',index=False); open(f'/workspace/out/q42/q42_modes_{D}.txt','w').write('\n'.join(txt))
pd.set_option('display.width',300); pd.set_option('display.max_columns',50)
print(df.drop(columns=['sec','ncall']).round(2).to_string()); print('\n'.join(txt))
