# Q16r2 (8): excitation contrast between N, S and Bridge from the MATCHED-BEAM ladders only (HCN, HCO+ jointly; 4 ladders: 3-2/2-1 in the Xb55 beam, 4-3/3-2 in the Xb4f beam), analysed with p0012's RADEX
# machinery (q16_joint_lib: classes A, B, F2 pooled; Tkin 10-150 K and log n 3-7 flat grids, 14% per ladder, Tbg 2.73 K).  The R32/1-0 data are dropped for all three regions (no J=1-0 values exist for the Bridge): a control test.
# Variants: 'aper' (observed matched-beam ladders) and 'rec' (x1.07 on 4-3/3-2 and x1.11 on 3-2/2-1: Q12 recovery corrections, common to the regions).  Output: marginal T and n of each region, and the joint distribution
# of Tk_X/Tk_Y and n_X/n_Y (regions independent, shared systematics ignored: the recovery variants bracket them).   usage: python3 -I s83_q16_contrast.py -> out/q16b/s83.txt/json
import sys,os,json; os.environ['Q2_TMIN']='0'; sys.path.insert(0,'/workspace/scripts')
import numpy as np
import q16_joint_lib as J
S72=json.load(open('/workspace/out/q16/s72.json')); OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
def apert(c):
    k='Xb4f native 5.72x3.54 PA82 (J3,J4)|%s|%s'; k2='Xb55 native 8.64x6.00 PA-90 (J2,J3,J4)|%s|%s'; out=[]
    for mol in ('HCN','HCO+'):
        a=S72[k%(mol,c)]['tb']; b=S72[k2%(mol,c)]['tb']; out+=[b['J3']/b['J2'],a['J4']/a['J3']]
    return np.array(out)
REC=np.array([1.11,1.07,1.11,1.07]); VAR='aper'
def data(c):
    v=apert(c)*(REC if VAR=='rec' else 1.0); d=np.zeros((8,2)); d[:4,0]=v; d[:4,1]=np.hypot(0.14*v,0.03*v); d[4:,0]=1.0; d[4:,1]=1e9; return d
J.data=data
names=['A both thick, kappa=1','B both thick, same Tex (0.8-1.25), kappa free','F2 separate n, kappa=1, both thick']
res={}
for VAR in ('aper','rec'):
    post={}
    for c in ('N','S','Bridge'):
        p=J.posterior(c,'all',names,prior='env',thr=1e-5)
        n=0.5*(p['logn_H']+p['logn_O']); w=p['w']
        # 2D marginal on (T, n) grid (n rounded to the grid)
        Ts=np.unique(p['T']); H={}
        for T in Ts:
            m=p['T']==T; H[T]=w[m].sum()
        tot=sum(H.values()); pT={int(T):H[T]/tot for T in Ts}
        post[c]=dict(T=p['T'],n=n,w=w/w.sum(),nH=p['logn_H'])
        q=lambda x,qq: J.wq(x,w,qq)
        P(f'[{VAR}] {c:6s} data {np.round(data(c)[:4,0],3)}  Tk 16/50/84 = {q(p["T"],.16):.0f}/{q(p["T"],.5):.0f}/{q(p["T"],.84):.0f} K (2.5-97.5%: {q(p["T"],.025):.0f}-{q(p["T"],.975):.0f})  log n(H,O mean) 16/50/84 = {q(n,.16):.2f}/{q(n,.5):.2f}/{q(n,.84):.2f}  P(T>=80K)={sum(v for k,v in pT.items() if k>=80):.2f}  P(T<=30K)={sum(v for k,v in pT.items() if k<=30):.2f}')
        res[f'{VAR}|{c}']=dict(pT=pT,T=[q(p['T'],x) for x in (.16,.5,.84)],logn=[q(n,x) for x in (.16,.5,.84)])
    rng=np.random.default_rng(1); NS=200000
    def draw(c):
        i=rng.choice(len(post[c]['w']),NS,p=post[c]['w']); return post[c]['T'][i],post[c]['n'][i]
    D={c:draw(c) for c in post}
    for a,b in (('N','S'),('N','Bridge'),('S','Bridge')):
        rT=D[a][0]/D[b][0]; rn=10**(D[a][1]-D[b][1])
        P(f'[{VAR}] {a}/{b}: Tk ratio 16/50/84 = {np.percentile(rT,16):.2f}/{np.percentile(rT,50):.2f}/{np.percentile(rT,84):.2f}; P(Tk_{a}>Tk_{b})={np.mean(rT>1):.2f}; P(ratio>1.4)={np.mean(rT>1.4):.2f}   n ratio 16/50/84 = {np.percentile(rn,16):.2f}/{np.percentile(rn,50):.2f}/{np.percentile(rn,84):.2f}; P(n_{a}>n_{b})={np.mean(rn>1):.2f}; P(n ratio>2)={np.mean(rn>2):.2f}; P(either Tk or n higher by >1.4/2)={np.mean((rT>1.4)|(rn>2)):.2f}; P(Tk ratio>1 and n ratio<1)={np.mean((rT>1)&(rn<1)):.2f}')
        # joint 2D histogram of (log Tk ratio, log n ratio)
        Hh,xe,ye=np.histogram2d(np.log10(rT),np.log10(rn),bins=[np.linspace(-1,1,9),np.linspace(-2,2,9)],weights=None)
        res[f'{VAR}|{a}/{b}']=dict(P_T_gt=float(np.mean(rT>1)),P_n_gt=float(np.mean(rn>1)),Tratio=np.percentile(rT,[16,50,84]).tolist(),nratio=np.percentile(rn,[16,50,84]).tolist(),
                                   joint_logTratio_edges=xe.tolist(),joint_lognratio_edges=ye.tolist(),joint=(Hh/NS).tolist())
json.dump(res,open(OUT+'s83.json','w'),indent=1,default=float); open(OUT+'s83.txt','w').write('\n'.join(log))
