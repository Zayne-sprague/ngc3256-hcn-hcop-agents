# Q16 (3): re-run team-a's p0012 joint non-LTE fit (their cube code, copied to q16_joint_lib.py with paths -> out/q16/in; grids Tbg 2.73/5/10 K) for the S and N nuclei with PERTURBED data vectors,
# and report the profile chi2(T) (min over n_H,n_O,N_H,N_O,kappa,g within classes A, B, F2 = p0012's pooled classes; prior T>=10 K grid edge, no floor).
# Data vector d0 HCN 3-2/2-1, d1 HCN 4-3/3-2, d2 HCO+ 3-2/2-1, d3 HCO+ 4-3/3-2 (deblended Tb ladders), d4 R32, d5 R10(Harada), d6 HCN 3-2/1-0, d7 HCO+ 3-2/1-0.
# Variants: base (reproduces p0012), leave-one-out, ladders replaced by our matched-beam aperture values (s72; Xb4f beam for 4-3/3-2, Xb55 beam for 3-2/2-1), mid-point between team-a cores and our aperture values,
# error inflation of the ladders (14% -> 25%), a common band-scale error, Tbg 5 and 10 K.   usage: python3 -I s73_q16_Sfit_variants.py <variant-set> -> out/q16/s73_<set>.csv
import sys,os,json; os.environ['Q2_TMIN']='0'
sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
which=sys.argv[1]
suf={'base':'','tbg5':'_tbg5','tbg10':'_tbg10'}.get(which.split('+')[0],'')
os.environ['Q2_GRIDSUF']=suf
import q16_joint_lib as J
orig=J.data
S72=json.load(open('/workspace/out/q16/s72.json'))
def apert(c):
    k='Xb4f native 5.72x3.54 PA82 (J3,J4)|%s|%s'; k2='Xb55 native 8.64x6.00 PA-90 (J2,J3,J4)|%s|%s'
    out=[]
    for mol in ('HCN','HCO+'):
        a=S72[k%(mol,c)]['tb']; b=S72[k2%(mol,c)]['tb']
        out+= [b['J3']/b['J2'], a['J4']/a['J3']]
    return out   # d0,d1,d2,d3 order: HCN 3-2/2-1, HCN 4-3/3-2, HCO+ 3-2/2-1, HCO+ 4-3/3-2
def make(c,var):
    d=orig(c).copy()
    if var=='base': return d
    if var.startswith('loo'):
        i=int(var[3:]); d[i,1]=1e6; return d       # effectively drops datum i (huge error)
    if var in ('aper','mid'):
        ap=apert(c)
        for i in range(4):
            tgt=ap[i] if var=='aper' else 0.5*(ap[i]+d[i,0]); d[i,0]=tgt; d[i,1]=np.hypot(0.14*tgt,abs(ap[i]-orig(c)[i,0])*(0.5 if var=='mid' else 0.0)+0.03*tgt)
        return d
    if var in ('aper43','aper32'):
        ap=apert(c)
        for i in ((1,3) if var=='aper43' else (0,2)): d[i,0]=ap[i]; d[i,1]=np.hypot(0.14*ap[i],0.03*ap[i])
        return d
    if var.startswith('infl'):
        f=float(var[4:]); d[:4,1]=np.hypot(d[:4,1],f*d[:4,0]); return d
    if var=='ladders_off': d[:4,1]=1e6; return d
    raise SystemExit(var)
rows=[]
for c in ('S','N'):
    for var in sys.argv[2:]:
        J.data=lambda cc,c=c,var=var: make(cc,var) if cc==c else orig(cc)
        d,use,chi,cls,aux=J.cube(c,'all','env'); sh=aux['sh']
        for nm in ('A','B','F2'):
            full={'A':'A both thick, kappa=1','B':'B both thick, same Tex (0.8-1.25), kappa free','F2':'F2 separate n, kappa=1, both thick'}[nm]
            mask,iks=cls[full]
            for iT,T in enumerate(J.TK):
                v=min(float(np.where(mask,chi[ik],1e9)[iT].min()) for ik in iks)
                rows.append(dict(region=c,variant=var,grid=which,cls=nm,T=T,chi2=v))
        dd=pd.DataFrame(rows); q=dd[(dd.region==c)&(dd.variant==var)].groupby('T').chi2.min()
        print(which,c,var,'data',np.round(d[:4,0],3),'| profile chi2(T):',' '.join(f'{int(T)}:{x:.1f}' for T,x in q.items()),flush=True)
pd.DataFrame(rows).to_csv(f'/workspace/out/q16/s73_{which}.csv',index=False)
