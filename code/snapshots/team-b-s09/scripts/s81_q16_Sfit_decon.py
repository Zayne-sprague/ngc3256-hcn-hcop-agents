# Q16r2 (6): RADEX profile chi2(T) for S (and N) with decontaminated / recovery-corrected ladder vectors, using p0012's joint fit (q16_joint_lib = team-a q4c_joint.py), grids Tbg 2.73 / 5 K.
# Ladders at the matched native beams (s72): HCN/HCO+ 3-2/2-1 (Xb55 beam), 4-3/3-2 (Xb4f beam). Core-ratio estimate for contaminant share w (fraction of the LOWER-J flux in the aperture that belongs
# to gas with ladder rho_c):  rho_core = (rho_ap - w rho_c)/(1-w).
# Variants (all four ladders replaced; error = hypot(14%, 3%) unless stated):
#   base (p0012 core ladders) | aper (aperture, w=0) | decB10/20/30: contaminant = Bridge aperture ladder | decN10/20/30: contaminant = N aperture ladder (hottest gas seen; unphysical share at S, N response 0.003/0.13: a bound)
#   rec: aper with Q12 recovery corrections (T43/T32 true/obs 1.07 [1.03-1.12]; T32/T21 true/obs 1/0.90 = 1.11) | decN30rec: both | brk: centre = mean(aper, decN30rec), error includes half the range
#   plus the ladder-pair table: base; aper32 (only the 3-2/2-1 pair swapped); aper43 (only 4-3/3-2 swapped); aper (both)
# usage: python3 -I s81_q16_Sfit_decon.py <grid:base|tbg5> <variants...>  -> out/q16b/s81_<grid>.csv
import sys,os,json; os.environ['Q2_TMIN']='0'
sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
which=sys.argv[1]; suf={'base':'','tbg5':'_tbg5','tbg10':'_tbg10'}[which]
os.environ['Q2_GRIDSUF']=suf
import q16_joint_lib as J
orig=J.data
S72=json.load(open('/workspace/out/q16/s72.json'))
def apert(c):
    k='Xb4f native 5.72x3.54 PA82 (J3,J4)|%s|%s'; k2='Xb55 native 8.64x6.00 PA-90 (J2,J3,J4)|%s|%s'; out=[]
    for mol in ('HCN','HCO+'):
        a=S72[k%(mol,c)]['tb']; b=S72[k2%(mol,c)]['tb']; out+=[b['J3']/b['J2'],a['J4']/a['J3']]
    return np.array(out)
REC=np.array([1.11,1.07,1.11,1.07])
def make(c,var):
    d=orig(c).copy(); ap=apert(c); err=lambda v: np.hypot(0.14*v,0.03*v)
    def setv(v,e=None):
        d[:4,0]=v; d[:4,1]=err(v) if e is None else e; return d
    if var=='base': return d
    if var=='aper': return setv(ap)
    if var=='aper32': d[[0,2],0]=ap[[0,2]]; d[[0,2],1]=err(ap[[0,2]]); return d
    if var=='aper43': d[[1,3],0]=ap[[1,3]]; d[[1,3],1]=err(ap[[1,3]]); return d
    if var=='rec': return setv(ap*REC)
    if var.startswith('dec'):
        w=int(var[4:6])/100; cont=apert('Bridge' if var[3]=='B' else 'N'); v=(ap-w*cont)/(1-w)
        if var.endswith('rec'): v=v*REC
        return setv(v)
    if var=='brk':
        lo=(ap-0.3*apert('N'))/0.7*REC; hi=ap*REC; v=0.5*(lo+hi); return setv(v,np.hypot(err(v),0.5*(hi-lo)))
    raise SystemExit(var)
rows=[]
for c in ('S','N'):
    for var in sys.argv[2:]:
        if c=='N' and var.startswith('dec'): continue
        J.data=lambda cc,c=c,var=var: make(cc,var) if cc==c else orig(cc)
        d,use,chi,cls,aux=J.cube(c,'all','env')
        for nm in ('A','B','F2'):
            full={'A':'A both thick, kappa=1','B':'B both thick, same Tex (0.8-1.25), kappa free','F2':'F2 separate n, kappa=1, both thick'}[nm]
            mask,iks=cls[full]
            for iT,T in enumerate(J.TK):
                v=min(float(np.where(mask,chi[ik],1e9)[iT].min()) for ik in iks)
                rows.append(dict(region=c,variant=var,grid=which,cls=nm,T=T,chi2=v))
        dd=pd.DataFrame(rows); q=dd[(dd.region==c)&(dd.variant==var)].groupby('T').chi2.min()
        print(which,c,var,'data',np.round(d[:4,0],3),'| profile:',' '.join(f'{int(T)}:{x:.1f}' for T,x in q.items()),flush=True)
pd.DataFrame(rows).to_csv(f'/workspace/out/q16b/s81_{which}.csv',index=False)
