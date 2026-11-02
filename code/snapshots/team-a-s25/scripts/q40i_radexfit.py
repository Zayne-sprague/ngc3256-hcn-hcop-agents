# Q40i: which (T,n,N(HCN)/dv,X=N(HCN)/N(HCO+)) reproduce R10,R21,R32 (matched-uv, common-beam, T_B scale) of N,S,E1,W; can the clump deficit be abundance (X), density (n), column/opacity (N) or T?
# Hypothesis test per pair (region a, clump b): share all params (H0) vs let one parameter differ between the two (H_X, H_n, H_N, H_T): min joint chi2 over the grid (6 data).
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
import q40h_radex as R
T=json.load(open('/workspace/out/q40b/q40f_table.json'))
H=R.table('hcn',R.LNH); O=R.table('hco+',R.LNO)
LX=np.arange(-1,1.01,0.25)
nT,nn,nH=H.shape[:3]
# model ratio cube RJ[J][T,n,NH,X] and tau
RJ=np.full((4,nT,nn,nH,len(LX)),np.nan); TAU=np.full((nT,nn,nH,len(LX)),np.nan)
for ix,lx in enumerate(LX):
    for ih,lh in enumerate(R.LNH):
        lo=lh-lx; io=np.where(np.abs(R.LNO-lo)<1e-6)[0]
        if len(io)==0: continue
        io=io[0]
        for J in range(4): RJ[J,:,:,ih,ix]=H[:,:,ih,J,0]/O[:,:,io,J,0]
        TAU[:,:,ih,ix]=H[:,:,ih,2,1]   # tau(HCN 3-2)
def obs(reg):
    d={}
    for J,nm in ((0,'Michiyama'),(1,'Xb55'),(2,'Xb4f')):
        c,e,s=T[f'R|{reg}|{nm}']; d[J]=(c,np.sqrt(e*e+s*s+(0.05*c)**2))
    return d
def chi2(reg):
    o=obs(reg); c=np.zeros(RJ.shape[1:])
    for J,(v,e) in o.items(): c+=((RJ[J]-v)/e)**2
    if reg=='E1': c+=(np.maximum(RJ[3]-0.15,0)/0.05)**2   # 4-3 HCN 3-sigma upper limit R43<0.15
    return np.where(np.isfinite(c),c,1e9)
CH={r:chi2(r) for r in ('N','S','N+S','E1','W')}
log=[]
def P(s=''): print(s); log.append(s)
P('R10, R21, R32 used (value, total sigma):'); 
for r in CH: P(f'  {r}: '+'; '.join(f'{v:.3f}+-{e:.3f}' for v,e in obs(r).values()))
P('\n## per region: allowed set chi2<=chi2min+4 and chi2<=7.8 (3 data, 4 params: underdetermined, so ranges not best fit)')
LT,LNr,LNHa=R.TK,R.LN,R.LNH
for r,c in CH.items():
    m=c.min(); ok=(c<=7.8)&(c<=m+4); idx=np.argwhere(ok)
    w=np.exp(-0.5*(c-m)); w/=w.sum()
    def marg(ax,vals):
        p=w.sum(axis=tuple(i for i in range(4) if i!=ax)); cdf=np.cumsum(p); return [vals[np.searchsorted(cdf,q)] for q in (0.16,0.5,0.84)]
    P(f' {r}: chi2min {m:.2f}; allowed cells {ok.sum()}; log X (HCN/HCO+ column) [16,50,84%] {marg(3,LX)} range {LX[idx[:,3]].min():+.2f}..{LX[idx[:,3]].max():+.2f}; log n {marg(1,LNr)} range {LNr[idx[:,1]].min():.2f}..{LNr[idx[:,1]].max():.2f}; T {marg(0,LT)}; log N(HCN)(dv150) {marg(2,LNHa)} range {LNHa[idx[:,2]].min():.2f}..{LNHa[idx[:,2]].max():.2f}; tau32(HCN) median {np.median(TAU[ok]):.1f} range {np.nanmin(TAU[ok]):.1f}-{np.nanmax(TAU[ok]):.1f}')
P('\n## pair test: min joint chi2 (6 data; E1 incl. R43 limit) when sharing everything (H0) or letting one parameter differ; larger = worse')
names={0:'T',1:'n',2:'N',3:'X'}
out={}
for a,b in (('N','E1'),('S','E1'),('N+S','E1'),('N','W'),('S','W'),('N+S','W')):
    ca,cb=CH[a],CH[b]; res={'H0 all shared':(ca+cb).min()}
    for ax in range(4):
        res[f'{names[ax]} differs']=(ca.min(axis=ax)+cb.min(axis=ax)).min()
    # which X ratio is required when only X differs: dlogX at best
    ca3=ca.min(axis=3); 
    bestX=None
    # ratio of best X in a vs b at the best joint (T,n,N)
    tot=ca.min(axis=3)+cb.min(axis=3); i=np.unravel_index(tot.argmin(),tot.shape); xa=LX[ca[i].argmin()]; xb=LX[cb[i].argmin()]
    P(f' {a} vs {b}: '+'; '.join(f'{k}: {v:.1f}' for k,v in res.items())+f'; H_X best dlogX({b}-{a}) = {xb-xa:+.2f} at T={LT[i[0]]:.0f} log n={LNr[i[1]]:.2f} logN={LNHa[i[2]]:.2f}')
    out[f'{a}|{b}']={**res,'dlogX':float(xb-xa)}
# which individual quantity is needed? marginal difference of logX posterior, and an opacity statement: fraction of E1 allowed cells with tau32(HCN)>1,>3,>10
for r in ('N','S','E1','W'):
    c=CH[r]; ok=(c<=c.min()+4)&(c<=7.8)
    t=TAU[ok]; P(f' {r}: allowed cells with tau32(HCN)>1: {(t>1).mean():.2f}, >3: {(t>3).mean():.2f}, >10: {(t>10).mean():.2f}; N(HCN)/dv range log {np.log10((10**LNHa[np.argwhere(ok)[:,2]]/150).min()):.2f}..{np.log10((10**LNHa[np.argwhere(ok)[:,2]]/150).max()):.2f} (cm-2 / km/s)')
# optically thin limit: tau limited by H13CN (N: HCN/H13CN(4-3)>8.3-9.5 -> tau43(HCN)<=3.2/5.5/10 for 12C/13C 30/50/90)
P('\n## opacity limit from H13CN (nuclei only; no clump H13CN limit exists): see report entries; tau(HCN 4-3)<=3-10 (N), looser (S).')
P(f'## a thin (tau32<0.3) model: allowed? chi2min over cells with tau32(HCN)<0.3: '+'; '.join(f'{r}: {np.where(TAU<0.3,CH[r],1e9).min():.1f}' for r in CH))
P('## a thick (tau32>3) model: chi2min over cells with tau32(HCN)>3: '+'; '.join(f'{r}: {np.where(TAU>3,CH[r],1e9).min():.1f}' for r in CH))
open('/workspace/out/q40b/q40i_radexfit.log','w').write('\n'.join(log)); json.dump(out,open('/workspace/out/q40b/q40i_radexfit.json','w'),indent=1)
# required offsets under each single-parameter hypothesis: range over shared-cell family with joint chi2 <= min+4 (dlog = clump - reference)
P('\n## required offset (clump minus reference) when exactly one parameter differs; range over the allowed shared family (joint chi2<=min+4); median in brackets')
for a,b in (('N','E1'),('S','E1'),('N+S','E1'),('N','W'),('S','W'),('N+S','W')):
    ca,cb=CH[a],CH[b]; line=[]
    for ax,vals,lab in ((3,LX,'dlogX'),(1,R.LN,'dlogn'),(2,R.LNH,'dlogN')):
        tot=ca.min(axis=ax)+cb.min(axis=ax); ok=np.argwhere(tot<=tot.min()+4)
        d=[]
        for i in ok:
            ia=np.take(ca,0,axis=0)  # placeholder to keep shape logic simple
            sl=[slice(None)]*4
            idx=list(i); full=list(range(3))
            # build index tuple for the remaining 3 axes -> insert slice at ax
            tup=[None]*4; k=0
            for dim in range(4):
                if dim==ax: tup[dim]=slice(None)
                else: tup[dim]=idx[k]; k+=1
            va=ca[tuple(tup)]; vb=cb[tuple(tup)]
            d.append(vals[vb.argmin()]-vals[va.argmin()])
        d=np.array(d); line.append(f'{lab} {d.min():+.2f}..{d.max():+.2f} [{np.median(d):+.2f}]')
    P(f' {b} - {a}: '+'; '.join(line))
open('/workspace/out/q40b/q40i_radexfit.log','w').write('\n'.join(log))
