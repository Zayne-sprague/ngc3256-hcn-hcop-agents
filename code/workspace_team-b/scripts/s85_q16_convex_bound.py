# Q16r2 (10): (a) share of the Xb4f-beam 3-2 flux at S (and N) that comes from the compact core (Xb4d template, r<rc, convolved to the Xb4f beam, >4 sigma pixels) ; (b) the lower bound on the S core T43/T32
# if all non-core gas has T43/T32 <= rho_max (N aperture value 0.450 for HCN, 0.565 HCO+, or a hotter 0.6) -> rho_core >= (rho_ap - (1-f) rho_max)/f ; (c) Bridge subtraction table (w = share of Bridge-like gas).
# usage: python3 -I s85_q16_convex_bound.py -> out/q16b/s85.txt
import sys; sys.path.insert(0,'/workspace/scripts')
from s79_q16_injection import *
OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
S72=json.load(open('/workspace/out/q16/s72.json')); k='Xb4f native 5.72x3.54 PA82 (J3,J4)|%s|%s'
for mol in ('HCN','HCO+'):
    D=setup(mol,mol,ring=None)
    ap={p:S72[k%(mol,p)]['tb']['J4']/S72[k%(mol,p)]['tb']['J3'] for p in ('N','S','Bridge')}
    fl=[]; 
    for rc in (1.0,1.5,2.0,2.5):
        e,n=D['e'],D['n']; msk=lambda p: np.hypot(e-POS[p][0],n-POS[p][1])<rc
        out={}
        for p in ('N','S'):
            core=to_coarse(conv_map(np.where(msk(p),D['X'],0),D['ps'],D['b0'],B0),'4d3',D['mol_t'],D['ec'],D['nc'])
            out[p]=float(sample_arr(core,D,p)/sample_arr(D['d3'],D,p))
        fl.append((rc,out))
        P(f'{mol:5s} compact-core share of the Xb4f-beam 3-2 flux, rc={rc}": N {out["N"]:.2f}  S {out["S"]:.2f}   (Xb4d recovers ~0.8-0.9 of the Xb4f flux at the peaks (Q12), so the share is a lower limit by ~that factor)')
    for rmax in (ap['N'],ap['N']+0.1,ap['N']+0.2):
        for rc,out in fl:
            f=out['S']/0.85   # correct for the Xb4d recovery 0.85 (Q12: 0.77-0.87 HCN/HCO+ at S) -> share of the S aperture flux in the compact core, bounded to <=1
            f=min(f,1.0); b=(ap['S']-(1-f)*rmax)/f
            P(f'   {mol:5s} S aperture {ap["S"]:.3f}; contaminants <= {rmax:.2f} ; core share f={f:.2f} (rc {rc}") -> S core ratio >= {b:.3f}')
    for w in (0.0,0.1,0.2,0.3):
        P(f'   {mol:5s} Bridge subtraction w={w:.1f} (Bridge ladder {ap["Bridge"]:.3f}): S core = {(ap["S"]-w*ap["Bridge"])/(1-w):.3f};  with N-like contaminant ({ap["N"]:.3f}): {(ap["S"]-w*ap["N"])/(1-w):.3f}')
    P(f'   {mol:5s} contaminant share needed to bring S from {ap["S"]:.3f} to team-a core: w=(ap-core)/(rho_c-core): '+', '.join(f'rho_c={r:.2f}: w={(ap["S"]-c)/(r-c):.2f}' for c in ((0.268 if mol=="HCN" else 0.442),) for r in (ap['N'],ap['N']+0.1,ap['N']+0.2)))
# (d) consistency of p0012's own component ladders with the observed S aperture: predicted rho_ap(S) = f rho_core(S) + (1-f) rho_rest, rho_rest in [Outer, Bridge] of p0012
P('\n# (d) p0012 component ladders predict the S aperture ratio (f = compact-core share /0.85 Xb4d recovery; rest = p0012 Outer .. Bridge ladder)')
TA={'HCN':dict(S=0.268,Bridge=0.363,Outer=0.280),'HCO+':dict(S=0.442,Bridge=0.464,Outer=0.299)}
for mol in ('HCN','HCO+'):
    ap=S72[k%(mol,'S')]['tb']['J4']/S72[k%(mol,'S')]['tb']['J3']; err=0.011 if mol=='HCN' else 0.005
    for f in (0.38,0.52,0.67,0.78):
        lo=f*TA[mol]['S']+(1-f)*min(TA[mol]['Bridge'],TA[mol]['Outer']); hi=f*TA[mol]['S']+(1-f)*max(TA[mol]['Bridge'],TA[mol]['Outer'])
        P(f'  {mol:5s} f={f:.2f}: predicted S aperture T43/T32 {lo:.3f}-{hi:.3f} vs observed {ap:.3f}+-{err:.3f} (stat): shortfall {ap-hi:.3f} to {ap-lo:.3f}')
open(OUT+'s85.txt','w').write('\n'.join(log))
