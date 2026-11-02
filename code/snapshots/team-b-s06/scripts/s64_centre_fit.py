# Q3 round 4, fix 2: centre ladder fit with SHARED calibration, one-sided limits, prior sensitivity and posterior/prior ratios.
# Nuisances shared by all positions in a fit and marginalised ANALYTICALLY (Gaussian; s62_lib.GEngine; the Monte-Carlo version with 256 Sobol draws agrees): band scales (Xb55, Xb4f, Xb47) ~ N(0,sig_b) in ln; Q12 recovery factors (centre: f2=true/obs(2-1) U[.85,.95], f4 U[1.03,1.12]) common to N,S,Bridge.
# Per-position nuisance: ln eta (flat, <=0, exact truncated-normal marginal). Per-line independent term s_ind=7% (independent per position) or, in 'all shared' mode, shared by the positions too (aperture overlap makes the positions near copies).
# Priors: flat log Tk 10-300 (baseline, as rounds 1-3), flat log Tk 10-500, flat Tk 10-500, flat n (1e3-1e8, 1e3-1e7, 1e3-1e6; Tk<=300), always flat in log N(HCO+) 11.5-19 and log A -1..1.4.
# usage: python3 -I s64_centre_fit.py -> out/q3d/s64_centre_fit.txt/json, s64_post.npz
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L, s62_lib as S
OUT='/workspace/out/q3d/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid()
def pos(k,cls='centre'):
    I,e=L.DATA[k]; return dict(y=np.log(I),rel=e,cls=cls)
PR=['logT300','logT500','linT500','linn6','linn7','linn']
T=10**F.lT
def brief(d):
    return (f"log n {d['log n'][2]:.2f} [{d['log n'][1]:.2f},{d['log n'][3]:.2f}]; Tk {d['Tk'][2]:.0f} [{d['Tk'][1]:.0f},{d['Tk'][3]:.0f}] (5%: {d['Tk5']:.0f}); log A {d['log A'][2]:.1f} [{d['log A'][1]:.1f},{d['log A'][3]:.1f}]; log N(HCO+)/dv {d['log N(HCO+)'][2]-np.log10(150):.2f} [{d['log N(HCO+)'][1]-np.log10(150):.2f},{d['log N(HCO+)'][3]-np.log10(150):.2f}]; "
            f"tau(HCN3-2) {10**d['tau(HCN3-2)'][2]:.1f} [{10**d['tau(HCN3-2)'][1]:.1f},{10**d['tau(HCN3-2)'][3]:.1f}]; log P/k {d['log P/k'][2]:.2f} [{d['log P/k'][1]:.2f},{d['log P/k'][3]:.2f}]")
def pbrief(d):
    p=d['P']; return '  '.join(f"P({k})={p[k][0]:.2f} (prior {p[k][1]:.2f})" for k in ('n>1e5','n<1e4','Tk>30','Tk>50','Tk>100','A>1','tauHCN>1','tauHCN>5'))
BT=[10,30,60,120,250,500.1]; BN=[3,4,4.5,5,6,8.01]
def ratios_tab(Ps):
    a,b=Ps.marg1d(T,BT); c,e=Ps.marg1d(F.ln,BN)
    s1=' '.join(f"[{BT[i]:.0f}-{min(BT[i+1],500):.0f} K: post {a[i]:.2f} / prior {b[i]:.2f} = {a[i]/max(b[i],1e-9):.2f}]" for i in range(len(a)))
    s2=' '.join(f"[logn {BN[i]:g}-{min(BN[i+1],8):g}: post {c[i]:.2f} / prior {e[i]:.2f} = {c[i]/max(e[i],1e-9):.2f}]" for i in range(len(c)))
    kl=lambda q,r:float(np.sum(q[q>0]*np.log(q[q>0]/r[q>0])))
    return s1,s2,kl(a,b),kl(c,e)
def run(name,positions,indep=False,**kw):
    if indep:   # old round-3 treatment: every position has its own band scales and recovery
        outs=[S.gfit(F,[p],keep_samples=False,**kw) for p in positions]; lnp=sum(o['lnpost'] for o in outs); o={'lnpost':lnp}
        o2=S.gfit(F,[positions[0]],**kw); o.update(sidx=o2['sidx'],seta=o2['seta'],sw=o2['sw'])   # eta of the first position only (indicative)
    else: o=S.gfit(F,positions,**kw)
    out={}; P(f'\n## {name}   [min chi2 of the best model (profile over eta, 1 dof) = {o.get("chi",np.array([np.nan])).min():.2f}]'); posts[name]=o['lnpost']
    for pr in PR:
        w=S.prior_w(F,pr); Ps=S.Post(F,o['lnpost'],w,o); d=Ps.summary(); d['Tk5']=Ps.wq(T,(5,))[0]; d['lnZ']=Ps.lnZ; d['edgeThi']=float(Ps.p[T>=280].sum())
        P(f' [{pr}] '+brief(d)); P('      '+pbrief(d)); P(f"      Tk>=280 K mass {d['edgeThi']:.2f} (prior {float(Ps.pr[T>=280].sum()):.2f}); log10 eta "+(f"{d['log10 eta[0]'][2]:.2f} [{d['log10 eta[0]'][1]:.2f},{d['log10 eta[0]'][3]:.2f}]" if 'log10 eta[0]' in d else ''))
        if pr in ('logT300','linn7'):
            s1,s2,k1,k2=ratios_tab(Ps); P('      posterior/prior (Tk): '+s1+f'  KL={k1:.2f} nats'); P('      posterior/prior (n):  '+s2+f'  KL={k2:.2f} nats'); d['KL']=(k1,k2)
        out[pr]=d
    return out
res={}; posts={}
for tag,kw,pl,ind in (
  ('N alone (representative), band 10/10/10%',dict(sig_b=(.1,.1,.1)),['N'],False),
  ('S alone',dict(sig_b=(.1,.1,.1)),['S'],False),
  ('Bridge alone',dict(sig_b=(.1,.1,.1)),['Bridge'],False),
  ('N alone, band 15/10/15%',dict(sig_b=(.15,.1,.15)),['N'],False),
  ('N alone, no recovery term',dict(sig_b=(.1,.1,.1),use_rec=False),['N'],False),
  ('joint N+S+Bridge, band scales and recovery SHARED (per-position 7% line term independent)',dict(sig_b=(.1,.1,.1)),['N','S','Bridge'],False),
  ('joint N+S+Bridge, half of the 7% line term (5%) shared across positions, 5% independent (apertures overlap)',dict(sig_b=(.1,.1,.1),s_ind=0.05,s_line=0.05),['N','S','Bridge'],False),
  ('joint N+S+Bridge, calibration INDEPENDENT per position (round-3 treatment, for comparison)',dict(sig_b=(.1,.1,.1)),['N','S','Bridge'],True)):
    res[tag]=run(tag,[pos(k) for k in pl],indep=ind,**kw)
json.dump(res,open(OUT+'s64_centre_fit.json','w'),indent=1,default=float); np.savez(OUT+'s64_post.npz',keys=np.array(list(posts)),**{f'k{i}':v for i,v in enumerate(posts.values())})
open(OUT+'s64_centre_fit.txt','w').write('\n'.join(log))
