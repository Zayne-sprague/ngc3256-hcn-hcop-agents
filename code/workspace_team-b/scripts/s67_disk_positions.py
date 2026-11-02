# Q3 round 4: single-position fits at the disk positions E1, E2, W with the converged grid (0.1 dex in N, A), Q12 recovery classes, PB extra terms (s55a), band 10/10/10%; posterior vs prior.
# W: HCN 4-3 (2.7 sigma) dropped (error set to 500%). usage: python3 -I s67_disk_positions.py -> out/q3d/s67_disk_positions.txt/json
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s53_ladder_fit as L, s62_lib as S
OUT='/workspace/out/q3d/'; log=[]; res={}
def P(s=''): print(s,flush=True); log.append(s)
F=S.Grid(); T=10**F.lT; W0=S.prior_w(F,'logT300')
PJ=json.load(open('/workspace/out/q3b/s55a_pb_snr.json')); PBJ=PJ['pb']; PT=PJ['pbtest']
def pbextra(pos):
    pb=np.array(PBJ[pos]); s=0.3*(1-pb); ab=[PT[f'{b}|{pos}'] for b in ('Xb55 (B5, J=2)','Xb4f (B6, J=3)','Xb47 (B7, J=4)')]
    return np.r_[np.maximum(s,0.5*np.abs(1-np.array([a[0] for a in ab]))),np.maximum(s,0.5*np.abs(1-np.array([a[1] for a in ab])))]
BT=[10,30,60,120,250,500.1]; BN=[3,4,4.5,5,6,8.01]
for pos,cls in (('N','centre'),('E1','E1'),('E2','E2'),('W','W')):
    I,e=L.DATA[pos]; e=e.copy()
    if pos=='W': e[2]=5.0
    p=dict(y=np.log(I),rel=e,cls=cls,extra=(pbextra(pos) if pos!='N' else None))
    for tag,kw in (('band 10/10/10%',dict(sig_b=(.1,.1,.1))),('ratios only (band scales free)',dict(sig_b=(10,10,10),use_rec=False))):
        o=S.gfit(F,[p],**kw); Ps=S.Post(F,o['lnpost'],W0,o); d=Ps.summary(); a,b=Ps.marg1d(T,BT); c,e2=Ps.marg1d(F.ln,BN)
        kl=lambda q,r:float(np.sum(q[q>0]*np.log(q[q>0]/r[q>0])))
        P(f'[{pos}] {tag}: log n {d["log n"][2]:.2f} [{d["log n"][1]:.2f},{d["log n"][3]:.2f}]; Tk {d["Tk"][2]:.0f} [{d["Tk"][1]:.0f},{d["Tk"][3]:.0f}] (5%: {Ps.wq(T,(5,))[0]:.0f}); tau(HCN3-2) {10**d["tau(HCN3-2)"][2]:.1f}; P(n>1e5)={d["P"]["n>1e5"][0]:.2f} (prior {d["P"]["n>1e5"][1]:.2f}); P(Tk>50)={d["P"]["Tk>50"][0]:.2f} (prior {d["P"]["Tk>50"][1]:.2f}); P(Tk<30)={d["P"]["Tk<30"][0]:.2f} (prior {d["P"]["Tk<30"][1]:.2f}); KL(Tk)={kl(a,b):.2f} KL(n)={kl(c,e2):.2f} nats; min chi2 {o["chi"].min():.2f}')
        res[f'{pos}|{tag}']=dict(summary=d,KL=(kl(a,b),kl(c,e2)),Tk5=Ps.wq(T,(5,))[0])
json.dump(res,open(OUT+'s67_disk_positions.json','w'),indent=1,default=float); open(OUT+'s67_disk_positions.txt','w').write('\n'.join(log))
