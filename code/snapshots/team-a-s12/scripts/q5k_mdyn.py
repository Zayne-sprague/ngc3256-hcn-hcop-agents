# Q5k: dynamical mass within R from the beam-forward-modelled rotation curve (q5i), MC over inclination (N 20-40 deg, S 70-90 deg; Sakamoto+2014: N nearly face-on, S almost edge-on)
# python3 -I scripts/q5k_mdyn.py -> out/q5h/q5k.log
import json,glob,numpy as np
G=4.30091e-3; PC=194.; rng=np.random.default_rng(5); log=open('/workspace/out/q5h/q5k.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
INC={'N':(20,40),'S':(70,90)}
for f in sorted(glob.glob('/workspace/out/q5h/q5i_*_*_*.json')):
    j=json.load(open(f)); src=j['src']; R=np.array(j['R']); sg=np.array(j['sig_int']); vl=np.array(j['vlos'])
    mode=f.split('q5i_')[1].split('_')[0]
    for Rq in (0.25,0.5,1.0):
        i=int(np.argmin(abs(R-Rq))); v=vl[1][i]; ev=0.5*(vl[2][i]-vl[0][i]); s=sg[1][i]; es=0.5*(sg[2][i]-sg[0][i])
        n=20000; inc=np.deg2rad(rng.uniform(*INC[src],n)); V=np.clip(rng.normal(v,ev,n),0,None)/np.sin(inc); S=np.clip(rng.normal(s,es,n),0,None)
        Rpc=Rq*PC; M1=V**2*Rpc/G; M2=(V**2+3*S**2)*Rpc/G   # rotation only; with pressure support (sigma isotropic)
        q=lambda m:np.percentile(m,[16,50,84])/1e8
        pr(f'{mode:5s} {src} {j["line"]:5s} R={Rq}" ({Rpc:.0f} pc): vlos={v:.0f}+-{ev:.0f} sigma_int={s:.0f}+-{es:.0f}  Vc(i-MC)={np.percentile(V,[16,50,84]).round(0)} km/s ; Mdyn(rot)={q(M1).round(2)}e8 ; Mdyn(rot+3sigma^2)={q(M2).round(2)}e8 Msun')
