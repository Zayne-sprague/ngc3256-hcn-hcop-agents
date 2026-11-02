# Q10 v2: 3-component (free-free + synchrotron + dust) decomposition of the 98/258 GHz continuum per annulus (1.05"), MC error propagation, R vs dust SB.
# python3 -I q10_decomp.py -> out/q10v2/decomp.log, decomp.json (needs out/q10/q10_alpha.json, out/q10v2/h40a_v2.json)
import json,numpy as np
rng=np.random.default_rng(7)
A=json.load(open('/workspace/out/q10/q10_alpha.json')); Hh=json.load(open('/workspace/out/q10v2/h40a_v2.json'))['RES']
Rr={'N':[0.805,0.726,0.599,0.474],'S':[0.390,0.364,0.358,0.390]}
ann=['0-0.5','0.5-1','1-1.5','1.5-2']
nu1,nu3=257.8,98.05; q=nu1/nu3
LCff=lambda Te,y=0.08: 6.985e3*98.11**1.1*Te**-1.15/(1+y)*1.0645
log=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s); log.append(s)
N=20000; out={}
def ffdraw(reg,scen):
    if scen=='ours':   # H40a L/C (1.05", r<1" aperture, matched filter) with Te ~ U(5000,12000)
        r=Hh[f'1.05|{reg}|0-1.0']; Te=rng.uniform(5000,12000,N); f=(r['LC']+r['eLC']*rng.standard_normal(N))/LCff(Te); return np.clip(f,0,1),None
    if scen=='mich':   # Michiyama+2020 Sect III.1: ~76% (N), ~90% (S) of 98 GHz continuum is free-free; +-0.12 adopted systematic
        mu={'N':0.76,'S':0.90}[reg]; return np.clip(mu+0.12*rng.standard_normal(N),0,1),None
    if scen=='two':    # previous 2-component: no synchrotron
        return None,None
res={}
P('3-component decomposition: S98 = ff + sync + dust fractions; ratio S258/S98 = f_ff q^-0.1 + f_s q^a_s + f_d q^a_d, q=%.3f; a_d~U(3.5,4.0) (beta 1.5-2), a_s~U(-0.9,-0.5); f_ff98 constant over annuli (assumption); flux ratio error: 5%% scale each band + stat'%q)
for scen,lab in (('ours','H40a(ours, Te 5000-12000 K, r<1")'),('mich','Michiyama+2020 f_ff (N 0.76, S 0.90, +-0.12)')):
    P(f'\n--- scenario f_ff98 from {lab}')
    for reg in 'NS':
        f_ff,_=ffdraw(reg,scen); P(f' {reg}: f_ff98 prior median {np.median(f_ff):.2f} [{np.percentile(f_ff,16):.2f},{np.percentile(f_ff,84):.2f}]')
        for i,a in enumerate(ann):
            d=A[f'{reg}_{a}']; ratio=d['s1']/d['s3']
            er=ratio*np.hypot(0.07,np.hypot(0.119/d['s1'],0.063/d['s3'])*1.0)  # 5%+5% scale, noise
            rd=ratio+er*rng.standard_normal(N)
            ad=rng.uniform(3.5,4.0,N); as_=rng.uniform(-0.9,-0.5,N)
            # ratio = f q^-0.1 + fs q^as + (1-f-fs) q^ad -> fs=(f q^-0.1+(1-f) q^ad - rd)/(q^ad - q^as)
            fs=(f_ff*q**-0.1+(1-f_ff)*q**ad-rd)/(q**ad-q**as_); fd=1-f_ff-fs
            ok=(fs>=0)&(fd>=0); rej=1-ok.mean()
            fd258=np.where(ok,fd*q**ad/rd,np.nan); fff258=np.where(ok,f_ff*q**-0.1/rd,np.nan); fs258=np.where(ok,fs*q**as_/rd,np.nan)
            dust=d['s1']*fd258
            pc=lambda x:(np.nanmedian(x),np.nanpercentile(x,16),np.nanpercentile(x,84))
            res[f'{scen}|{reg}|{a}']=dict(fd258=pc(fd258),fs98=pc(np.where(ok,fs,np.nan)),dust=pc(dust),rej=rej,fff258=pc(fff258),fs258=pc(fs258))
            P(f'  {reg} {a}": S98={d["s3"]:.2f} S258={d["s1"]:.2f} mJy/b  f_s98={pc(np.where(ok,fs,np.nan))[0]:.2f}  at 258: f_dust={pc(fd258)[0]:.2f} [{pc(fd258)[1]:.2f},{pc(fd258)[2]:.2f}] f_ff={pc(fff258)[0]:.2f} f_sync={pc(fs258)[0]:.2f}  dust SB(1.1mm)={pc(dust)[0]:.2f} [{pc(dust)[1]:.2f},{pc(dust)[2]:.2f}] mJy/b  rejected {rej:.2f}')
# two-component reference
P('\n--- R vs 1.1 mm continuum SB and vs dust SB (annulus level; the 4 annuli are not independent: descriptive)')
from scipy.stats import spearmanr
for scen in ('ours','mich'):
    for reg in 'NS':
        tot=[A[f'{reg}_{a}']['s1'] for a in ann]; du=[res[f'{scen}|{reg}|{a}']['dust'][0] for a in ann]
        # linear slope dR/d(SB) per mJy/b
        st=np.polyfit(tot,Rr[reg],1)[0]; sd=np.polyfit(du,Rr[reg],1)[0]
        P(f' {scen} {reg}: R={Rr[reg]}  total 1.1mm SB={np.round(tot,2).tolist()}  dust SB={np.round(du,2).tolist()}  dR/dSB_tot={st:+.3f}, dR/dSB_dust={sd:+.3f} per mJy/b')
    # cross-nucleus: R at equal dust SB
    dN=[res[f'{scen}|N|{a}']['dust'][0] for a in ann]; dS=[res[f'{scen}|S|{a}']['dust'][0] for a in ann]
    P(f' {scen}: dust SB overlap N {min(dN):.2f}-{max(dN):.2f} vs S {min(dS):.2f}-{max(dS):.2f} mJy/b; R(N) interpolated at S core dust SB {dS[0]:.2f}: {np.interp(dS[0],sorted(dN),[x for _,x in sorted(zip(dN,Rr["N"]))]):.2f} vs R(S core) {Rr["S"][0]:.2f}')
json.dump(res,open('/workspace/out/q10v2/decomp.json','w'),indent=1); open('/workspace/out/q10v2/decomp.log','w').write('\n'.join(log))
