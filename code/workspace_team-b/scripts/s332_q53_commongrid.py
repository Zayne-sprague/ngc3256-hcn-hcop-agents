"""Q53 round 2: velocity-resolved X/HCO+ N/S on a COMMON velocity grid (cumulative-flux interpolation, fractional channels),
profiles with >5 sigma masking, sensitivity (v_sys, edges, baseline), leave-one-bin null test, joint narrow+broad Gaussian decomposition,
3-2/1-0 and CN hyperfine opacity checks.  usage: python3 -I s332_q53_commongrid.py -> out/q53/s332.txt/json"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
import s320_q52_chem_ratios as m
from s180_q29_load import to_tb_conv
from scipy.optimize import least_squares
OTHER={'10':{'HCN':['HCO+'],'HCO+':['HCN'],'CN_G2':['CN_G1'],'CN_G1':['CN_G2'],'HNC':[]},'32':{'HCN':[],'HCO+':[],'HNC':[]}}
TB={}
def cubeTb(J,ln,th):
    k=(J,ln,th)
    if k not in TB:
        path,f0=m.LINES[J][ln]; cut,nu,pix,bm,pxy=m.getcube(path,th)
        Tb=to_tb_conv(cut,nu,pix,bm,th,float(np.mean(nu))); TB[k]=(Tb,nu,pxy)
    return TB[k]
def raw(J,ln,th,rn):
    """pixel spectrum on the native velocity axis (ascending), v relative to nominal v_sys"""
    Tb,nu,pxy=cubeTb(J,ln,th); path,f0=m.LINES[J][ln]; pos,v0=m.POS[rn]
    x,y=pxy(pos); sp=Tb[:,int(round(y)),int(round(x))].astype(float); v=m.vel(nu,f0,v0)
    o=np.argsort(v); return v[o],sp[o]
def prof(J,ln,th,rn,order=1,rng=(600,1200),skip=None):
    v,sp=raw(J,ln,th,rn); path,f0=m.LINES[J][ln]; v0=m.POS[rn][1]
    free=(abs(v)>rng[0])&(abs(v)<rng[1])
    for o in OTHER[J][ln]:
        g0=m.LINES[J][o][1]; vg=m.vel(m.C*0+ (m.C*(1-v/m.C*0)) ,1,0) if False else None
        vo=v+ (f0-g0)/f0*m.C*0   # placeholder
    # velocity of channel frequencies relative to the other line
    nu=f0*(1-(v+v0)/m.C)
    for o in OTHER[J][ln]:
        vg=m.vel(nu,m.LINES[J][o][1],v0); free&=~((vg>-650)&(vg<650))
    if skip is not None: free&=~((v>skip[0])&(v<skip[1]))
    p=np.polyfit(v[free],sp[free],order); s=sp-np.polyval(p,v); sig=np.std(s[free])
    dv=np.gradient(v); return v,s,sig,dv,free
def integ(v,s,dv,a,b):
    """integral of channel-boxcar spectrum over [a,b] km/s, fractional channels"""
    e=np.concatenate([[v[0]-dv[0]/2],(v[1:]+v[:-1])/2,[v[-1]+dv[-1]/2]]); F=np.concatenate([[0],np.cumsum(s*dv)])
    return np.interp(b,e,F)-np.interp(a,e,F)
def binI(J,ln,th,rn,a,b,vc=0,**kw):
    v,s,sig,dv,_=prof(J,ln,th,rn,**kw); I=integ(v,s,dv,a+vc,b+vc); return I,sig*np.sqrt(abs(b-a)*np.mean(dv[abs(v)<300]))
out=[]; res={}
def P(*a):
    t=' '.join(str(x) for x in a); out.append(t); print(t,flush=True)
def ratio(J,th,rn,X,a,b,vc=0,**kw):
    x,ex=binI(J,X,th,rn,a,b,vc,**kw); h,eh=binI(J,'HCO+',th,rn,a,b,vc,**kw)
    q=x/h; return q,abs(q)*np.hypot(ex/x,eh/h),h,eh,x,ex
# A. profiles
EDG=[-350,-200,-100,-50,0,50,100,200,350]
for J,th in (('10',2.5),('32',2.5),('32',1.4)):
    for X in ['HCN','HNC']+(['CN_G2'] if J=='10' else []):
        P(f'\n== A. profile J={J} {th}" {X}/HCO+ per bin (|bins| with both lines >5 sigma shown, else limit on the ratio from 3 sigma of the weaker)')
        for rn in 'NS':
            for a,b in zip(EDG[:-1],EDG[1:]):
                q,e,h,eh,x,ex=ratio(J,th,rn,X,a,b)
                ok=(h>5*eh) and (x>5*ex)
                s=f'{q:.2f}+-{e:.2f}' if ok else (f'<{(max(x,0)+3*ex)/h:.2f} (lim)' if h>5*eh else 'HCO+<5sig')
                P(f'  {rn} {a:5d}..{b:4d} HCO+={h:7.1f}+-{eh:.1f} {X}={x:7.1f}+-{ex:.1f} ratio={s}')
                res[f'A|{J}|{th}|{X}|{rn}|{a}|{b}']=[q,e,h,eh,x,ex,bool(ok)]
# B. summary bins and N/S
BIN={'core150':[(-150,150)],'core100':[(-100,100)],'core50':[(-50,50)],'wings100_350':[(-350,-100),(100,350)],'blue':[(-350,-100)],'red':[(100,350)]}
def ratio_multi(J,th,rn,X,segs,vc=0,**kw):
    x=h=ex=eh=0
    for a,b in segs:
        i,e=binI(J,X,th,rn,a,b,vc,**kw); x+=i; ex=np.hypot(ex,e); i,e=binI(J,'HCO+',th,rn,a,b,vc,**kw); h+=i; eh=np.hypot(eh,e)
    q=x/h; return q,abs(q)*np.hypot(ex/x,eh/h),h,eh
def ns(J,th,X,segs,**kw):
    r={rn:ratio_multi(J,th,rn,X,segs,**kw) for rn in 'NS'}
    n=r['N'][0]/r['S'][0]; return r,n,n*np.hypot(r['N'][1]/r['N'][0],r['S'][1]/r['S'][0])
P('\n== B. common-grid summary bins (centre on nominal v_sys)')
for J,th in (('10',2.5),('32',2.5),('32',1.4)):
    for X in ['HCN','HNC']+(['CN_G2'] if J=='10' else []):
        for bn,segs in BIN.items():
            r,n,en=ns(J,th,X,segs)
            sg=min(r['N'][2]/r['N'][3],r['S'][2]/r['S'][3])
            P(f'J{J} {th}" {X:5s}/HCO+ {bn:12s} N={r["N"][0]:.3f}+-{r["N"][1]:.3f} S={r["S"][0]:.3f}+-{r["S"][1]:.3f} N/S={n:.2f}+-{en:.2f} HCO+ N={r["N"][2]:.1f} S={r["S"][2]:.1f} min SNR(HCO+)={sg:.1f}')
            res[f'B|{J}|{th}|{X}|{bn}']=[r['N'][0],r['S'][0],n,en]
# C. sensitivity of core100 N/S HCN/HCO+ (and HNC, CN)
P('\n== C. sensitivity of N/S (J=1-0, 2.5", core |v|<100 unless noted)')
for X in ['HCN','HNC','CN_G2']:
    vals=[]
    for lab,kw,segs in [('nominal',{},[(-100,100)]),('vsys+20',{'vc':20},[(-100,100)]),('vsys-20',{'vc':-20},[(-100,100)]),
        ('edges+-50 (150)',{},[(-150,150)]),('edges-50 (50)',{},[(-50,50)]),('edges 50..150 both',{},[(-150,-50),(50,150)]),
        ('baseline order2',{'order':2},[(-100,100)]),('baseline 400-1000',{'rng':(400,1000)},[(-100,100)]),('baseline 800-1500',{'rng':(800,1500)},[(-100,100)]),
        ('baseline order0',{'order':0},[(-100,100)])]:
        try:
            r,n,en=ns('10',2.5,X,segs,**kw); vals.append(n); P(f'  {X:5s} {lab:20s} N={r["N"][0]:.3f} S={r["S"][0]:.3f} N/S={n:.3f}+-{en:.3f}')
        except Exception as ex_: P('  ',X,lab,'ERR',ex_)
    P(f'  {X}: N/S range {min(vals):.2f}-{max(vals):.2f}'); res[f'C|{X}']=[min(vals),max(vals)]
P('  wings 100-350 (HCN) sensitivity:')
for lab,kw in [('nominal',{}),('vsys+20',{'vc':20}),('vsys-20',{'vc':-20}),('order2',{'order':2}),('base 400-1000',{'rng':(400,1000)}),('base 800-1500',{'rng':(800,1500)}),('order0',{'order':0})]:
    for X in ['HCN','CN_G2']:
        r,n,en=ns('10',2.5,X,BIN['wings100_350'],**kw); P(f'  {X:5s} {lab:14s} N={r["N"][0]:.2f}+-{r["N"][1]:.2f} S={r["S"][0]:.2f}+-{r["S"][1]:.2f} N/S={n:.2f}+-{en:.2f}  HCO+ N={r["N"][2]:.1f}')
# D. null test: leave-one-bin-out wing-width bins in baseline region
P('\n== D. null test: bins of 250 km/s width at line-free velocities (baseline fit with the null bin masked), J=1-0 2.5"')
for rn in 'NS':
    for X in ['HCN','HCO+','CN_G2']:
        vals=[]
        for c in [-1000,-900,-800,-700,700,800,900,1000]:
            if X=='CN_G2' and False: pass
            v,s,sig,dv,free=prof('10',X,2.5,rn,skip=(c-125,c+125)); I=integ(v,s,dv,c-125,c+125); vals.append(I)
        form=sig*np.sqrt(250*np.mean(dv[abs(v)<300]))
        P(f'  {rn} {X:6s} null-bin integrals (K km/s): rms={np.std(vals):.2f} mean={np.mean(vals):.2f} formal={form:.2f}  ratio rms/formal={np.std(vals)/form:.1f}')
        res[f'D|{rn}|{X}']=[float(np.std(vals)),float(form)]
# E. joint decomposition: 3 Gaussians, shared kinematics across lines; component 3 = broad (sigma>=150 km/s, FWHM>=350) 'outflow' candidate
from math import erf
def G(v,a,c,s): return a*np.exp(-0.5*((v-c)/s)**2)
def joint(J,th,lines,rn,lim=700):
    D={}
    for X in lines:
        v,s,sig,dv,_=prof(J,X,th,rn); w=abs(v)<lim; D[X]=(v[w],s[w],sig)
    n=len(lines)
    def mod(p,X,i):
        v=D[X][0]; a=p[9+3*i:12+3*i]; return sum(G(v,a[k],p[3*k],p[3*k+1]) for k in range(3))
    def f(p): return np.concatenate([(mod(p,X,i)-D[X][1])/D[X][2] for i,X in enumerate(lines)])
    lo=[-80,15,0, -80,45,0, -250,150,0]+[]; 
    lo=[-80,15,0,-80,45,0,-250,150,0]; hi=[80,45,0,80,110,0,250,400,0]
    lo=[-80,15,-80,45,-250,150]; hi=[80,45,80,110,250,400]
    # parameter layout: centres/sigmas interleaved as (c,s)*3 -> 6 values; amplitudes follow
    def mod(p,X,i):
        v=D[X][0]; a=p[6+3*i:9+3*i]; return sum(G(v,a[k],p[2*k],p[2*k+1]) for k in range(3))
    def f(p): return np.concatenate([(mod(p,X,i)-D[X][1])/D[X][2] for i,X in enumerate(lines)])
    best=None
    for p0c in ([0,30,0,70,50,200],[-20,25,20,80,0,300],[0,40,0,90,100,250]):
        p0=p0c+[0.3]*(3*n); r=least_squares(f,p0,bounds=(lo+[0]*(3*n),hi+[1e3]*(3*n)))
        if best is None or r.cost<best.cost: best=r
    r=best; cov=np.linalg.pinv(r.jac.T@r.jac); e=np.sqrt(np.diag(cov)); chi=2*r.cost/(len(r.fun)-len(r.x))
    return r.x,e,chi
def frac(c,s,a=-100,b=100): return 0.5*(erf((b-c)/(s*2**.5))-erf((a-c)/(s*2**.5)))
LN=['HCN','HCO+','HNC']
P('\n== E. joint 3-Gaussian fit (narrow sig15-45, mid sig45-110, broad sig150-400) shared kinematics HCN,HCO+,HNC, J=1-0 2.5"')
dec={}
for rn in 'NS':
    p,e,chi=joint('10',2.5,LN,rn)
    P(f' {rn}: red.chi2={chi:.1f}; comps (centre,sigma km/s): '+'; '.join(f'({p[2*k]:.0f},{p[2*k+1]:.0f})' for k in range(3)))
    F={X:[p[6+3*i+k]*np.sqrt(2*np.pi)*p[2*k+1] for k in range(3)] for i,X in enumerate(LN)}
    fr=[frac(p[2*k],p[2*k+1]) for k in range(3)]
    for X in LN: P(f'    {X:5s} flux narrow/mid/broad = '+'/'.join(f'{x:.1f}' for x in F[X])+' K km/s')
    P(f'    broad-comp fraction of HCO+ flux: {F["HCO+"][2]/sum(F["HCO+"]):.2f}; of it inside |v|<100: {fr[2]:.2f} -> {F["HCO+"][2]*fr[2]:.1f} K km/s of the core integral')
    P('    ratios narrow+mid (core) HCN/HCO+=%.3f HNC/HCO+=%.3f ; broad HCN/HCO+=%.2f'%(sum(F['HCN'][:2])/sum(F['HCO+'][:2]),sum(F['HNC'][:2])/sum(F['HCO+'][:2]),F['HCN'][2]/F['HCO+'][2] if F['HCO+'][2]>0.5 else float('nan')))
    dec[rn]=dict(F=F,fr=fr,p=p)
    res[f'E|{rn}']=dict(F=F,fr=fr,chi=chi)
cr={X:sum(dec['N']['F'][X][:2])/sum(dec['N']['F']['HCO+'][:2])/(sum(dec['S']['F'][X][:2])/sum(dec['S']['F']['HCO+'][:2])) for X in ['HCN','HNC']}
P(f' core(narrow+mid) N/S: HCN/HCO+ {cr["HCN"]:.2f}; HNC/HCO+ {cr["HNC"]:.2f}')
P('\n== E2. brackets: remove the broad-component HCO+ inside |v|<100 from the 1-0 core integral, assuming outflow HCN/HCO+ = rho')
core={rn:{X:binI('10',X,2.5,rn,-100,100)[0] for X in LN} for rn in 'NS'}
for rho in (0.0,0.43,0.6,1.0):
    q={}
    for rn in 'NS':
        Hb=dec[rn]['F']['HCO+'][2]*dec[rn]['fr'][2]; q[rn]=(core[rn]['HCN']-rho*Hb)/(core[rn]['HCO+']-Hb)
    P(f'  rho={rho:.2f}: HCN/HCO+ N={q["N"]:.3f} S={q["S"]:.3f} N/S={q["N"]/q["S"]:.2f}  (uncorrected {core["N"]["HCN"]/core["N"]["HCO+"]/(core["S"]["HCN"]/core["S"]["HCO+"]):.2f})'); res[f'E2|{rho}']=q['N']/q['S']
nh=core['S']['HCN']/0.88; P(f'  HCO+ at S that would have to be removed from the core (|v|<100, {core["S"]["HCO+"]:.1f} K km/s) to equalise HCN/HCO+ to the N value 0.88 if it carries no HCN: {core["S"]["HCO+"]-nh:.1f} K km/s ({100*(core["S"]["HCO+"]-nh)/core["S"]["HCO+"]:.0f}% of core); HCO+ in |v|>100 at S total: {19.2:.1f}')
# F. opacity/excitation
P('\n== F. 3-2/1-0 at matched 2.5" (K km/s ratio, core |v|<100) and CN hyperfine')
for rn in 'NS':
    for X in ['HCN','HCO+','HNC']:
        a,ea=binI('32',X,2.5,rn,-100,100); b,eb=binI('10',X,2.5,rn,-100,100)
        P(f'  {rn} {X:5s} I32={a:.1f}+-{ea:.1f} I10={b:.1f}+-{eb:.1f} R32/10={a/b:.3f}+-{abs(a/b)*np.hypot(ea/a,eb/b):.3f}'); res[f'F|{rn}|{X}']=a/b
for th in (2.5,3.5):
  for rn in 'NS':
    for (a,b) in [(-100,100),(-50,50),(-150,150)]:
        g1,e1=binI('10','CN_G1',th,rn,a,b); g2,e2=binI('10','CN_G2',th,rn,a,b)
        P(f'  CN 1-0 {th}" {rn} |v| {a}..{b}: G1(113.12-113.19 grp)/G2(113.49 grp)={g1/g2:.3f}+-{abs(g1/g2)*np.hypot(e1/g1,e2/g2):.3f}  (thin LTE 0.50; tau-> ratio<0.5 needs G1 weaker, thick->1)')
        res[f'F|CN|{th}|{rn}|{b}']=g1/g2
json.dump(res,open('/workspace/out/q53/s332.json','w'),indent=1,default=float)
open('/workspace/out/q53/s332.txt','w').write('\n'.join(out))
