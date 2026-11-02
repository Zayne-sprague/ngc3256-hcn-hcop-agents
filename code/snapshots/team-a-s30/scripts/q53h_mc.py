"""Q53 r3b: error budget + Monte Carlo class probabilities for HCN/HCO+ at N and S (4-3 at 1.0"; 4-3 core extrapolation; 3-2 core). Inputs: out/q53/fixedwin_all.csv, fixedwin43.csv,
posoff43.csv, spec43.npz. Output out/q53/mc.json"""
import numpy as np, pandas as pd, json
O='/workspace/out/q53/'; A=pd.read_csv(O+'fixedwin_all.csv'); F=pd.read_csv(O+'fixedwin43.csv'); PO=pd.read_csv(O+'posoff43.csv'); d=np.load(O+'spec43.npz')
g=lambda J,beam,reg,win='all',ous=None: A[(A.J==J)&(abs(A.beam-beam)<1e-6)&(A.reg==reg)&(A.win==win)&((A.ous==ous) if ous else True)].iloc[0]
res={}
# (1) window sensitivity of 4-3 at 1.0" (half-width HW, baseline from 450..950 fixed)
VC={'N':2785.,'S':2750.}
def ratio_hw(reg,hw,bl=(450,950)):
    out=[]
    for ln in ('HCN','HCO+'):
        v=d[f'{ln}_{reg}_1.0_v']; y=d[f'{ln}_{reg}_1.0_y']; vc=VC[reg]
        b=(abs(v-vc)>bl[0])&(abs(v-vc)<bl[1])&np.isfinite(y); c=np.polyfit(v[b]-vc,y[b],1); w=abs(v-vc)<hw
        out.append(np.sum(y[w]-np.polyval(c,v[w]-vc))*abs(np.median(np.diff(v))))
    return out[0]/out[1]
res['window_sens_1.0']={reg:{f'HW{hw}_BL{bl[0]}':round(ratio_hw(reg,hw,bl),3) for hw in (250,300,350,400,450) for bl in ((450,950),(500,900)) if bl[0]>hw} for reg in 'NS'}
# (2) pointing: rms of ratio over |offset|<=0.2" grid and full +-0.4"
res['pointing']={}
for reg in 'NS':
    t=PO[PO.reg==reg]; c=t[(abs(t.dx)<=0.21)&(abs(t.dy)<=0.21)]; res['pointing'][reg]={'rms_0.2':round(c.ratio.std(),3),'rms_0.4':round(t.ratio.std(),3),'center':round(t[(t.dx==0)&(t.dy==0)].ratio.iloc[0],3)}
# (3) uv / missing flux: smooth additive DeltaT_B (Xb47 - new cubes at 4.5") applied at 1.0"
Fr=lambda s,reg: F[(F.beam==4.5)&(F.reg==reg)&(F.set==s)].iloc[0]
uv={}
for reg in 'NS':
    n=Fr('HCN/HCO+',reg); x=Fr('HCN47/HCO47',reg); a=g(4,1.0,reg)
    dH=x.IHCN-n.IHCN; dO=x.IHCO-n.IHCO; uv[reg]=dict(dHCN=dH,dHCO=dO,frac_HCN=n.IHCN/x.IHCN,frac_HCO=n.IHCO/x.IHCO,ratio45_new=n.ratio,ratio45_Xb47=x.ratio,
        ratio10_raw=a.ratio,ratio10_corr=(a.IHCN+dH)/(a.IHCO+dO))
res['uv_test']={k:{kk:round(float(vv),3) for kk,vv in v.items()} for k,v in uv.items()}
# (4) R43/R32 vs beam
rr={}
for reg in 'NS':
    bs=[1.0,1.4,1.9,2.5]; y=[]; e=[]
    for b in bs:
        a=g(4,b,reg); o='Xb4b' if b<=1.9 else 'Xb4d'; c=g(3,b,reg,ous=o); y.append(a.ratio/c.ratio); e.append(y[-1]*np.hypot(a.err/a.ratio,c.err/c.ratio))
    p=np.polyfit(bs,y,1,w=1/np.array(e)); rr[reg]=dict(beams=bs,R43_R32=[round(v,3) for v in y],slope_per_as=round(p[0],3),at1=round(np.polyval(p,1.0),3),at03=round(np.polyval(p,0.3),3))
res['R43_over_R32']=rr
# (5) Monte Carlo
rng=np.random.default_rng(53); n=200000; out={}
for reg in 'NS':
    a=g(4,1.0,reg); u=uv[reg]
    r10=a.ratio*rng.normal(1,a.err/a.ratio,n)                                  # statistical
    r10*=np.exp(rng.normal(0,np.hypot(.10,.10),n))                              # per-project flux scale 10% each (two projects)
    cf=u['ratio10_corr']/u['ratio10_raw']; r10*=rng.normal(1+(cf-1)/2,abs(cf-1)/2+1e-6,n)  # missing-flux correction: half-way, +-half
    r10*=rng.normal(1,res['pointing'][reg]['rms_0.2']/res['pointing'][reg]['center'],n)     # pointing 0.2"
    hw=list(res['window_sens_1.0'][reg].values()); r10*=rng.normal(1,np.std(hw)/np.mean(hw),n)  # window choice
    # core: 3-2 at 0.3" x R43/R32 extrapolation
    c32=g(3,0.3,reg,ous='Xb4b'); r32=c32.ratio*rng.normal(1,c32.err/c32.ratio,n)*np.exp(rng.normal(0,.05,n))*rng.normal(1,.05,n)  # stat; 5% 3-2 relative scale; 5% pointing at 0.3"
    f=rr[reg]; rat1=g(4,1.0,reg).ratio/g(3,1.0,reg,ous='Xb4b').ratio
    ex=rng.uniform(1.0,f['at03']/f['at1'],n)                                    # flat (1) ... linear-fit extrapolation
    core=r32*rat1*ex*np.exp(rng.normal(0,np.hypot(.10,.10),n))*rng.normal(1+(cf-1)/2,abs(cf-1)/2+1e-6,n)
    for nm,x in (('4-3 @1.0" (188 pc)',r10),('4-3 core 0.3" (56 pc), flat..linear R43/R32',core)):
        q=np.percentile(x,[16,50,84]); out[f'{reg} {nm}']=dict(p16_50_84=[round(v,3) for v in q],P_le_1p03=round((x<=1.03).mean(),3),P_ge_1p11=round((x>=1.11).mean(),3),P_gap=round(((x>1.03)&(x<1.11)).mean(),3),P_lt_0p48=round((x<=0.48).mean(),3),P_AGN_median_2p58_excl=round((x<1.11).mean(),3))
res['MC']=out
json.dump(res,open(O+'mc.json','w'),indent=1); print(json.dumps(res,indent=1))
