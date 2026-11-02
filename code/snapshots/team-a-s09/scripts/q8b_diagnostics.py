# Q8b: recovery-robust diagnostics and total error budget of the N/S T43/T32 contrast (reviewer round 1).
#  A  HCN-vs-HCO+ double ratio of the N/S contrast, per beam, per 3-2 array (cancels common recovery/pointing/calibration)
#  B  three-step ladder 2-1/3-2/4-3 at the Xb55 beam 9.5x6.4" (N,S,Bridge), N/S per step, no size model
#  C  line/continuum double ratios (continuum Xb47/Xb4f at the same positions/beams)
#  D  empirical null: T43/T32 (Xb47 4-3 / Xb4f 3-2 at the Xb4f beam) over a grid of control positions; |ln(R1/R2)| of position pairs
#     at the N-S separation (5.1"), same S/N cut -> percentile of the observed N/S
#  E  Monte Carlo error budget of N/S (beam in {Xb4f 6.2x3.7, circ 6.5, Xb55 9.5x6.4}, pos/peak, N and S position jitter 0.3", Xb47-vs-Xb4f
#     astrometric shift 0.3" per axis, thermal noise) + the 16 deblend variants of Q1 (3-component fits: N,S,Bridge,Outer) with Xb4f 3-2
import sys; sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd
from scipy.ndimage import map_coordinates
O='/workspace/out/q8b/'
rng=np.random.default_rng(81)
M=pd.read_csv('/workspace/out/q8/matched_m0.csv'); Cn=pd.read_csv('/workspace/out/q8/T43_T32_contrast.csv')
out={}
def NSerr(a,b): return np.hypot(a[1]/a[0],b[1]/b[0])
# ---------- A: double ratio HCN/HCO+ of N/S contrast
rows=[]
for beam in Cn.beam.unique():
  for arr in Cn.arr32.unique():
    for meas in ('pos_','peak'):
      g=Cn[(Cn.beam==beam)&(Cn.arr32==arr)&(Cn.measure==meas)&(Cn.region=='N/S')]
      h=g[g.line=='HCN']; c=g[g.line=='HCO+']
      if len(h)==1 and len(c)==1:
        a,b=h.iloc[0],c.iloc[0]; d=a.T43_T32/b.T43_T32; e=d*np.hypot(a.err/a.T43_T32,b.err/b.T43_T32)
        rows.append(dict(beam=beam,arr32=arr,measure=meas,NS_HCN=a.T43_T32,NS_HCOp=b.T43_T32,double_ratio=d,err_stat=e))
A=pd.DataFrame(rows); A.to_csv(O+'double_ratio_HCN_HCOp.csv',index=False); print('A: double ratio HCN/HCO+ of N/S'); print(A.round(3).to_string())
# ---------- B: ladder at the Xb55 beam
B=[]; bm='Xb55(9.5x6.4,-88)'
def T(sp,tr,reg,col='Tb_pos_Kkms'):
    r=M[(M.line==sp)&(M.transition==tr)&(M.beam==bm)&(M.region==reg)].iloc[0]; return r[col],r.eTb
for sp in ('HCN','HCO+'):
  for a32 in ('3-2 Xb4f','3-2 Xb4d','3-2 Xb4b'):
    d={}
    for reg in ('N','S','Bridge'):
        t21,e21=T(sp,'2-1 Xb55',reg); t32,e32=T(sp,a32,reg); t43,e43=T(sp,'4-3 Xb47',reg)
        r1=t32/t21; r2=t43/t32; d[reg]=(r1,r1*np.hypot(e21/t21,e32/t32),r2,r2*np.hypot(e32/t32,e43/t43))
        B.append(dict(line=sp,arr32=a32,region=reg,T21=t21,T32=t32,T43=t43,T32_T21=r1,e1=d[reg][1],T43_T32=r2,e2=d[reg][3]))
    for (a,b) in (('N','S'),('N','Bridge'),('S','Bridge')):
        B.append(dict(line=sp,arr32=a32,region=f'{a}/{b}',T32_T21=d[a][0]/d[b][0],e1=d[a][0]/d[b][0]*np.hypot(d[a][1]/d[a][0],d[b][1]/d[b][0]),T43_T32=d[a][2]/d[b][2],e2=d[a][2]/d[b][2]*np.hypot(d[a][3]/d[a][2],d[b][3]/d[b][2])))
B=pd.DataFrame(B); B.to_csv(O+'ladder_Xb55beam.csv',index=False); print('B: ladder at Xb55 beam'); print(B.round(3).to_string())
# ---------- C: line/continuum double ratio at the Xb4f beam (continuum Xb47/Xb4f at positions)
cont=pd.read_csv('/workspace/out/q8/continuum_matched.csv'); cb='Xb4f(6.2x3.7,82)'
pc=cont[cont.beam==cb].pivot_table(index='array',columns='region',values='mJy_beam')
cNS=(pc.loc['Xb47','N']/pc.loc['Xb4f','N'])/(pc.loc['Xb47','S']/pc.loc['Xb4f','S'])
rows=[]
for sp in ('HCN','HCO+'):
    g=Cn[(Cn.line==sp)&(Cn.beam==cb)&(Cn.arr32=='3-2 Xb4f')&(Cn.measure=='pos_')&(Cn.region=='N/S')].iloc[0]
    rows.append(dict(line=sp,NS_line=g.T43_T32,err=g.err,NS_continuum_Xb47_over_Xb4f=cNS,line_over_cont=g.T43_T32/cNS,err_stat=g.err/cNS))
C=pd.DataFrame(rows); C.to_csv(O+'line_continuum_double_ratio.csv',index=False); print('C'); print(C.round(3).to_string())
# ---------- D/E need maps: Xb47 4-3 and Xb4f 3-2 moment-0 maps at three beams
BEAMS={'Xb4f(6.2x3.7,82)':(6.2,3.7,82.),'circ6.5':(6.5,6.5,0.),'Xb55(9.5x6.4,-88)':(9.5,6.4,-88.)}
def sm(d,tb): 
    k=kernel(d['hd'],tb,d['beam']); return fftconvolve(np.nan_to_num(d['m0']),k,mode='same')*OM(tb)/OM(d['beam'])
def noise(d,tb):
    s=[]
    for _,n in d['nulls']:
        k=kernel(d['hd'],tb,d['beam']); x=fftconvolve(np.nan_to_num(n),k,mode='same')*OM(tb)/OM(d['beam']); r=np.hypot(d['dx'],d['dy']); v=x[(r<12)&(x!=0)]; s.append(1.4826*np.median(np.abs(v-np.median(v))))
    return np.percentile(s,25)
FRQ={'HCN':(265.886434,354.505478),'HCO+':(267.557626,356.734223)}
MAP={}
for sp in ('HCN','HCO+'):
    d43=load(f'{sp} 4-3'); d32=load(f'{sp} 3-2@Xb4f')
    for bn,tb in BEAMS.items():
        MAP[(sp,bn)]=dict(a=sm(d43,tb),b=sm(d32,tb),d43=d43,d32=d32,s43=noise(d43,tb),s32=noise(d32,tb),fac=(FRQ[sp][0]/FRQ[sp][1])**2)  # Tb(4-3)/Tb(3-2) = (S43/S32)*(nu32/nu43)^2
def R(sp,bn,x43,y43,x32,y32,peak=False,noisy=True):
    m=MAP[(sp,bn)]
    if peak:
        for img,d,x,y,k in ((m['a'],m['d43'],x43,y43,'p4'),(m['b'],m['d32'],x32,y32,'p3')):
            mm=np.hypot(d['dx']-x,d['dy']-y)<0.8; v=np.nanmax(np.where(mm,img,np.nan)); m[k]=v
        s4,s3=m['p4'],m['p3']
    else: s4=sample(m['a'],m['d43'],x43,y43); s3=sample(m['b'],m['d32'],x32,y32)
    if noisy: s4+=rng.normal(0,m['s43']); s3+=rng.normal(0,m['s32'])
    return s4/s3*m['fac']
POS={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15)}
# ---------- D: control positions
rows=[]
for sp in ('HCN','HCO+'):
  for bn in ('Xb4f(6.2x3.7,82)','circ6.5'):
    m=MAP[(sp,bn)]; d=m['d32']
    xs=np.arange(-13,13.01,1.4); pts=[]
    for x in xs:
      for y in xs:
        s3=sample(m['b'],m['d32'],x,y); s4=sample(m['a'],m['d43'],x,y)
        pts.append((x,y,s3,s4))
    P=np.array(pts); 
    for snr in (10,20):
        ok=(P[:,2]>snr*m['s32'])&(P[:,3]>snr*m['s43']); Q=P[ok]; r=Q[:,3]/Q[:,2]*m['fac']
        for lo,hi in ((4.6,5.6),):
            dx=Q[:,0][:,None]-Q[:,0][None]; dy=Q[:,1][:,None]-Q[:,1][None]; dd=np.hypot(dx,dy); i,j=np.where((dd>lo)&(dd<hi)&(np.arange(len(Q))[:,None]<np.arange(len(Q))[None]))
            if len(i)==0: continue
            lr=np.abs(np.log(r[i]/r[j]))
            # restrict also to roughly N-S oriented pairs (within 35 deg of the N-S axis)
            ang=np.degrees(np.arctan2(np.abs(dx[i,j]),np.abs(dy[i,j]))); ns=ang<35
            obs=np.log(1.333 if sp=='HCN' else 1.133)
            rows.append(dict(line=sp,beam=bn,snr_cut=snr,npos=int(ok.sum()),npairs=len(i),med_absln=np.median(lr),p84=np.percentile(lr,84),p95=np.percentile(lr,95),frac_ge_obsNS=float(np.mean(lr>=abs(obs))),
                 pairs_NSoriented=int(ns.sum()),frac_ge_obsNS_nsor=float(np.mean(lr[ns]>=abs(obs))) if ns.sum()>0 else np.nan,R_range=f'{r.min():.2f}-{r.max():.2f}',R_median=np.median(r)))
D=pd.DataFrame(rows); D.to_csv(O+'control_positions_null.csv',index=False); print('D'); print(D.round(3).to_string())
# double ratio null: HCN contrast/HCO+ contrast for the same pairs
rows=[]
for bn in ('Xb4f(6.2x3.7,82)',):
    mh,mc=MAP[('HCN',bn)],MAP[('HCO+',bn)]; xs=np.arange(-13,13.01,1.4); pts=[]
    for x in xs:
      for y in xs:
        pts.append((x,y,sample(mh['b'],mh['d32'],x,y),sample(mh['a'],mh['d43'],x,y),sample(mc['b'],mc['d32'],x,y),sample(mc['a'],mc['d43'],x,y)))
    P=np.array(pts)
    for snr in (10,20):
        ok=(P[:,2]>snr*mh['s32'])&(P[:,3]>snr*mh['s43'])&(P[:,4]>snr*mc['s32'])&(P[:,5]>snr*mc['s43']); Q=P[ok]
        rh=Q[:,3]/Q[:,2]*mh['fac']; rc=Q[:,5]/Q[:,4]*mc['fac']
        dd=np.hypot(Q[:,0][:,None]-Q[:,0][None],Q[:,1][:,None]-Q[:,1][None]); i,j=np.where((dd>4.6)&(dd<5.6)&(np.arange(len(Q))[:,None]<np.arange(len(Q))[None]))
        dr=np.abs(np.log((rh[i]/rh[j])/(rc[i]/rc[j]))); obs=np.log(1.333/1.133)
        rows.append(dict(beam=bn,snr_cut=snr,npos=int(ok.sum()),npairs=len(i),med_absln_double=np.median(dr),p84=np.percentile(dr,84),p95=np.percentile(dr,95),frac_ge_obs=float(np.mean(dr>=obs)),obs_double=1.333/1.133))
Dd=pd.DataFrame(rows); Dd.to_csv(O+'control_positions_double_ratio_null.csv',index=False); print(Dd.round(3).to_string())
# ---------- E: Monte Carlo
NMC=4000; res={'HCN':[],'HCO+':[]}; recs=[]
for k in range(NMC):
    bn=list(BEAMS)[rng.integers(3)]; peak=rng.random()<0.5; sh=rng.normal(0,0.3,2)   # 4-3 map shifted by sh relative to 3-2
    jn=rng.normal(0,0.3,4)                                                          # jitter of the N and S sample positions (shared by both maps)
    row={}
    for sp in ('HCN','HCO+'):
        rN=R(sp,bn,POS['N'][0]+jn[0]+sh[0],POS['N'][1]+jn[1]+sh[1],POS['N'][0]+jn[0],POS['N'][1]+jn[1],peak)
        rS=R(sp,bn,POS['S'][0]+jn[2]+sh[0],POS['S'][1]+jn[3]+sh[1],POS['S'][0]+jn[2],POS['S'][1]+jn[3],peak)
        row[sp]=rN/rS
    recs.append(dict(beam=bn,peak=peak,NS_HCN=row['HCN'],NS_HCOp=row['HCO+'],double=row['HCN']/row['HCO+']))
E=pd.DataFrame(recs); E.to_csv(O+'mc_NS_samples.csv',index=False)
pc=lambda x: np.percentile(x,[2.5,16,50,84,97.5])
rows=[]
for col in ('NS_HCN','NS_HCOp','double'):
    p=pc(E[col]); rows.append(dict(quantity=col,sel='all',p2_5=p[0],p16=p[1],p50=p[2],p84=p[3],p97_5=p[4],P_gt1=float(np.mean(E[col]>1)),sigma_MC=E[col].std()))
    for bn in BEAMS:
        p=pc(E[E.beam==bn][col]); rows.append(dict(quantity=col,sel=bn,p2_5=p[0],p16=p[1],p50=p[2],p84=p[3],p97_5=p[4],P_gt1=float(np.mean(E[E.beam==bn][col]>1)),sigma_MC=E[E.beam==bn][col].std()))
    for pk in (False,True):
        p=pc(E[E.peak==pk][col]); rows.append(dict(quantity=col,sel=f'peak={pk}',p2_5=p[0],p16=p[1],p50=p[2],p84=p[3],p97_5=p[4],P_gt1=float(np.mean(E[E.peak==pk][col]>1)),sigma_MC=E[E.peak==pk][col].std()))
Ms=pd.DataFrame(rows); Ms.to_csv(O+'mc_NS_summary.csv',index=False); print('E: MC'); print(Ms.round(3).to_string())
# ---------- E2: 16 deblend variants of Q1 (N,S,Bridge,Outer fitted simultaneously with 1" templates), Xb4f 3-2 and Xb47 4-3 (6.5" fit beam)
a=pd.read_csv('/workspace/out/q1/q1_deblend_amplitudes.csv'); tags=[t for t in a.variant.unique() if t.endswith('/nomiss')]
rows=[]
for sp in ('HCN','HCO+'):
  for t in tags:
    g=lambda l,c: a[(a.line==l)&(a.template==sp)&(a.variant==t)&(a.comp==c)].a
    try:
        v={c:(g(f'{sp} 4-3',c).iloc[0]/g(f'{sp} 3-2@Xb4f',c).iloc[0], g(f'{sp} 3-2@Xb4f',c).iloc[0]/g(f'{sp} 2-1',c).iloc[0]) for c in ('N','S','Bridge')}
    except Exception as e: continue
    rows.append(dict(line=sp,variant=t,T43_T32_N=v['N'][0],T43_T32_S=v['S'][0],T43_T32_B=v['Bridge'][0],NS=v['N'][0]/v['S'][0],T32_T21_N=v['N'][1],T32_T21_S=v['S'][1],NS_32_21=v['N'][1]/v['S'][1]))
V=pd.DataFrame(rows); V.to_csv(O+'deblend_variants_NS.csv',index=False)
print('E2'); print(V.round(3).to_string()); print(V.groupby('line')[['T43_T32_N','T43_T32_S','NS','NS_32_21']].agg(['median','min','max']).round(3))
