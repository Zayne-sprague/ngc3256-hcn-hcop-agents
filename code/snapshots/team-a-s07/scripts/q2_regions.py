# Q2 reviewer fix 3: ratio-based (n,T) constraints for every region we can separate (N, S, Bridge, outer disk, N wings, S-outflow blue lobe), with the honest degeneracy statement.
# (1) per-molecule ladder (HCN, HCO+ 3-2/2-1 and 4-3/3-2, Q1b Xb4f deblend, error = hypot(14%, deblend spread)) -> chi2 on the RADEX grid (Tbg 2.73 K, dv 150) over (T,n,N); marginal 16-50-84% of T and n (flat priors, T=10-150 K, N<=n+12.5)
#     and the cells within dchi2<4 of the minimum ('allowed band'); 2 data, 3 parameters => band, not a solution.
# (2) R32 = HCN/HCO+ 3-2 at matched beam (Xb4d 1.05"; Q1 ratios, Q3 S-outflow, Q9) with common n, T and column ratio X=N(HCN)/N(HCO+) in {0.3,1,3,10}, both dv=150: cells with |R_model-R|<2 sigma, range of n and T.
# Output out/q2/q2_regions_ladder.csv, q2_regions_R32.csv
import numpy as np, pandas as pd, sys
sys.path.insert(0,'/workspace/scripts')
O='/workspace/out/'; G={m:np.load(O+f'q4c/grid4_{m}.npz') for m in ('hcn','hco+')}
TK,LN,LNC=G['hcn']['TK'],G['hcn']['LN'],G['hcn']['LNC']
def wq(x,w,q): i=np.argsort(x); c=np.cumsum(w[i]); c/=c[-1]; return float(np.interp(q,c,x[i]))
TT,NN,CC=np.meshgrid(TK,LN,LNC,indexing='ij')
phys=(CC<=NN+12.5)
L=pd.read_csv(O+'q1/q1_deblend_ladders_xb4f.csv')
rows=[]
for comp in ('N','S','Bridge','Outer'):
    for sp,k in (('HCN','hcn'),('HCO+','hco+')):
        TR=G[k]['TR']; tau=G[k]['tau'][...,2]
        d=[];m=[]
        for q,(i,j) in (('3-2 [Xb4f]/2-1 (Tb)',(2,1)),('4-3/3-2 [Xb4f 3-2] (Tb)',(3,2))):
            r=L[(L.comp==comp)&(L.species==sp)&(L.quantity==q)].iloc[0]; d.append((r.value,np.hypot(0.14*r.value,(r.spread_hi-r.spread_lo)/2))); m.append(TR[...,i]/TR[...,j])
        ok=phys&np.isfinite(m[0])&np.isfinite(m[1])
        chi=sum(((m[i]-d[i][0])/d[i][1])**2 for i in range(2)); chi=np.where(ok,chi,1e9)
        for sel,mask in (('all columns',np.ones_like(ok)),('tau32>1',tau>1)):
            c2=np.where(mask,chi,1e9); cm=c2.min(); w=np.exp(-0.5*(c2-cm)).ravel()
            band=c2<=cm+4
            r=dict(comp=comp,species=sp,selection=sel,d_32_21=d[0][0],e_32_21=d[0][1],d_43_32=d[1][0],e_43_32=d[1][1],chi2min=cm)
            for nm,a in (('T',TT),('logn',NN),('logN',CC)):
                for p in (16,50,84): r[f'{nm}_{p}']=wq(a.ravel(),w,p/100)
                r[nm+'_band_lo']=float(a[band].min()); r[nm+'_band_hi']=float(a[band].max())
            r['frac_grid_in_band']=float(band.sum()/ok.sum()); rows.append(r)
pd.DataFrame(rows).round(3).to_csv(O+'q2/q2_regions_ladder.csv',index=False)
# (2) R32
REG=[('N core 1.05" (Xb4d)',0.840,0.042),('S core 1.05" (Xb4d)',0.408,0.021),('Bridge 1.05" (Xb4d)',0.467,0.025),('N wings 100-260 km/s 1.05" (Xb4d)',1.35,0.32),
     ('S-outflow blue lobe (Q3, 2475-2575 km/s, Xb4d)',1.2,0.3),('outer disk 6.5" (ring zero level; 0.06-0.27 spread)',0.17,0.11)]
TH,TO=G['hcn']['TR'][...,2],G['hco+']['TR'][...,2]; tH,tO=G['hcn']['tau'][...,2],G['hco+']['tau'][...,2]
rr=[]
for nm,R,e in REG:
    for lx in (-0.5,0,0.5,1.0):
        s=int(round(lx/0.25)); 
        # N(HCN)=X N(HCO+): index of HCN column = index of HCO+ column + s
        a=np.full_like(TH,np.nan)
        if s>=0: a[...,s:]=TH[...,:TH.shape[-1]-s]; ta=np.full_like(TH,np.nan); ta[...,s:]=tH[...,:TH.shape[-1]-s]
        else: a[...,:s]=TH[...,-s:]; ta=np.full_like(TH,np.nan); ta[...,:s]=tH[...,-s:]
        Rm=a/TO; ok=np.isfinite(Rm)&(np.abs(Rm-R)<2*e)&phys&(tO>1)
        r=dict(region=nm,R=R,eR=e,logX=lx,n_cells=int(ok.sum()),frac_of_valid=float(ok.sum()/np.isfinite(Rm).sum()))
        if ok.sum():
            r.update(logn_min=float(NN[ok].min()),logn_max=float(NN[ok].max()),logn_med=float(np.median(NN[ok])),T_min=float(TT[ok].min()),T_max=float(TT[ok].max()),T_med=float(np.median(TT[ok])),
                     frac_T_le_50=float((TT[ok]<=50).mean()),frac_logn_ge5=float((NN[ok]>=5).mean()))
        rr.append(r)
pd.DataFrame(rr).round(3).to_csv(O+'q2/q2_regions_R32.csv',index=False)
pd.set_option('display.width',250,'display.max_columns',40); print(pd.DataFrame(rows).round(2).to_string()); print(pd.DataFrame(rr).round(2).to_string())
