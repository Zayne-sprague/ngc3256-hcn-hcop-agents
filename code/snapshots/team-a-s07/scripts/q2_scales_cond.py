# Q2 reviewer fix 4: physical scales as explicit functions of the ASSUMED X(HCN)=N(HCN)/N(H2), per class (A,B,F2) and pooled, T floor removed (Tmin=10 K = grid edge), Tbg 2.73 K.
# N(HCN, cloud) = N(HCN, dv=150)/3 (per-cloud dv=50 km/s; RADEX depends on N/dv).  1 Msun/pc2 = 4.6e19 H2 cm-2 (incl. He); L = N(HCN)/(X n_H2) (coherent path through the thick-line column, filling 1);
# Sigma tension expressed as a LOWER LIMIT on X(HCN): X_min = N(HCN,cloud)/(Sigma_max 4.6e19) for Sigma_max = 470 (CO GMC median, Brunetti+2022), 1e3, 3e3, 1e4 Msun/pc2 (nuclear CO surface densities of LIRG nuclei span ~1e3-1e4;
# the value for the NGC 3256 nuclei at 0.3" is not measured here).  Output out/q2/q2_scales_cond.csv
import os,sys; os.environ['Q2_TMIN']='0'
sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q4c_joint as J
POOL={'A':'A both thick, kappa=1','B':'B both thick, same Tex (0.8-1.25), kappa free','F2':'F2 separate n, kappa=1, both thick'}
rows=[]
for c in ('N','S'):
    d,use,chi,cls,aux=J.cube(c,'all','env'); sh=aux['sh']; B=lambda a: np.broadcast_to(a,sh)
    nH=B(J.LN[None,:,None,None,None]); NH=B(J.LNC[None,None,None,:,None])-np.log10(3.)
    gbest=min(float(np.where(cls[n][0],chi[ik],1e9).min()) for n in POOL.values() for ik in cls[n][1]); Wc={}
    for key,nm in POOL.items():
        mask,iks=cls[nm]; w=sum(np.exp(-0.5*(np.where(mask,chi[ik],1e9)-gbest)) for ik in iks)/len(iks); Wc[key]=w
    Wc['pooled']=sum(Wc.values())/3
    for key,w in Wc.items():
        w=w.ravel(); s=w>1e-6; w=w[s]; N_=NH.ravel()[s]; n_=nH.ravel()[s]
        def add(q,v,**kw): rows.append(dict(region=c,cls=key,quantity=q,p16=J.wq(v,w,.16),p50=J.wq(v,w,.5),p84=J.wq(v,w,.84),**kw))
        add('log N(HCN,cloud dv50) [cm-2]',N_,logX_HCN=np.nan)
        for lx in (-9,-8,-7):
            NH2=N_-lx; add('log N(H2) [cm-2]',NH2,logX_HCN=lx); add('log Sigma_H2 [Msun/pc2]',NH2-np.log10(4.6e19),logX_HCN=lx)
            L=NH2-n_-np.log10(3.086e18); add('log L [pc]',L,logX_HCN=lx); add('L/(2R), R=88 pc (fill length fraction)',L-np.log10(176.),logX_HCN=lx)
        for S in (470,1e3,3e3,1e4): add('log X_min(HCN) for Sigma_max=%d Msun/pc2'%S,N_-np.log10(S*4.6e19),logX_HCN=np.nan)
R=pd.DataFrame(rows); R.round(3).to_csv('/workspace/out/q2/q2_scales_cond.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',400); print(R[(R.cls!='pooled')|True][R.quantity.str.contains('X_min|Sigma|cloud')].round(2).to_string())
