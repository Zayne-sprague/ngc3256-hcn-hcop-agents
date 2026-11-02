# Q13c: independent-transition check of the S-outflow blue-lobe HCN excess: HCN/HCO+ 4-3 (Xb47, no HDO line near HCN 4-3) vs 3-2 (Xb4d), both convolved to 4.5" circular.
# E_J = HCN_J(blue window) - k_J*HCO+_J(blue window); k_J from the core window (2575-2875) of the same position/beam. Noise = scatter of the same window in
# random apertures (>=4.5" from N,S,Soutflow,Bridge) after per-pixel linear baseline (line-free channels). Tb ratio r43 = (S43/S32)*(nu32/nu43)^2 (same beam, same window).
# Outputs: out/q13/q13c.log, q13c_table.csv, q13c_spectra.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
import pandas as pd
O='/workspace/out/q13/'
FW=4.5
L={ # name: (ous, spw, rest GHz, free velocity ranges (radio, own frame))
 'HCN32':('Xb4d','spw27',265.886434,[(1380,2100)]),
 'HCO32':('Xb4d','spw25',267.557626,[(1280,1680),(1880,2100)]),
 'HCN43':('Xb47','spw16',354.505477,None),
 'HCO43':('Xb47','spw18',356.734223,None)}
WIN={'blue':(2375,2575),'blueB':(2300,2575),'core':(2575,2875),'red':(2975,3175)}
def load(n):
    ous,spw,f0,free=L[n]; c=fcube(ous,spw,FW,16.); v=vel(c['fr'],f0); c['v']=v; c['dv']=abs(np.median(np.diff(v))); return c
def base(c,free):
    v=c['v']; m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    X=np.vander(v-v[m].mean(),2); A=c['a'][m].reshape(m.sum(),-1); coef=np.linalg.lstsq(X[m],A,rcond=None)[0]
    return c['a']-(X@coef).reshape(len(v),*c['a'].shape[1:]),m
if __name__=='__main__':
    log=open(O+'q13c.log','w')
    def pr(*a):
        s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
    D={}
    for n in L:
        c=load(n); v=c['v']; pr(f'{n}: v range {v.min():.0f}..{v.max():.0f}, dv {c["dv"]:.2f} km/s, beam {c["beam"]}')
        D[n]=c
    # line-free ranges for the 4-3 cubes: use the channels farthest from the line, both sides, excluding +-650 km/s around 2775 and known neighbours
    FREE43={'HCN43':[(1790,2100),(3100,3280)],'HCO43':[(1600,2100),(3100,3400)]}
    free={'HCN32':L['HCN32'][3],'HCO32':L['HCO32'][3],**FREE43}
    R={}
    for n,c in D.items():
        r,m=base(c,free[n]); R[n]=r; c['free']=m
        sig=1.4826*np.median(np.abs(r[m]-np.median(r[m],0)),0)    # per-pixel channel rms
        c['rms']=sig
    # spectra at positions (point sample mean r=1")
    pos=['Soutflow','S','N','Bridge']
    rows=[];spec=[]
    for p in pos:
        x,y=POS[p]
        for n,c in D.items():
            m=np.hypot(c['dx']-x,c['dy']-y)<=1.0; s=R[n][:,m].mean(1); v=c['v']
            for vv,ss in zip(v,s): spec.append(dict(pos=p,line=n,v=vv,I=ss))
    pd.DataFrame(spec).round(5).to_csv(O+'q13c_spectra.csv',index=False)
    # window maps
    W={}
    for n,c in D.items():
        v=c['v']
        for w,(a,b) in WIN.items():
            m=(v>=a)&(v<b); W[(n,w)]=(R[n][m].sum(0)*c['dv'],m.sum())
        pr(f'{n}: channels in blue window {W[(n,"blue")][1]}, in red {W[(n,"red")][1]}')
    rng=np.random.default_rng(5); ex=[POS[p] for p in pos]
    c0=D['HCN32'];
    P=[]
    while len(P)<200:
        x,y=rng.uniform(-9,9,2)
        if np.hypot(x,y)>9 or min(np.hypot(x-a,y-b) for a,b in ex)<4.5: continue
        P.append((x,y))
    def ap(n,w,x,y,rad=1.0):
        c=D[n]; m=np.hypot(c['dx']-x,c['dy']-y)<=rad; return W[(n,w)][0][m].mean()
    rob=lambda a: 1.4826*np.median(np.abs(a-np.median(a)))
    out=[]
    for p in pos:
        x,y=POS[p]
        for J,(h,o,nu_ratio) in {'3-2':('HCN32','HCO32',1.0),'4-3':('HCN43','HCO43',1.0)}.items():
            row=dict(pos=p,J=J)
            for w in WIN:
                row[f'HCN_{w}']=ap(h,w,x,y); row[f'HCO_{w}']=ap(o,w,x,y)
                row[f'sHCN_{w}']=rob(np.array([ap(h,w,*q) for q in P])); row[f'sHCO_{w}']=rob(np.array([ap(o,w,*q) for q in P]))
            k=row['HCN_core']/row['HCO_core']
            for w in ('blue','blueB','red'):
                Eo=row[f'HCN_{w}']-k*row[f'HCO_{w}']
                nul=np.array([ap(h,w,*q)-k*ap(o,w,*q) for q in P]); sE=rob(nul)
                row[f'E_{w}']=Eo; row[f'sE_{w}']=sE; row[f'snrE_{w}']=Eo/sE
            row['k_core']=k; out.append(row)
    T=pd.DataFrame(out); T.round(4).to_csv(O+'q13c_table.csv',index=False)
    pd.set_option('display.width',250); pd.set_option('display.max_columns',60)
    pr(T[['pos','J','k_core','HCN_core','HCN_blue','sHCN_blue','snrE_blue','E_blue','sE_blue','HCO_blue','sHCO_blue','E_blueB','sE_blueB','E_red','sE_red']].round(3).to_string(index=False))
    # Tb ratios 4-3/3-2 (same 4.5" beam): Tb ~ S/nu^2
    nu32=265.886434/1.0093; nu43=354.505477/1.0093
    for p in pos:
        a=T[(T.pos==p)&(T.J=='3-2')].iloc[0]; b=T[(T.pos==p)&(T.J=='4-3')].iloc[0]
        for line,key in (('HCN','HCN'),('HCO+','HCO')):
            for w in ('core','blue','blueB'):
                s32=a[f'{key}_{w}']; s43=b[f'{key}_{w}']; e32=a[f's{key}_{w}']; e43=b[f's{key}_{w}']
                r=s43/s32*(nu32/nu43)**2 if s32>0 else np.nan
                er=abs(r)*np.hypot(e32/s32,e43/s43) if (s32>0 and s43>0) else np.nan
                pr(f'{p:9s} {line:5s} {w:6s}: S32={s32:.3f}+-{e32:.3f} S43={s43:.3f}+-{e43:.3f} Jy/b km/s -> Tb r43/32={r:.2f}+-{er:.2f}')
