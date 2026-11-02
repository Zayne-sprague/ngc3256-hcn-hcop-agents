# Q26: observed R_ext(3-2)/R_ext(2-1) from the parametric fit with EB-jackknife + baseline-row bootstrap, profile variants.
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts'); import q26_uvana as U
rng=np.random.default_rng(1)
def Rext(L,fw,ul):
    return U.fit(L['HCN'],fw,ul)[0][2]/U.fit(L['HCOp'],fw,ul)[0][2]
def boot(L):
    out={}
    for sp in L:
        out[sp]=[]
        for d in L[sp]:
            i=rng.integers(0,len(d['u']),len(d['u'])); out[sp].append({k:(v[i] if isinstance(v,np.ndarray) else v) for k,v in d.items()})
    return out
A={sp:U.load('Xb55',sp) for sp in('HCN','HCOp')}; B={sp:U.load('Xb4f',sp) for sp in('HCN','HCOp')}; D={sp:U.load('Xb4d',sp) for sp in('HCN','HCOp')}
res={}
for fw in (6,10,15):
  for ul in (8,10):
    a=Rext(A,fw,ul); b=Rext(B,fw,ul); bs=[]
    for _ in range(100):
        bs.append(Rext(boot(B),fw,ul)/Rext(boot(A),fw,ul))
    # EB jackknife (drop one EB at a time)
    jk=[]
    for i in range(len(B['HCN'])):
        Bj={sp:[x for j,x in enumerate(B[sp]) if j!=i] for sp in B}; jk.append(Rext(Bj,fw,ul))
    for i in range(len(A['HCN'])):
        Aj={sp:[x for j,x in enumerate(A[sp]) if j!=i] for sp in A}; jk.append(Rext(Aj,fw,ul))
    d=Rext(D,fw,max(ul,10))
    print(f'fw={fw} u>{ul}: R_ext21={a:.2f} R_ext32(Xb4f)={b:.2f} ratio={b/a:.2f} boot sd={np.std(bs):.2f} | jackknife range of Rext over dropped EB={min(jk):.2f}-{max(jk):.2f} | Xb4d R_ext32={d:.2f} ratio={d/a:.2f}')
    res[f'{fw}|{ul}']=dict(R21=a,R32=b,ratio=b/a,boot_sd=float(np.std(bs)),jk=[min(jk),max(jk)],Xb4d=d)
json.dump(res,open('/workspace/out/q26/q26_boot.json','w'),indent=1)
