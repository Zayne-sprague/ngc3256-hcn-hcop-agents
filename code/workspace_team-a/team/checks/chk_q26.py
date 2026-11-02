import numpy as np,glob
def load(p):
    a=np.load(p); u,v,re,im,s=a[:5]; return np.hypot(u,v)/1e3,re,s
def binr(tag,lines,edges):
    out=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        r=[]
        for l in lines:
            R=[];W=[]
            for p in glob.glob(f'out/q26/uv/{tag}*_{l}_*.npy'):
                k,re,s=load(p); m=(k>=lo)&(k<hi)&(s>0)&np.isfinite(re)
                R+=list(re[m]);W+=list(1/s[m]**2)
            R=np.array(R);W=np.array(W)
            r.append(np.sum(R*W)/np.sum(W) if len(R) else np.nan)
        out.append((lo,hi,r[0],r[1],r[0]/r[1]))
    return out
e=[6,8,10,14,20]
for t in ['Xb55','Xb4f','Xb4d']: print(t,[tuple(round(x,2) for x in o) for o in binr(t,['HCN','HCOp'],e)])
