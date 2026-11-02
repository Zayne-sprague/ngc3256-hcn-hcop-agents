# Q11 forward model B fit: constrain the recovery scale f and the extra intrinsic HCN fall eps (relative to the 2-1 sky proxy) from the 3-2 annulus ratios (and, mode A, HCO+ flux profile).
# eps = fractional HCN/HCO+ fall, ramping from 0 at r=4" to eps at r>=10", ADDITIONAL to the 2-1 gradient. Sigma: R 0.01,0.01,0.06 (8-12": zero level/beam), HCO+ F(8-12)/F(4-8) 0.06 (mode A).
import sys,json; sys.path.insert(0,'/workspace/scripts')
src=open('/workspace/scripts/q11_proxy.py').read().split("res=[]")[0]; exec(src)
sig=np.array([0.01,0.01,0.06]); out={}
for lab,pbfw,Bm,tag in (('Xb4f PB37 (7m-like)',37.,8.7,'Xb4f'),('Xb4f PB22 (12m-like)',22.,8.7,'Xb4f'),('Xb4d PB22',22.,8.7,'Xb4d')):
    fs=np.arange(0.05,0.41,0.025); es=np.arange(-0.2,0.81,0.05); C=np.zeros((len(fs),len(es))); CA=np.zeros_like(C); Rs={}
    # observed HCO+ F(8-12)/F(4-8)
    kk='HCO+ 3-2@'+tag; o,so=mp(kk); r=np.hypot(o['dx']-XM,o['dy']); om=OM(tb)/o['pix']**2
    fo=np.array([so[(r>=x)&(r<y)].sum()/om for x,y in ANN]); fobs=fo[2]/fo[1]
    for i,f in enumerate(fs):
        ua=f*Bm/(c/265.9e9)/1e3; o_=filt(sO,pbfw,ua); b=ann(o_); fm=b[2]/b[1]
        for j,e in enumerate(es):
            a=ann(filt(sH*(1-e*np.clip((Rc-4)/6,0,1)),pbfw,ua)); Rm=a/b; Rm=Rm*obs[tag][0]/Rm[0]
            C[i,j]=(((Rm-obs[tag])/sig)**2).sum(); CA[i,j]=C[i,j]+((fm-fobs)/0.06)**2
    for mode,CC in (('B f free [0.05-0.4]',C),('A f from HCO+ flux profile',CA)):
        i,j=np.unravel_index(CC.argmin(),CC.shape); prof=CC.min(axis=0); ok=es[prof<=prof.min()+1.0]
        ok2=es[prof<=prof.min()+4.0]
        print(f'{lab:22s} {mode:28s} chi2min {CC.min():5.1f} at f={fs[i]:.3f} eps={es[j]:+.2f}; eps 1sig [{ok.min():+.2f},{ok.max():+.2f}] 2sig [{ok2.min():+.2f},{ok2.max():+.2f}]  (obs HCO+ F8-12/F4-8 {fobs:.3f})')
        out[lab+'|'+mode]=dict(chi2=float(CC.min()),f=float(fs[i]),eps=float(es[j]),eps1=[float(ok.min()),float(ok.max())],eps2=[float(ok2.min()),float(ok2.max())])
json.dump(out,open('/workspace/out/q11/q11_proxy_fit.json','w'),indent=1)
