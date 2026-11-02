# Q13c2: HCN blue-window (2375-2575) maps at 4.5" for 3-2 (Xb4d) and 4-3 (Xb47): peak position/significance vs field rms, and spectra binned at the lobe
import sys; sys.path.insert(0,'/workspace/scripts')
from q13c_hcn43 import *
O='/workspace/out/q13/'
log=open(O+'q13c2.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
free={'HCN32':L['HCN32'][3],'HCN43':[(1790,2100),(3100,3280)]}
X0,Y0=POS['Soutflow']; XS,YS=POS['S']
for n in ('HCN32','HCN43'):
    c=load(n); r,m=base(c,free[n]); v=c['v']
    for w in ('blue','red'):
        a,b=WIN[w]; k=(v>=a)&(v<b); M=r[k].sum(0)*c['dv']
        rr=np.hypot(c['dx'],c['dy']); fld=(rr>5)&(rr<14)
        # field noise from the same window in empty field (|dx|,|dy| away from sources) : robust sigma of pixel values
        sig=1.4826*np.median(np.abs(M[fld]-np.median(M[fld])))
        near=np.hypot(c['dx']-XS,c['dy']-YS)<6; i=np.unravel_index(np.argmax(np.where(near,M,-9)),M.shape)
        j=np.unravel_index(np.argmin(np.hypot(c['dx']-X0,c['dy']-Y0)),M.shape)
        pr(f'{n} {w} window: peak within 6" of S {M[i]:.3f} Jy/b km/s ({M[i]/sig:.1f} sigma, field pixel rms {sig:.3f}) at ({c["dx"][i]:+.2f},{c["dy"][i]:+.2f})"; value at Soutflow {M[j]:.3f}; peak-Soutflow distance {np.hypot(c["dx"][i]-X0,c["dy"][i]-Y0):.2f}"')
    # binned spectrum at Soutflow (r=1")
    mm=np.hypot(c['dx']-X0,c['dy']-Y0)<=1.0; s=r[:,mm].mean(1)
    edges=np.arange(2175,2975,100.); out=[]
    for a,b in zip(edges[:-1],edges[1:]):
        k=(v>=a)&(v<b); out.append(f'{a:.0f}-{b:.0f}:{s[k].mean()*1e3:.1f}')
    pr(f'{n} spectrum at Soutflow (mJy/beam, 100 km/s bins): '+' '.join(out))
