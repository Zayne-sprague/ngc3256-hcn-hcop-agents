import sys,os,numpy as np
os.environ['Q12_NP']='3000'
reg,model,dset=sys.argv[1:4]
sys.argv=['x',model,reg,dset,'g273','77']
ns={}
src=open('/workspace/scripts/q12_twocomp.py').read()
# stop before writing outputs
src=src.split("os.makedirs('/workspace/out/q12/runs'")[0]
exec(src,ns)
d=ns['d'];print('data',np.round(d,3))
# residuals at best: recompute prediction pieces
xb=ns['xbest'][None,:]; p=ns['unpack'](xb)
tr,Tb57,Tb17,fc,we=ns['model_obs'](p)
H,O=Tb57['hcn'],Tb57['hco+']
pr=[H[2]/H[1],H[3]/H[2],O[2]/O[1],O[3]/O[2]]
c17={m:Tb17[m][0] for m in Tb17};e17={m:Tb17[m][1] for m in Tb17}
best=None
import itertools
for g in ns['GG']:
  for h in ns['HH']:
    T={m:h*c17[m]+e17[m] for m in c17}
    q=pr+[T['hcn'][2]/T['hco+'][2],T['hcn'][0]/T['hco+'][0],g*T['hcn'][2]/T['hcn'][0],g*T['hco+'][2]/T['hco+'][0]]
    c2=sum(((q[i]-d[i,0])/d[i,1])**2 for i in range(8))+(np.log(g)/.1)**2+(np.log(h)/.3)**2
    if best is None or c2<best[0]: best=(c2,g,h,np.array(q).ravel())
print('chi2min',ns['cm'],'logZ',ns['logZ'],'best g,h',best[1:3]);print('pred',np.round(best[3],3));print('pull',np.round((best[3]-d[:,0])/d[:,1],2))
print('Tc pct',np.percentile(10**ns['unpack'](ns['X'])['lTc'],[16,50,84]))
