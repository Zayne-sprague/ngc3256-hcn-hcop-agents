# Q22: continuum-scaled edge template: stack (N,S,Bridge)/continuum low-edge profiles in spw22 and spw18 (10-ch bins), slope of ramp amplitude vs continuum over 7 regions.
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
J=json.load(open('/workspace/out/q22/s132.json')); cont=J['cont']
regs=['N','S','Bridge','E1','E2','W','Sblue']
prof={}
for s in (22,18,16,20):
    z=np.load(f'/workspace/out/regspec/Xb55.spw{s}.cube.pbcor.npz')
    prof[s]={r:z[r+'_pix'][:120].reshape(12,10).mean(1)*1e3 for r in regs}
    prof[s]['hi']={r:z[r+'_pix'][::-1][:120].reshape(12,10).mean(1)*1e3 for r in regs}
e1=2.1*1.3/np.sqrt(10)   # mJy per 10-ch bin (rms 2.1 per 1.95 MHz ch, 1.3 correlation inflation)
print('stacked ramp/continuum (mean of N,S,Bridge) per 10-ch bin from low edge, +-',round(e1/np.sqrt(3)/14,3))
for s in (22,18,16,20):
    st=np.mean([prof[s][r]/cont[r] for r in ('N','S','Bridge')],0); sh=np.mean([prof[s]['hi'][r]/cont[r] for r in ('N','S','Bridge')],0)
    print(f' spw{s} low ',np.round(st,3)); print(f' spw{s} high',np.round(sh,3))
# slope vs continuum, spw22 mean of ch 0-100 and ch 0-30
for lab,sl in (('ch0-30',slice(0,3)),('ch0-100',slice(0,10))):
    y=np.array([prof[22][r][sl].mean() for r in regs]); x=np.array([cont[r] for r in regs]); e=np.full(len(regs),2.1*1.3/np.sqrt(10*(sl.stop-sl.start)))
    A=np.vstack([x,np.ones_like(x)]).T; W=1/e**2
    cov=np.linalg.inv(A.T@(A*W[:,None])); p=cov@(A.T@(W*y)); pe=np.sqrt(np.diag(cov))
    A0=x[:,None]; c0=1/np.sum(W*x**2); p0=c0*np.sum(W*x*y)
    chi=np.sum(W*(y-p0*x)**2)
    print(f'spw22 {lab}: y=a*cont+b: a={p[0]:.3f}+-{pe[0]:.3f} b={p[1]:.2f}+-{pe[1]:.2f}; through origin a={p0:.3f}+-{np.sqrt(c0):.3f} chi2={chi:.1f}/6')
    print('   data',np.round(y,2),'cont',np.round(x,1))
json.dump({s:{k:(v.tolist() if hasattr(v,'tolist') else {kk:vv.tolist() for kk,vv in v.items()}) for k,v in prof[s].items()} for s in prof},open('/workspace/out/q22/s133.json','w'))
