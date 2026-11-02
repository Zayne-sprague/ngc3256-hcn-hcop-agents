import pickle,numpy as np
r=pickle.load(open('/workspace/out/q30/spectra.pkl','rb'))
for b in ['native','1.0']:
  for c in 'SN':
    for l in ['12','13']:
        d=r[b][(l,c)]; v=d['v']; sp=d['sp']; bm=d['bm']
        off=(np.abs(v-2780)>600)
        s=np.nanstd(sp[off]); k=np.nanargmax(sp)
        T=lambda S: 1.222e3*S*1e3/((d['nu'][k]/1e9)**2*bm[0]*bm[1])
        print(b,c,l,'peak %.1f mJy at v=%.0f rms %.2f mJy Tb=%.2f K rms=%.3f K nfree=%d'%(sp[k]*1e3,v[k],s*1e3,T(sp[k]),T(s),off.sum()), 'vrange',v.min(),v.max())
