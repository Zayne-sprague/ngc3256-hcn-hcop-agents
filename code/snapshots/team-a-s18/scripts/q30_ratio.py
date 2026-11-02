# Q30 step 2: velocity-matched ratios, noise, nulls, tau, Tk table. Needs out/q30/spectra.pkl from q30_13co.py
import pickle,numpy as np,json,csv
r=pickle.load(open('/workspace/out/q30/spectra.pkl','rb'))
hcn={x['src']:x for x in csv.DictReader(open('/workspace/out/q5/q5e_width.csv')) if x['line']=='HCN' and x['aperture_r']=='0.0'}
def T(S,nu,bm): return 1.222e3*S*1e3/((nu/1e9)**2*bm[0]*bm[1])
out={}
for b in ['native','1.0']:
  for c in 'SN':
    d12=r[b][('12',c)]; d13=r[b][('13',c)]
    v12,v13=d12['v'],d13['v']; bm12,bm13=d12['bm'],d13['bm']
    T12=T(d12['sp'],228.4e9,bm12); T13=T(d13['sp'],218.4e9,bm13)
    # noise: line-free channels. 12CO: |v-2780|>550 (22 ch at native); 13CO: >600
    s12=np.std(T12[np.abs(v12-2780)>550]); s13=np.std(T13[np.abs(v13-2780)>600])
    # 13CO regridded to 12CO channel centres by boxcar-averaging over 12CO channel width
    dv12=d12['dv']; T13r=np.array([np.mean(T13[np.abs(v13-v)<=dv12/2]) if (np.abs(v13-v)<=dv12/2).any() else np.nan for v in v12])
    n13=np.array([ (np.abs(v13-v)<=dv12/2).sum() for v in v12]); s13r=s13/np.sqrt(np.maximum(n13,1))
    k=np.nanargmax(T12); kpk13=np.nanargmax(np.where(np.abs(v12-v12[k])<60,T13r,-9))
    v0=float(hcn[c]['v0']); fw=float(hcn[c]['fwhm'])
    m=np.abs(v12-v0)<=fw/2
    mean12=T12[m].mean(); mean13=T13r[m].mean(); e13=s13r[m].mean()/np.sqrt(m.sum())  # channels independent (not hanning-corr.)
    e12=s12/np.sqrt(m.sum())
    # integrated over full line (12CO ±(v0±400))
    mi=np.abs(v12-v0)<=400
    I12=T12[mi].sum()*dv12; I13=np.nansum(T13r[mi])*dv12; eI13=s13r[mi].mean()*np.sqrt(mi.sum())*dv12; eI12=s12*np.sqrt(mi.sum())*dv12
    # null: off positions (0.4") at same velocity window
    nl=[]
    for s12n,s13n in zip(d12['nulls'],d13['nulls']):
        t13n=T(s13n,218.4e9,bm13); t13nr=np.array([np.mean(t13n[np.abs(v13-v)<=dv12/2]) if (np.abs(v13-v)<=dv12/2).any() else np.nan for v in v12])
        nl.append([float(T(s12n,228.4e9,bm12)[m].mean()),float(np.nanmean(t13nr[m]))])
    D=dict(beam=bm12,beam13=bm13,rms12_K=float(s12),rms13_K=float(s13),
      pk12=float(T12[k]),v_pk12=float(v12[k]),T13_at_pk12=float(T13r[k]),pk13_near=float(T13r[kpk13]),v_pk13=float(v12[kpk13]),
      snr13_pk=float(T13r[kpk13]/s13r[kpk13]),
      hcn_v0=v0,hcn_fwhm=fw,nch=int(m.sum()),mean12=float(mean12),mean13=float(mean13),e12=float(e12),e13=float(e13),snr13_mean=float(mean13/e13),
      ratio_peak=float(T12[k]/T13r[k]),ratio_mean=float(mean12/mean13),ratio_mean_lo=float(mean12/(mean13+e13)),ratio_mean_hi=float(mean12/(mean13-e13)) if mean13>e13 else None,
      I12=float(I12),I13=float(I13),eI13=float(eI13),ratio_int=float(I12/I13),nulls=nl)
    out[f'{b}_{c}']=D
    print(b,c,{a:(round(x,3) if isinstance(x,float) else x) for a,x in D.items() if a!='nulls'}); print('  nulls (12,13 mean K at 0.4" offsets):',np.round(nl,3).tolist())
    np.savez(f'/workspace/out/q30/spec_{b}_{c}.npz',v=v12,T12=T12,T13=T13r,v13=v13,T13nat=T13)
json.dump(out,open('/workspace/out/q30/ratios.json','w'),indent=1)
