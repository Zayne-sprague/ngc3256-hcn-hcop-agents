# CASA script: dump continuum-subtracted, line-integrated (VSYS+-HW km/s) visibilities of HCN/HCO+ from the restored split MSs to npy.
# columns: u,v (wavelengths at the line freq), Re, Im (Jy km/s), sigma (Jy km/s, from off-line scatter), Re_null, Im_null (same window width shifted to
# +NULLOFF km/s, line-free control), time(s). Also prints/saves per-DD diagnostics and the baseline-length coverage of each EB.
import glob, json, os, numpy as np
from casatools import table
C=299792.458; VS=2780.; HW=350.; NULLOFF=1000.
LINES={'HCN':{'2-1':177.2611,'3-2':265.8864},'HCO+':{'2-1':178.3751,'3-2':267.5576}}
OUT='/workspace/out/q26/uv/'; os.makedirs(OUT,exist_ok=True)
groups={}
for s in glob.glob('/archive/derived/*/*/summary.json'):
    j=json.load(open(s)); d=os.path.dirname(s); o=j['ous'].split('/')[-1]
    if o in ('Xb55','Xb4f','Xb4d'): groups.setdefault(o,[]).extend(sorted(glob.glob(d+'/*.split.ms')))
J={'Xb55':'2-1','Xb4f':'3-2','Xb4d':'3-2'}
cov={}
for o,mss in groups.items():
    for ms in sorted(set(mss)):
        eb=os.path.basename(ms).split('.')[0]
        tb=table(); tb.open(ms+'/SPECTRAL_WINDOW'); fr=[tb.getcell('CHAN_FREQ',i) for i in range(tb.nrows())]; tb.close()
        tb.open(ms+'/DATA_DESCRIPTION'); spw_of=tb.getcol('SPECTRAL_WINDOW_ID'); tb.close()
        tb.open(ms)
        print(o,eb,'nrows',tb.nrows(),'DD->spw',list(spw_of),flush=True)
        for sp,fr0 in LINES.items():
            f0=fr0[J[o]]*1e9; res=[]
            for dd_i,spw in enumerate(spw_of):
                v=(1-fr[spw]/f0)*C
                if v.min()>VS-HW-50 or v.max()<VS+HW+50:
                    print('  skip',sp,'dd',dd_i,'v range',v.min(),v.max()); continue
                ln=np.where(abs(v-VS)<=HW)[0]; nl=np.where(abs(v-VS-NULLOFF)<=HW)[0]
                off=np.where((abs(v-VS)>600)&(abs(v-VS-NULLOFF)>HW+100))[0]
                print('  use',sp,'dd',dd_i,'spw',spw,'nline',len(ln),'nnull',len(nl),'noff',len(off),flush=True)
                if len(ln)<5 or len(off)<10: continue
                sub=tb.query(f'DATA_DESC_ID=={dd_i}')
                if sub.nrows()==0: sub.close(); print('  no rows'); continue
                data=sub.getcol('DATA'); flag=sub.getcol('FLAG'); uv=sub.getcol('UVW'); t=sub.getcol('TIME'); sub.close()
                dv=abs(np.median(np.diff(v)))
                dm=np.where(flag,np.nan,data)                       # (pol,chan,row)
                cont=np.nanmean(dm[:,off,:],axis=1)                 # (pol,row)
                def integ(idx):
                    x=dm[:,idx,:]-cont[:,None,:]
                    nok=np.isfinite(x).sum(1)
                    s=np.nansum(x,axis=1)*dv*len(idx)/np.maximum(nok,1)
                    return s,nok>=0.8*len(idx)
                lin,ok=integ(ln); nul,okn=integ(nl) if len(nl)>=5 else (lin*np.nan,ok)
                sig=np.nanstd(dm[:,off,:].real,axis=1)*np.sqrt(len(ln))*dv   # per pol, per row
                ok=ok&np.isfinite(cont)
                V=np.where(ok,lin,np.nan); N=np.where(ok&okn,nul,np.nan)
                npol=np.maximum(ok.sum(0),1)
                Vm=np.nanmean(V,axis=0); Nm=np.nanmean(N,axis=0)
                sm=np.sqrt(np.nansum(np.where(ok,sig,0)**2,axis=0))/npol
                lam=C*1e3/(f0*(1-VS/C))
                res.append(np.array([uv[0]/lam,uv[1]/lam,Vm.real,Vm.imag,sm,Nm.real,Nm.imag,t]))
                b=np.hypot(uv[0],uv[1])/lam/1e3
                cov[f'{o}|{eb}|{sp}']=np.histogram(b,bins=[0,2,4,6,8,10,14,20,30,60,200])[0].tolist()
            if res:
                a=np.hstack(res); np.save(f'{OUT}{o}_{eb}_{sp.replace("+","p")}_{J[o]}.npy',a); print(o,eb,sp,a.shape,flush=True)
        tb.close()
json.dump({'bins_kl':[0,2,4,6,8,10,14,20,30,60,200],'counts':cov},open(OUT+'coverage.json','w'),indent=1)
