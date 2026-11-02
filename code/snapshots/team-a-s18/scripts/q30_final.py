# Q30 step 3: 2.5" matched 12CO(2-1)/(1-0)/(3-2), tau(12CO) from 12/13 ratio, Tk_min tables. Inputs: out/q30/spectra*.pkl
import pickle,numpy as np,json,csv
from scipy.optimize import brentq
r=pickle.load(open('/workspace/out/q30/spectra.pkl','rb')); o=pickle.load(open('/workspace/out/q30/spectra_other_2p5.pkl','rb'))
hcn={x['src']:x for x in csv.DictReader(open('/workspace/out/q5/q5e_width.csv')) if x['line']=='HCN' and x['aperture_r']=='0.0'}
def T(S,nu,bm): return 1.222e3*S*1e3/((nu/1e9)**2*bm[0]*bm[1])
def regrid(v,y,vt,dv): return np.array([np.nanmean(y[np.abs(v-x)<=dv/2]) if (np.abs(v-x)<=dv/2).any() else np.nan for x in vt])
res={}
B='2.5'; BS=(2.5,2.5)
for c in 'SN':
    d12=r[B][('12',c)]; d13=r[B][('13',c)]; v=d12['v']; dv=d12['dv']
    v0=float(hcn[c]['v0']); fw=float(hcn[c]['fwhm']); m=np.abs(v-v0)<=fw/2
    T12=T(d12['sp'],228.4e9,BS); T13=regrid(d13['v'],T(d13['sp'],218.4e9,BS),v,dv)
    def sp(k,nu0): 
        return regrid(o[k]['v'],T(o[k]['conv'][c],nu0,BS),v,dv), o[k]['v'], T(o[k]['conv'][c],nu0,BS)
    T10,v10,T10n=sp('10',112.9e9*1.0); T32,v32,T32n=sp('32',339.5e9)
    # rms from line-free channels (|v-v0|>550) on native-sampled convolved spectra
    def rms(vv,tt): return np.nanstd(tt[np.abs(vv-2780)>600])
    s={'12':rms(v,T12),'10':rms(v10,T10n),'32':rms(v32,T32n)}
    D={'v0':v0,'fwhm':fw}
    for k,t in [('12',T12),('13',T13),('10',T10),('32',T32)]:
        D['peak_'+k]=float(np.nanmax(np.where(np.abs(v-v0)<250,t,-9)));D['mean_'+k]=float(np.nanmean(t[m]))
    D['rms_K']={k:float(x) for k,x in s.items()}
    D['ratio_32_21_mean']=D['mean_32']/D['mean_12']; D['ratio_10_21_mean']=D['mean_10']/D['mean_12']; D['ratio_21_13_mean']=D['mean_12']/D['mean_13']
    D['ratio_32_21_peak']=D['peak_32']/D['peak_12']; D['ratio_10_21_peak']=D['peak_10']/D['peak_12']; D['ratio_21_13_peak']=D['peak_12']/D['peak_13']
    D['snr_mean_32']=float(D['mean_32']/(s['32']/np.sqrt(m.sum()))); D['snr_mean_10']=float(D['mean_10']/(s['10']/np.sqrt(m.sum())))
    res[c]=D; print(c,{a:(round(x,3) if isinstance(x,float) else x) for a,x in D.items()})
    np.savez(f'/workspace/out/q30/spec2p5_{c}.npz',v=v,T12=T12,T13=T13,v10=v10,T10=T10n,v32=v32,T32=T32n)
# ---- tau and Tk
T0=11.0716; Tbg=2.725*(1)  # 2-1: nu=228.4/230.538 GHz; z=0.0093 -> Tbg=2.75
def J(T): return T0/np.expm1(T0/T)
Jbg=J(2.75)
def tau_from_ratio(rat,R):
    f=lambda t:(1-np.exp(-t))/(1-np.exp(-t/R))-rat
    if rat>=R*(1-1e-6): return np.inf
    try: return brentq(f,1e-4,2000)
    except: return np.nan
def Tex_min(Tb,tau,ff):
    Jx=Tb/(ff*(1-np.exp(-tau)))+Jbg; return T0/np.log1p(T0/Jx)
tab=[]
Tb={'S':{'peak':31.16,'HCNmean':20.07},'N':{'peak':13.35,'HCNmean':11.56}}   # Q23 0.26" T_b (out/q23b/tk_ff.json)
ffs={'S':{'peak':0.88,'HCNmean':0.89},'N':{'peak':0.95,'HCNmean':0.96}}      # Q23 Gaussian-fit filling at 0.26"
rat1=json.load(open('/workspace/out/q30/ratios.json'))
for c in 'SN':
    for rk,rv in [('peak_ch',rat1[f'1.0_{c}']['ratio_peak']),('HCNFWHM_mean',rat1[f'1.0_{c}']['ratio_mean'])]:
        for fac,fl in [(0.8,'-20%'),(1.0,'central'),(1.2,'+20%')]:
            for R in [25,40,60,90]:
                tau=tau_from_ratio(rv*fac,R)
                for tk in ['peak','HCNmean']:
                    for ffn,ff in [('1',1.0),('fit',ffs[c][tk]),('0.5',0.5),('0.3',0.3)]:
                        tab.append(dict(core=c,ratio_kind=rk,ratio=rv*fac,ratio_variant=fl,R=R,tau12=tau,Tb_kind=tk,Tb=Tb[c][tk],ff=ffn,Tk_min=float(Tex_min(Tb[c][tk],tau,ff)) if np.isfinite(tau) else float(Tex_min(Tb[c][tk],50,ff))))
with open('/workspace/out/q30/tau_tk_table.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(tab[0])); w.writeheader(); w.writerows(tab)
json.dump(res,open('/workspace/out/q30/other_co_2p5.json','w'),indent=1)
# compact print
for c in 'SN':
  for rk in ['peak_ch','HCNFWHM_mean']:
    for R in [25,40,60,90]:
      row=[t for t in tab if t['core']==c and t['ratio_kind']==rk and t['R']==R and t['ratio_variant']=='central' and t['Tb_kind']=='peak']
      lo=[t for t in tab if t['core']==c and t['ratio_kind']==rk and t['R']==R and t['ratio_variant']=='-20%'][0]['tau12']
      hi=[t for t in tab if t['core']==c and t['ratio_kind']==rk and t['R']==R and t['ratio_variant']=='+20%'][0]['tau12']
      print(c,rk,'R',R,'ratio %.1f tau %.2f [%.2f,%.2f]'%(row[0]['ratio'],row[0]['tau12'],lo,hi),'Tk_min(peakTb) ff1/fit/.5/.3:',[round(t['Tk_min'],1) for t in row])
# ---- Tex implied by 2.5" T_b ratios (optically thick, thermal-limit J(T) ratios; same ff)
def Jn(T,t0): return t0/np.expm1(t0/T)
print('thick-LTE T_b ratios vs T:')
for T_ in [10,15,20,30,50]: print(T_,'J32/J21=%.3f J10/J21=%.3f'%(Jn(T_,16.6)-Jn(2.75,16.6) if False else (Jn(T_,16.60)-Jn(2.75,16.60))/(Jn(T_,11.07)-Jn(2.75,11.07)),(Jn(T_,5.53)-Jn(2.73,5.53))/(Jn(T_,11.07)-Jn(2.75,11.07))))
imp={}
for c in 'SN':
    D=res[c]
    for lab,t0,key in [('32',16.60,'ratio_32_21_mean'),('10',5.53,'ratio_10_21_mean')]:
        f=lambda T_:(Jn(T_,t0)-Jn(2.75,t0))/(Jn(T_,11.07)-Jn(2.75,11.07))-D[key]
        try: imp[f'{c}_{lab}']=brentq(f,3,500)
        except Exception as e: imp[f'{c}_{lab}']=None
print('implied Tex(thick,LTE,same ff) from 2.5" mean ratios:',imp)
res['implied_Tex']=imp; json.dump(res,open('/workspace/out/q30/other_co_2p5.json','w'),indent=1)
print('HCN-mean-Tb Tk_min, central ratio, R=40/60/90, ff 1/fit:')
for c in 'SN':
  for R in [25,40,60,90]:
    row=[t for t in tab if t['core']==c and t['ratio_kind']=='HCNFWHM_mean' and t['R']==R and t['ratio_variant']=='central' and t['Tb_kind']=='HCNmean']
    print(c,R,[round(t['Tk_min'],1) for t in row])
