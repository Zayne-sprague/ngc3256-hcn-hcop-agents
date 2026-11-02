# Q16r2 (7): (i) astrometry/centring of J3 (Xb4f) vs J4 (Xb47 convolved to the Xb4f beam): centroid of the N and S peaks (moment-0, 3" radius) and the effect of shifting the J4 map by up to 0.5" on T43/T32 at N, S, Bridge;
# (ii) beam-mismatch test: T43/T32 at N, S, Bridge for Xb4f-native 5.72x3.54 PA82, Xb55-native 8.64x6.00, circular 6.5", circular 10" (s30 maps), and the N/S contrast in each;
# (iii) N/S contrast significance with systematics: stat, differential recovery (+-0.05 on the ratio of ratios), centring, calibration cancels (same two bands for N and S).  usage: python3 -I s82_q16_beamtests.py -> out/q16b/s82.txt/json
import sys,json; sys.path.insert(0,'/workspace/scripts')
from q16b_lib import *
from scipy.ndimage import shift as ndshift
OUT='/workspace/out/q16b/'; log=[]
def P(s=''): print(s,flush=True); log.append(s)
res={}
def val(tag,mol,p,dshift=(0,0)):
    m=M[f'{tag}|{mol}|m0']; e,n=grid(tag,mol)
    return float(sample(tag,mol,np.array([POS[p][0]+dshift[0]]),np.array([POS[p][1]+dshift[1]]))[0])
P('# (i) centroids (E,N offset arcsec) of the moment-0 peaks within 3" of N and S, at the common Xb4f beam')
for mol in ('HCN','HCO+'):
    for p in ('N','S'):
        out={}
        for tag in ('4f3','47_4'):
            e,n=grid(tag,mol); m=np.nan_to_num(M[f'{tag}|{mol}|m0']); w=np.hypot(e-POS[p][0],n-POS[p][1])<3.0; mm=np.where(w,np.clip(m-0.5*m[w].max(),0,None),0)
            out[tag]=((mm*e).sum()/mm.sum(),(mm*n).sum()/mm.sum())
        d=(out['47_4'][0]-out['4f3'][0],out['47_4'][1]-out['4f3'][1]); P(f'  {mol:5s} {p}: J3 centroid ({out["4f3"][0]:.2f},{out["4f3"][1]:.2f})  J4 ({out["47_4"][0]:.2f},{out["47_4"][1]:.2f})  J4-J3 = ({d[0]:+.2f},{d[1]:+.2f})"')
        res[f'centroid|{mol}|{p}']=dict(J3=out['4f3'],J4=out['47_4'])
P('\n# T43/T32 (Tb) sensitivity to a shift of the J4 map relative to J3 (arcsec; E,N): reading at the fixed position')
for mol in ('HCN','HCO+'):
    for p in ('N','S','Bridge'):
        r0=ratio43(val('47_4',mol,p),val('4f3',mol,p),mol); s=[]
        for sh in ((0.5,0),(-0.5,0),(0,0.5),(0,-0.5)):
            s.append(ratio43(val('47_4',mol,p,sh),val('4f3',mol,p),mol)-r0)
        P(f'  {mol:5s} {p:6s} T43/T32 {r0:.3f}; changes for J4 shifted 0.5" E/W/N/S: '+' '.join(f'{x:+.3f}' for x in s)+f'  -> |max| {max(abs(x) for x in s):.3f}')
        res[f'shift|{mol}|{p}']=dict(r0=r0,d=s)
P('\n# (ii) beam-mismatch test')
S30={r['tag']:r for r in json.load(open('/workspace/out/s30_ratios.json'))}
S72=json.load(open('/workspace/out/q16/s72.json'))
tab={}
for mol,key in (('HCN','HCN'),('HCO+','HCOp')):
    for p in ('N','S','Bridge'):
        row={}
        k72='Xb4f native 5.72x3.54 PA82 (J3,J4)|%s|%s'%(mol,p); t=S72[k72]['tb']; row['Xb4f native 5.7x3.5']=t['J4']/t['J3']
        k72='Xb55 native 8.64x6.00 PA-90 (J2,J3,J4)|%s|%s'%(mol,p); t=S72[k72]['tb']; row['Xb55 native 8.6x6.0']=t['J4']/t['J3']
        for fw in ('6.5','10'):
            a=S30[f'Xb47_J4_{fw}as']['rows'][p+'_pos'][key]; b=S30[f'Xb4f_J3_{fw}as']['rows'][p+'_pos'][key]
            row[f'circular {fw}"']=a/b*(NU[mol][3]/NU[mol][4])**2
        tab[(mol,p)]=row
for mol in ('HCN','HCO+'):
    for p in ('N','S','Bridge'): P(f'  {mol:5s} {p:6s} '+'  '.join(f'{k} {v:.3f}' for k,v in tab[(mol,p)].items()))
    P(f'  {mol:5s} N/S contrast: '+'  '.join(f'{k} {tab[(mol,"N")][k]/tab[(mol,"S")][k]:.3f}' for k in tab[(mol,'N')])+'   N/Bridge: '+'  '.join(f'{k} {tab[(mol,"N")][k]/tab[(mol,"Bridge")][k]:.3f}' for k in tab[(mol,'N')])+'   S/Bridge: '+'  '.join(f'{k} {tab[(mol,"S")][k]/tab[(mol,"Bridge")][k]:.3f}' for k in tab[(mol,'N')]))
res['beam_table']={f'{m}|{p}':v for (m,p),v in tab.items()}
# (iii) significance of N/S (HCN, Xb4f beam) and N/Bridge, S/Bridge
st={'HCN':dict(N=0.006,S=0.011,Bridge=0.008),'HCO+':dict(N=0.005,S=0.005,Bridge=0.004)}
P('\n# (iii) contrast significance at the Xb4f native beam; stat from s72 errors; systematic terms added in quadrature on the ratio of ratios: differential recovery 0.05 (Q12: T43/T32 true/obs 1.03-1.12 at the centre; the N-S differential is taken as +-0.05), centring (max shift effect above), calibration cancels (same Xb4f and Xb47 data at both positions)')
for mol in ('HCN','HCO+'):
    for a,b in (('N','S'),('N','Bridge'),('S','Bridge')):
        ra,rb=tab[(mol,a)]['Xb4f native 5.7x3.5'],tab[(mol,b)]['Xb4f native 5.7x3.5']; c=ra/rb
        es=c*np.hypot(st[mol][a]/ra,st[mol][b]/rb)
        sh=max(max(abs(x) for x in res[f'shift|{mol}|{a}']['d']),0)/ra; sh2=max(abs(x) for x in res[f'shift|{mol}|{b}']['d'])/rb
        ec=c*np.hypot(sh,sh2)*0.4   # effect scaled linearly from the 0.5" test shift to a 0.2" rms registration error (measured J4-J3 centroid differences are <=0.21")
        es_tot=np.hypot(np.hypot(es,0.05),ec)
        P(f'  {mol:5s} {a}/{b}: {c:.3f} +- {es:.3f} (stat) +- 0.05 (recovery) +- {ec:.3f} (centring 0.2" rms) = +- {es_tot:.3f} -> ({c-1:.3f})/{es_tot:.3f} = {(c-1)/es_tot:.1f} sigma from 1')
        res[f'contrast|{mol}|{a}/{b}']=dict(c=c,stat=es,total=es_tot,sigma=(c-1)/es_tot)
json.dump(res,open(OUT+'s82.json','w'),indent=1,default=float); open(OUT+'s82.txt','w').write('\n'.join(log))
