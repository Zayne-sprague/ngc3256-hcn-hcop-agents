# Q13a: reconcile Xb4d vs Xb4b for the S-outflow blue-lobe HCN excess.
# Both cubes are convolved channel-wise to the same circular 1.05" beam, baseline-subtracted per pixel (linear, line-free channels in each line's frame),
# window maps made in the same velocity windows (radio, own rest frequency), E = HCN - k*HCO+ (k from the lobe's own core window or fixed from Xb4d).
# Apertures: r=0.5", 1.0" (as team-b) and a box. Noise: random apertures of the same size over the field (bright core emission masked).
# Outputs: out/q13/q13a_maps_<ous>.npz, q13a.log, q13a_table.csv
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
import pandas as pd
O='/workspace/out/q13/'
WIN={'blue':(2375,2575),'blueB':(2300,2575),'core':(2575,2875),'red':(2975,3175)}   # HCN/HCO+ own frames, radio km/s; blue = team-b/Q3b excess window
def pixbase(c,free):
    v=c['v']; m=np.zeros(len(v),bool)
    for a,b in free: m|=(v>=a)&(v<=b)
    X=np.vander(v-v[m].mean(),2); nz=c['a'].shape[1:]
    A=c['a'][m].reshape(m.sum(),-1); coef=np.linalg.lstsq(X[m],A,rcond=None)[0]
    r=c['a']-(X@coef).reshape(len(v),*nz)
    return r
def wmaps(ous,fw=1.05,crop=9.):
    out={}
    for sp in ('HCN','HCO+'):
        c=vcube(sp,ous,fw,crop); r=pixbase(c,FREE[sp]); v=c['v']; dv=c['dv']
        for k,(a,b) in WIN.items():
            m=(v>=a)&(v<b); out[f'{sp}_{k}']=r[m].sum(0)*dv; out[f'{sp}_{k}_n']=m.sum()
        # per-pixel channel-rms in line-free channels -> formal noise of window sum
        fm=np.zeros(len(v),bool)
        for a,b in FREE[sp]: fm|=(v>=a)&(v<=b)
        out[f'{sp}_rms']=1.4826*np.median(np.abs(r[fm]-np.median(r[fm],0)),0)
        out['dx']=c['dx']; out['dy']=c['dy']; out['dv']=dv; out['pix']=c['pix']
        out[f'{sp}_spec']=r; out['v_'+sp]=v
    return out
def apmean(M,dx,dy,x,y,rad):
    m=np.hypot(dx-x,dy-y)<=rad; return M[m].mean()
if __name__=='__main__':
    log=open(O+'q13a.log','w')
    def pr(*a):
        s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
    maps={}
    for ous in ('Xb4d','Xb4b'):
        maps[ous]=wmaps(ous); np.savez_compressed(O+f'q13a_maps_{ous}.npz',**{k:v for k,v in maps[ous].items() if not k.endswith('_spec')})
    dx=maps['Xb4d']['dx']; dy=maps['Xb4d']['dy']
    pr('# both arrays convolved to 1.05" circular; windows (radio km/s, own frame):',WIN)
    rows=[]
    # k per array from lobe core window in aperture r=1": k=HCN/HCO+ (core), plus fixed k=0.574 (Xb4d fit of Q3b)
    X0,Y0=POS['Soutflow']
    # noise: random apertures; exclude positions within 3.5" of N, S, Bridge, Soutflow, W clump and any aperture with core HCO+ > 3x median-of-field
    rng=np.random.default_rng(13)
    for ous in ('Xb4d','Xb4b'):
        Mp=maps[ous]; dxo=Mp['dx']; dyo=Mp['dy']
        ex=[POS['N'],POS['S'],POS['Bridge'],POS['Soutflow'],(-1.8,-2.4)]
        for rad in (0.5,1.0):
            # candidate centres
            cands=[]
            for _ in range(4000):
                x,y=rng.uniform(-7.5,7.5,2)
                if min(np.hypot(x-a,y-b) for a,b in ex)<3.0: continue
                if np.hypot(x,y)>7.5: continue
                cands.append((x,y))
            Ob=Mp['HCO+_core'];
            # mask bright: HCO+ core ap mean > 0.25 Jy/beam km/s (~5 sigma-ish) excluded to avoid real wings
            good=[(x,y) for x,y in cands if apmean(Ob,dxo,dyo,x,y,rad)<0.25]
            good=good[:300]
            Mp[f'nullpos_{rad}']=np.array(good)
        pr(f'{ous}: null apertures r=0.5": {len(Mp["nullpos_0.5"])}, r=1.0": {len(Mp["nullpos_1.0"])}')
    for rad in (0.5,1.0):
        for ous in ('Xb4d','Xb4b'):
            Mp=maps[ous]; P=Mp[f'nullpos_{rad}']
            def ap(key,x,y): return apmean(Mp[key],Mp['dx'],Mp['dy'],x,y,rad)
            kcore=ap('HCN_core',X0,Y0)/ap('HCO+_core',X0,Y0)
            for kname,k in (('k_own_core',kcore),('k=0.574(Xb4d fit)',0.574)):
                # null distribution of E_blue and of each line's window sums at random apertures
                nul={}
                for key in ('HCN_blue','HCO+_blue','HCN_blueB','HCO+_blueB','HCN_red','HCO+_red'):
                    nul[key]=np.array([ap(key,x,y) for x,y in P])
                Enul=nul['HCN_blue']-k*nul['HCO+_blue']; EnulB=nul['HCN_blueB']-k*nul['HCO+_blueB']; Ered=nul['HCN_red']-k*nul['HCO+_red']
                sE=1.4826*np.median(np.abs(Enul-np.median(Enul))); sEB=1.4826*np.median(np.abs(EnulB-np.median(EnulB)))
                Eb=ap('HCN_blue',X0,Y0)-k*ap('HCO+_blue',X0,Y0); EbB=ap('HCN_blueB',X0,Y0)-k*ap('HCO+_blueB',X0,Y0)
                Er=ap('HCN_red',X0,Y0)-k*ap('HCO+_red',X0,Y0)
                rows.append(dict(ous=ous,rad=rad,k=kname,kval=k,HCN_blue=ap('HCN_blue',X0,Y0),HCO_blue=ap('HCO+_blue',X0,Y0),
                    sHCN=1.4826*np.median(np.abs(nul['HCN_blue']-np.median(nul['HCN_blue']))),sHCO=1.4826*np.median(np.abs(nul['HCO+_blue']-np.median(nul['HCO+_blue']))),
                    HCN_core=ap('HCN_core',X0,Y0),HCO_core=ap('HCO+_core',X0,Y0),
                    E_blue=Eb,sig_E=sE,snr_E=Eb/sE,E_blueB=EbB,snr_EB=EbB/sEB,E_red=Er,snr_Ered=Er/sE,
                    null_med_E=np.median(Enul),null_n=len(P),n_null_ge=int((Enul>=Eb).sum()),n_null_le_neg=int((Enul<=-Eb).sum())))
    T=pd.DataFrame(rows); T.round(4).to_csv(O+'q13a_table.csv',index=False)
    pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
    pr(T.round(3).to_string(index=False))
