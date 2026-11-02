"""Q44 round 2: one-zone RADEX fits (s261 machinery) to the CORE (r<1.1") and RING (1.1-3.3") data vectors (s264/s265; same 2.2" cubes, +-150 km/s; no CS, no Harada).
pos tags: Nc,Nr,Sc,Sr,Cc,Cr.  12C 1-0 errors 2%, 3-2 formal (the fit adds 5% floor); 13C 1-0 amps from s265.  usage: python3 -I s266_q44_core_ring_fit.py POS[,POS] VAR [k list]; env XSTEP/XEXT as s261"""
import sys,os,json,time,numpy as np
sys.path.insert(0,'/workspace/scripts'); _a=sys.argv; sys.argv=['x']
import s261_q44_matched_fit as F
sys.argv=_a
S235=F.S235; m=F.m
O='/workspace/out/q44/'
v264=json.load(open(O+'s264.json')); v265=json.load(open(O+'s265.json')); K=F.K
RG={'N':'N','S':'S','C':'Bridge'}; AP={'c':'core','r':'ring'}
def vec(pos):
    a=v264[f'{RG[pos[0]]}|{AP[pos[1]]}']
    I=np.array([a['HCN10'],a['HCN32'][0],a['HCOp10'],a['HCOp32'][0]])*K; e=np.array([0.02*a['HCN10'],a['HCN32'][1],0.02*a['HCOp10'],a['HCOp32'][1]])*K
    return I,e,[1.,1.],[1.,1.],(1,1)
def iso13(pos):
    r=v265[f'{RG[pos[0]]}|{AP[pos[1]]}']; return dict(O=[(0,r['H13COp']['amp']*K,r['H13COp']['err']*K)],H=[(0,r['H13CN']['amp']*K,r['H13CN']['err']*K)])
F.vec=vec; F.iso13=iso13
if __name__=='__main__':
    poss=sys.argv[1].split(','); var=sys.argv[2]; ks=[int(x) for x in sys.argv[3].split(',')] if len(sys.argv)>3 else [0,1]
    for pos in poss:
      for k in ks:
        LX=[2.05,2.15,2.25] if os.environ.get('XEXT') else (S235.LOGX if not os.environ.get('XSTEP') else list(S235.LOGX[::2])+[S235.LOGX[-1]])
        for lx in LX:
            fn=O+f'fit_{pos}_{var}_s{k}_lX{int(round(lx*100))}.npz'
            if os.path.exists(fn): continue
            t=time.time(); o=F.run(pos,var,k,lx); np.savez(fn,**o); cm=o['chimin']; Tk=m.TK
            print(pos,var,S235.SIGS[k],lx,f"chi2 {o['chibest']:.2f} cold-hot {cm[Tk<=30].min()-cm[Tk>=100].min():.2f} {time.time()-t:.0f}s",flush=True)
