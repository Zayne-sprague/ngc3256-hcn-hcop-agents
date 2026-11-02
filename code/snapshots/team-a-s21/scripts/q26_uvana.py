# Q26: uv-plane HCN/HCO+ comparison at matched |uv| (kilo-lambda) from out/q26/uv/*.npy (made by q26_uvdump.py).
# (1) binned Re/Im/|V| per EB and inverse-variance combined per set; ratio HCN/HCO+ per bin, 2-1 vs 3-2; null (+1000 km/s window) per bin.
# (2) linear complex LSQ: N + S point sources (positions from Q22/Xb4b) + one extended circular Gaussian (FWHM fixed 6/10/15") at the N-S midpoint.
import glob, json, numpy as np
D='/workspace/out/q26/uv/'; ARC=np.pi/180/3600
BINS=[6,8,10,14,20,30]; PK={'N':(0.,2.55),'S':(-0.15,-2.55)}; MID=(-0.075,0.)
def load(o,sp):
    L=[]
    for f in sorted(glob.glob(f'{D}{o}_*_{sp}_*.npy')):
        a=np.load(f); u,v,re,im,sg,nr,ni=a[:7]; k=np.hypot(u,v)>0
        k&=np.isfinite(re)&np.isfinite(im)&np.isfinite(sg)&(sg>0)&(np.hypot(u,v)>1e3)
        L.append(dict(eb=f.split('_uid___')[1][:12],u=u[k],v=v[k],re=re[k],im=im[k],sg=sg[k],nr=nr[k],ni=ni[k]))
    return L
def binstat(d,lo,hi,key='re'):
    k=np.hypot(d['u'],d['v'])/1e3; m=(k>=lo)&(k<hi)
    if m.sum()<20: return None
    # phase-rotate to the N-S midpoint (shift of phase centre)
    ph=np.exp(2j*np.pi*(d['u']*MID[0]+d['v']*MID[1])*ARC)
    z=(d['re']+1j*d['im'])*ph; w=1/d['sg'][m]**2
    zn=(d['nr'][m]+1j*d['ni'][m])
    # error: empirical scatter of the mean (robust: includes real structure) and formal
    mre=np.sum(w*z.real[m])/w.sum(); mim=np.sum(w*z.imag[m])/w.sum(); fe=1/np.sqrt(w.sum())
    se=np.sqrt(np.sum(w**2*(z.real[m]-mre)**2))/w.sum()
    amp=np.sum(w*np.abs(z[m]))/w.sum()
    nn=np.isfinite(zn.real); mn=np.sum(w[nn]*zn.real[nn])/w[nn].sum() if nn.sum()>20 else np.nan
    return dict(n=int(m.sum()),re=mre,im=mim,err=max(fe,se),fe=fe,amp=amp,null=mn)
def comb(L,lo,hi):
    s=[binstat(d,lo,hi) for d in L]; s=[x for x in s if x]
    if not s: return None,[]
    w=np.array([1/x['err']**2 for x in s]); g=lambda k:np.sum([wi*x[k] for wi,x in zip(w,s)])/w.sum()
    ebs=[x['re'] for x in s]; sp=np.std(ebs,ddof=1) if len(ebs)>1 else np.nan
    return dict(re=g('re'),im=g('im'),err=1/np.sqrt(w.sum()),eb_spread=sp,amp=g('amp'),null=np.nanmean([x['null'] for x in s]),n=sum(x['n'] for x in s)),s
def design(d,fw):
    u,v=d['u']*ARC,d['v']*ARC; cols=[]
    for p in PK.values(): cols.append(np.exp(-2j*np.pi*(u*p[0]+u*0+v*p[1])))
    ph=np.exp(-2j*np.pi*(u*MID[0]+v*MID[1])); s=fw*ARC/(2*np.sqrt(2*np.log(2)))
    cols.append(ph*np.exp(-2*np.pi**2*s**2*(d['u']**2+d['v']**2)))
    return np.array(cols).T
def fit(L,fw,ulo=6,uhi=30):
    A=[];b=[];w=[]
    for d in L:
        k=np.hypot(d['u'],d['v'])/1e3; m=(k>=ulo)&(k<uhi)
        dd={x:d[x][m] for x in('u','v')}; X=design(dd,fw)
        A.append(X); b.append(d['re'][m]+1j*d['im'][m]); w.append(1/d['sg'][m])
    A=np.vstack(A);b=np.concatenate(b);w=np.concatenate(w)
    # real-valued flux parameters: stack Re/Im
    Ar=np.vstack([A.real*w[:,None],A.imag*w[:,None]]); br=np.concatenate([b.real*w,b.imag*w])
    p,res,_,_=np.linalg.lstsq(Ar,br,rcond=None); chi2=np.sum((Ar@p-br)**2); dof=len(br)-3
    cov=np.linalg.inv(Ar.T@Ar)*max(chi2/dof,1)
    return p,np.sqrt(np.diag(cov)),chi2/dof
if __name__=='__main__':
    SET=[('Xb55','2-1'),('Xb4f','3-2'),('Xb4d','3-2')]; out={'bins':{},'fit':{}}
    for o,J in SET:
        for lo,hi in zip(BINS[:-1],BINS[1:]):
            r={}
            for sp in('HCN','HCOp'):
                c,s=comb(load(o,sp),lo,hi); r[sp]=c; r[sp+'_eb']=[x['re'] for x in s]
            if r['HCN'] and r['HCOp']:
                h,q=r['HCN'],r['HCOp']; R=h['re']/q['re']; eR=abs(R)*np.hypot(h['err']/h['re'],q['err']/q['re'])
                print(f'{o} {J} {lo}-{hi}kl n={h["n"]:6d} HCN Re={h["re"]:7.3f}+-{h["err"]:.3f} (Im {h["im"]:.3f}, null {h["null"]:.3f}, EBspread {h["eb_spread"]:.3f}) HCO+ Re={q["re"]:7.3f}+-{q["err"]:.3f} (Im {q["im"]:.3f}, null {q["null"]:.3f}) R={R:.3f}+-{eR:.3f}')
                out['bins'][f'{o}|{J}|{lo}-{hi}']=dict(HCN=h,HCO=q,R=R,eR=eR,HCN_eb=r['HCN_eb'],HCO_eb=r['HCOp_eb'])
    for o,J in SET:
        for fw in (6,10,15):
            for ulo in (6,8,10):
                if o=='Xb4d' and ulo<10: continue
                res={}
                for sp in('HCN','HCOp'):
                    L=load(o,sp)
                    p,e,c=fit(L,fw,ulo); res[sp]=(p,e,c)
                (ph,eh,ch),(pq,eq,cq)=res['HCN'],res['HCOp']
                print(f'FIT {o} {J} extFWHM={fw}" u>{ulo}: HCN N,S,ext={ph.round(2)}+-{eh.round(2)} chi2r {ch:.1f} | HCO+ {pq.round(2)}+-{eq.round(2)} chi2r {cq:.1f} | R(ext)={ph[2]/pq[2]:.2f} R(N)={ph[0]/pq[0]:.2f} R(S)={ph[1]/pq[1]:.2f}')
                out['fit'][f'{o}|{J}|{fw}|{ulo}']=dict(HCN=[ph.tolist(),eh.tolist(),ch],HCO=[pq.tolist(),eq.tolist(),cq])
    json.dump(out,open('/workspace/out/q26/q26_uvana.json','w'),indent=1,default=float)
