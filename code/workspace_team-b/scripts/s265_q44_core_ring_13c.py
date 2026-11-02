"""Q44 round 2: 13C/12C core (r<1.1") vs ring (1.1-3.3") at N, S, Bridge directly (independent, non-overlapping apertures; null-calibrated errors as s260), plus
 a noise-injection null test: the null apertures (s185 'nulls', no line emission) give the empirical noise of core-minus-ring ratio difference (null amplitudes divided by the target's 12C refs);
 false-gradient rate = fraction of null positions with Delta_null >= Delta_obs (one-sided), per line and for the sum of the two lines.  usage: python3 -I s265_q44_core_ring_13c.py -> out/q44/s265.json,s265.txt"""
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s185_q29_weaklines as S
O='/workspace/out/q44/'; L=[]
def P(s=''): print(s,flush=True); L.append(s)
BOX=(-150,150)
def ap_ann(Tb,pix,cen,r0,r1):
    ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]; d=np.hypot(X-cen[0],Y-cen[1])*pix; m=(d<=r1)&(d>=r0) if r0>0 else (d<=r1); return Tb[:,m].mean(axis=1)
def spectra_r(pos,r0,r1):
    out={}
    for k,(Tb,nu,pix,bm,pxy) in S.DATA.items():
        x,y=pxy(pos); out[k]=ap_ann(Tb,pix,(x,y),r0,r1)
    return out
SEL={'H13COp':'comb','H13CN':'X95d_spw0'}
APS={'core':(0,1.1),'ring':(1.1,3.3)}
nul={a:{} for a in APS}; 
for an,(r0,r1) in APS.items():
    for (rr,ang,p,mi) in S.nulls:
        for kk,(a,e) in S.boxfit(spectra_r(p,r0,r1),BOX).items():
            if kk[0] in SEL and kk[1]==SEL[kk[0]]: nul[an].setdefault(kk[0],[]).append((a,e))
P(f'{len(S.nulls)} null apertures')
res={}
for rn in ('N','S','Bridge'):
    pos=S.positions(rn)
    for an,(r0,r1) in APS.items():
        sp=spectra_r(pos,r0,r1); R=S.boxfit(sp,BOX); IO=S.refint(sp['X961_spw0'],S.F0['HCOp'],BOX); IH=S.refint(sp['X961_spw0'],S.F0['HCN'],BOX)
        row={'HCOp_int':IO,'HCN_int':IH}
        for kk,(a,e) in R.items():
            l=kk[0]
            if l not in SEL or kk[1]!=SEL[l]: continue
            Ln=np.array(nul[an][l]); c=max(1.0,float(np.sqrt(np.mean((Ln[:,0]/Ln[:,1])**2)))); e2=max(e*c,float(np.std(Ln[:,0])))
            ref=IH if l=='H13CN' else IO; row[l]=dict(amp=a,err=e2,ratio=a/ref,ratio_err=e2/ref)
        res[f'{rn}|{an}']=row
        P(f'{rn:6s} {an:5s} HCO+ {IO:6.1f} HCN {IH:6.1f} | H13CO+ {row["H13COp"]["amp"]:.2f}+-{row["H13COp"]["err"]:.2f} ratio {row["H13COp"]["ratio"]:.4f}+-{row["H13COp"]["ratio_err"]:.4f} | H13CN {row["H13CN"]["amp"]:.2f}+-{row["H13CN"]["err"]:.2f} ratio {row["H13CN"]["ratio"]:.4f}+-{row["H13CN"]["ratio_err"]:.4f}')
P('\ncore-minus-ring ratio difference (independent apertures; 8% scale cancels), formal significance and noise-injection false-gradient rate (null positions)')
rng=np.random.default_rng(1)
for rn in ('N','S','Bridge'):
    tot=0; ev=0; nulD=0; D={}
    for l,ref in (('H13COp','HCOp_int'),('H13CN','HCN_int')):
        c=res[f'{rn}|core']; r=res[f'{rn}|ring']; d=c[l]['ratio']-r[l]['ratio']; e=np.hypot(c[l]['ratio_err'],r[l]['ratio_err'])
        nc=np.array(nul['core'][l])[:,0]; nr=np.array(nul['ring'][l])[:,0]
        dn=nc/c[ref]-nr/r[ref]   # null differences with this target's reference intensities (zero true gradient, zero signal)
        pfa=float((dn>=d).mean()); D[l]=(d,e,pfa,float(np.std(dn)))
        P(f'  {rn:6s} {l:7s}: core-ring = {d:+.4f} +- {e:.4f} ({d/e:.1f} sigma); empirical null sigma {np.std(dn):.4f}; P(null >= obs) = {pfa:.3f} ({len(dn)} nulls)')
        tot+=d/e**2; ev+=1/e**2; nulD=nulD+dn/e**2 if not isinstance(nulD,int) else dn/e**2
    dw=tot/ev; ew=ev**-0.5; dnw=nulD/ev
    P(f'  {rn:6s} weighted two-line mean: {dw:+.4f} +- {ew:.4f} ({dw/ew:.1f} sigma); null P = {(dnw>=dw).mean():.3f} (lines share apertures; the two H13C amplitudes are independent measurements of different species)')
    res[f'diff|{rn}']={k:dict(d=v[0],e=v[1],pnull=v[2],snull=v[3]) for k,v in D.items()}; res[f'diff|{rn}']['wmean']=dict(d=dw,e=ew,pnull=float((dnw>=dw).mean()))
json.dump(res,open(O+'s265.json','w'),indent=1,default=float); open(O+'s265.txt','w').write('\n'.join(L))
