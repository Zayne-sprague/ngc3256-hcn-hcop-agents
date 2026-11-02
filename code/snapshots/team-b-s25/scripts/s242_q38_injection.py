"""Q38 fix 3: injection test of the s58d NNLS de-contamination (clump totals + nuclear totals) on SIMULATED skies passed through assumed ACA-like uv coverage
(s58f machinery: gridded uv, Hogbom clean, PB, restored native-beam map) for Xb55 (J=2), Xb4f (J=3), Xb47 (J=4).  Truth: HCO+ sky = s58d V0 Xb4f fit components (+ optional unmodelled
extended Gaussian), HCN sky = R_comp(J) x HCO+ sky.  Truth ladders: 'flat' (clumps 0.50, nucleus 0.80 at every J) and 'decline' (clump 0.50/0.35/0.20, nucleus 0.80/0.65/0.50).
Output: recovered clump/nuclear R(J), R32/R21, R43/R21 vs truth, and the 10" position-value R at E1 for comparison.  usage: python3 -I s242_q38_injection.py -> out/q38/s242.txt/json"""
import sys,json; sys.path.insert(0,'/workspace/scripts')
src=open('/workspace/scripts/s58f_simulate.py').read().split("OUT={}")[0]
exec(src)      # defines sky, observe, read10, true_sky, LAMS, POS, ...
import s58d_clumpfit as D
log2=[]
def P2(s=''): print(s,flush=True); log2.append(s)
RC={'flat':{2:0.5,3:0.5,4:0.5},'decline':{2:0.5,3:0.35,4:0.20}}; RN={'flat':{2:0.8,3:0.8,4:0.8},'decline':{2:0.8,3:0.65,4:0.50}}
ARR={2:('Xb55',8.7,45,'B5',56.4),3:('Xb4f',8.7,45,'B6',37.6),4:('Xb47',8.7,45,'B7',28.2)}
def true_sky_r(fext,wext,rc,rn):
    def g(p,fw,fl,r): return gauss(*sky_pos(p),fw,fl*r)
    I=g('N',1.0,CF['N'],rn)+g('S',1.0,CF['S'],rn)+g('Bridge',2.0,CF['Bridge'],rn)+gauss(0,0,5,CF.get('haloA',0)*rn)+gauss(0,0,14,CF['haloB']*rn)
    for c in ('E1','E2','W'): I=I+g(c,2.5,CF[c],rc)
    if fext>0: I=I+gauss(0,0,wext,fext*CF['nuc_total']*rn)
    return I
def sky_hco(fext,wext): return true_sky_r(fext,wext,1.0,1.0)
sub=2; c0=NG//2; n19=int(19/PX); sl=slice(c0-n19,c0+n19+1,sub)
dxg=(np.arange(NG)-c0)*PX
DX,DY=np.meshgrid(dxg[sl],dxg[sl])
OUT={}
for sn,fe,we in (('no extended',0,40),('ext 25" f=0.3',0.3,25),('ext 40" f=0.5',0.5,40)):
  for tn in ('flat','decline'):
    rec={}
    for J,(a,Lmin,Lmax,b,pbf) in ARR.items():
        for ln,rc,rn in (('HCO',1.0,1.0),('HCN',RC[tn][J],RN[tn][J])):
            s=true_sky_r(fe,we,rc,rn); img,bm=observe(s,Lmin,Lmax,LAMS[b],pbf)
            im=np.nan_to_num(img)[sl,sl]
            sa,sb=bm; fx,fy=sa*FW,sb*FW; beam=(max(fx,fy),min(fx,fy),90.0 if fx>fy else 0.0)
            m0=im*(2*np.pi*sa*sb)         # Jy/pixel-area surface brightness [Jy km/s/arcsec^2] x beam solid angle -> Jy/beam km/s
            m0=im/PX**2*(2*np.pi*sa*sb)
            names,A=D.design(DX.ravel(),DY.ravel(),beam,D.VARS['V0']); x=D.fit(A,m0); d={n:float(x[i]) for i,n in enumerate(names)}
            d['nuc_total']=sum(d[n] for n in ('N','S','Bridge','haloA','haloB') if n in d); rec[(J,ln)]=d
            if ln=='HCN': rec[(J,'pos')]=None
        # 10" position value at E1 from the same simulated images (re-observe cheap: reuse read10)
    P2(f'\n### sky {sn}; truth ladder {tn}: clump R {RC[tn]}, nucleus R {RN[tn]}')
    res={}
    for c in ('E1','E2','W','nuc_total'):
        R={J:rec[(J,'HCN')][c]/rec[(J,'HCO')][c] for J in (2,3,4)}
        tr=RN[tn] if c=='nuc_total' else RC[tn]
        P2(f' {c:9s} recovered R21 {R[2]:.3f} R32 {R[3]:.3f} R43 {R[4]:.3f} | truth {tr[2]:.2f} {tr[3]:.2f} {tr[4]:.2f} | R32/R21 rec {R[3]/R[2]:.3f} truth {tr[3]/tr[2]:.3f}; R43/R21 rec {R[4]/R[2]:.3f} truth {tr[4]/tr[2]:.3f} | HCO+ flux rec/true J3: {rec[(3,"HCO")][c]/(CF[c] if c!="nuc_total" else CF["nuc_total"]):.2f}')
        res[c]=R
    OUT[f'{sn}|{tn}']=res
    # position values at 10" (E1) from truth-smoothed sky and observed sky
json.dump(OUT,open('/workspace/out/q38/s242.json','w'),indent=1,default=float); open('/workspace/out/q38/s242.txt','w').write('\n'.join(log2))
