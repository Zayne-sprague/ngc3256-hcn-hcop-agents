"""Q52: independent reproduction of team-a p0056 chemistry ratios at N and S in matched circular beams.
Lines (cube -> rest GHz): 1-0: HCN, HCO+, HNC (2015.1.00993 X961; 00412 X145 for cross-check), CN 1-0 (2018.1.00223 X4d spw27, N=1-0 hyperfine groups);
3-2: HCN, HCO+ (Xb4d), HNC (00412 X169), CN 3-2 (2018.1.00493 X1df1 spw25); HCN 4-3 (X1df1 spw29 not used here).
Method: pbcor cube -> RJ Tb (native beam) -> unit-area kernel to circular theta -> pixel at fixed N/S position (beam-weighted value) ->
linear baseline from line-free channels (|v|>600 km/s of line group edges) -> integral over +-350 km/s (radio LSRK) about v_sys (N 2785, S 2750);
noise = rms of line-free channels * sqrt(Nch) * dv. Also +-250 km/s and a systemic-fixed (2780) window for window sensitivity.
usage: python3 -I s320_q52_chem_ratios.py  -> out/q52/s320.json, s320.txt"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import load, to_tb_conv, C
from regions import REGIONS
A='/archive/files/'
B='/data/member.uid___A001_X1273_Xb4d.ngc_3256_sci.%s.cube.I.pbcor.fits'
# name: (path, [rest GHz groups (lo,hi) of line complex for integration], reference rest)
CN10=[(113.1200,113.1950)]; CN10b=[(113.4850,113.5250)]
LINES={
 '10':{ 'HCN':(A+'2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',88.6318),
        'HCO+':(A+'2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits',89.1885),
        'HNC':(A+'2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw1_89515MHz.12m.cube.I.pbcor.fits',90.6636),
        'HNC_412':(A+'2015.1.00412.S/member.uid___A001_X33e_X149.NGC_3256_sci.spw25.cube.I.pbcor.fits',90.6636),
        'HCN_412':(A+'2015.1.00412.S/member.uid___A001_X33e_X145.ari_l.NGC_3256_sci.spw1_88144MHz.12m.cube.I.pbcor.fits',88.6318),
        'HCO+_412':(A+'2015.1.00412.S/member.uid___A001_X33e_X145.ari_l.NGC_3256_sci.spw1_88144MHz.12m.cube.I.pbcor.fits',89.1885),
        'CN_G2':(A+'2018.1.00223.S/member.uid___A001_X133f_X4d.NGC3256_sci.spw27.cube.I.pbcor.fits',113.4910),
        'CN_G1':(A+'2018.1.00223.S/member.uid___A001_X133f_X4d.NGC3256_sci.spw27.cube.I.pbcor.fits',113.1441)},
 '32':{ 'HCN':(B%'spw27',265.8864),'HCO+':(B%'spw25',267.5576),
        'HNC':(A+'2015.1.00412.S/member.uid___A001_X33e_X169.ari_l.ngc_3256_sci.spw3_268777MHz.12m.cube.I.pbcor.fits',271.9811),
        'CN_G2':(A+'2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw25.cube.I.pbcor.fits',340.2478),
        'CN_G1':(A+'2018.1.00493.S/member.uid___A001_X133d_X1df1.NGC3256_sci.spw25.cube.I.pbcor.fits',340.0315)},
}
POS={'N':(REGIONS['N'][:2],2785.),'S':(REGIONS['S'][:2],2750.)}
# integration windows (km/s relative to v0, radio, w.r.t. the reference rest frequency); extra windows are masked from the baseline
WIN={'HCN':(-350,350),'HCO+':(-350,350),'HNC':(-350,350),'HNC_412':(-350,350),'HCN_412':(-350,350),'HCO+_412':(-350,350),
     'CN_G2':(-350,350),'CN_G1':(-650,350)}
WIN32={'CN_G2':(-350,350),'CN_G1':(-350,350)}   # 3-2: groups are only 190 km/s apart -> use one window on G2 reference
cache={}
def vel(nu,f0,v0): return C*(1-nu/f0)-v0
def getcube(path,theta):
    key=(path,theta)
    if key not in cache:
        cut,nu,pix,bm,pxy=load(path,half=max(8.0,3*theta))
        cache[key]=(cut,nu,pix,bm,pxy)
    return cache[key]
def measure(path,f0,theta,pos,v0,win,masks,margin=250):
    """masks: list of (f0,lo,hi) other windows to exclude from baseline"""
    cut,nu,pix,bm,pxy=getcube(path,theta)
    if theta<bm[0]*0.999: return None
    v=vel(nu,f0,v0)
    lo,hi=win
    sel=(v>lo-1400)&(v<hi+1400)
    Tb=to_tb_conv(cut[sel],nu[sel],pix,bm,theta,float(np.mean(nu)))
    vs=v[sel]; nus=nu[sel]
    x,y=pxy(pos); sp=Tb[:,int(round(y)),int(round(x))]; dv=abs(vs[1]-vs[0])
    free=(vs<lo-margin)|(vs>hi+margin)
    for (g0,glo,ghi) in masks:
        vg=vel(nus,g0,v0); free&=~((vg>glo-margin)&(vg<ghi+margin))
    free&=(vs>lo-1200)&(vs<hi+1200)
    p=np.polyfit(vs[free],sp[free],1); sp=sp-np.polyval(p,vs); sig=float(np.std(sp[free]))
    out={}
    for wn,h in (('350',1.0),('250',250/350.)):
        w=(vs>=(lo*h if lo<0 else lo))&(vs<=hi*h) if wn=='250' else (vs>=lo)&(vs<=hi)
        out[wn]=[float(sp[w].sum()*dv),float(sig*np.sqrt(w.sum())*dv)]
    w=(vs>=lo)&(vs<=hi)
    out['sig_ch_K']=sig; out['dv']=float(dv); out['Tpk']=float(sp[w].max()); out['beam']=bm; out['nfree']=int(free.sum())
    return out
def run(J,ln,th,rn):
    path,f0=LINES[J][ln]; pos,v0=POS[rn]
    base=ln.split('_')[0]+('_'+ln.split('_')[1] if ln.startswith('CN') else '')
    win=(WIN32 if J=='32' else WIN)[ln] if ln.startswith('CN') or J=='10' else (-350,350)
    masks=[]
    if J=='10' and ln.startswith('CN'):
        other='CN_G1' if ln=='CN_G2' else 'CN_G2'; masks=[(LINES[J][other][1],)+WIN[other]]
    if J=='32' and ln.startswith('CN'): win=(-350,350) if ln=='CN_G2' else None
    if J=='10' and ln in ('HCN','HCO+') : masks=[(LINES[J]['HCO+' if ln=='HCN' else 'HCN'][1],-350,350)]
    return measure(path,f0,th,pos,v0,win,masks)
if __name__=='__main__':
    res={}; txt=[]
    for J,thetas in (('10',[2.5,3.5]),('32',[1.4,2.5])):
        for th in thetas:
            for rn in POS:
                for ln in LINES[J]:
                    if J=='32' and ln=='CN_G1': continue
                    try:
                        if J=='32' and ln=='CN_G2':
                            path,f0=LINES[J][ln]; pos,v0=POS[rn]; r=measure(path,f0,th,pos,v0,(-350,450),[(f0,640,740)],margin=60)   # whole CN 3-2 complex (G2 -100..+0, G1 +190..+330 km/s); unidentified line at ~+700 masked
                            ln2='CN'
                        else: r=run(J,ln,th,rn); ln2=ln
                    except Exception as e: r=None; ln2=ln; print('ERR',J,ln,e)
                    res[f'{J}|{th}|{rn}|{ln2}']=r
                    if r: txt.append(f'{J} th={th} {rn} {ln2}: I={r["350"][0]:.1f}+-{r["350"][1]:.1f} K km/s  I(250)={r["250"][0]:.1f} sig_ch={r["sig_ch_K"]*1e3:.1f} mK Tpk={r["Tpk"]:.2f} nfree={r["nfree"]}')
    json.dump(res,open('/workspace/out/q52/s320.json','w'),indent=1)
    open('/workspace/out/q52/s320.txt','w').write('\n'.join(txt)); print('\n'.join(txt))
