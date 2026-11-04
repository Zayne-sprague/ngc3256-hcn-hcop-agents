# Free-amplitude fits of individual CH3OH lines / blends (fixed rest frequencies, common v, FWHM per region; linear baseline) in all OUS and regions.
# Output out/ch3oh_table.csv (Jy/beam km/s, K km/s, 1-sigma, 3-sigma limits), figure out/ch3oh_q_branch.png
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from ch3oh import d,toK,beam,C,c
LT=pd.read_csv('/workspace/out/line_table.csv')
def ref(ous,reg):
    ln={'Xb55':'HCO+ 2-1','Xb47':'HCO+ 4-3'}.get(ous,'HCO+ 3-2')
    r=LT[(LT.ous==ous)&(LT.region==reg)&(LT.line==ln)&(LT.kind=='peakpix')]
    return (float(r.g_v0.iloc[0]),float(r.g_fwhm.iloc[0])) if len(r) and np.isfinite(r.g_v0.iloc[0]) else (2780.,150.)
# name -> rest freqs, list of (weight at Trot) handled by explicit relative weights (LTE at 10 K for blends; only matters within a blend)
def lte_w(fr,T=10.):
    l=c[(c.f0-fr).abs()<2e-4].iloc[0]; return l.f0*l.sijmu2*np.exp(-l.Eu/T),l
Q=[(165.050175,165.061130,165.099240),(165.190475,),(165.369341,),(165.678649,),(166.169098,),(166.898566,)]
LINES={'CH3OH Q J=1-3 (-1,J)-(0,J) E blend':Q[0],'CH3OH 4(-1,3)-4(0,4) E':Q[1],'CH3OH 5(-1,4)-5(0,5) E':Q[2],'CH3OH 6(-1,5)-6(0,6) E':Q[3],
 'CH3OH 7(-1,6)-7(0,7) E':Q[4],'CH3OH 8(-1,7)-8(0,8) E':Q[5],
 'CH3OH 2(0,2)-1(-1,1) E':(254.015377,),'CH3OH 5(-2,4)-4(-1,3) E':(266.838148,),'CH3OH 6(-1,5)-5(-2,4) E':(265.289562,),'CH3OH 9(-0,9)-8(-1,7) E':(267.403471,),
 'CH3OH 7(3)-7(2) A- ':(251.641787,), 'CH3OH 3(3)-3(2) A+- (cluster 4,5,6 3_k)':(251.866524,251.890886,251.900452,251.905729,251.917065)}
def fitline(ous,reg,name,fr0s,T=10.,half=450.):
    v0,fw=ref(ous,reg); fw=min(max(fw,90),200)
    out=[]
    ws=[lte_w(f)[0] for f in fr0s]; ws=np.array(ws)/max(ws)
    for spw in [k.split('__')[1] for k in d.files if k.startswith(f'{ous}__') and k.endswith('__freq')]:
        fr=d[f'{ous}__{spw}__freq']; kk=f'{ous}__{spw}__{reg}__pk'
        if kk not in d.files: continue
        s=d[kk]; fc=fr0s[0]*(1-v0/C)
        if not (fr.min()+0.012*fc/250<fc<fr.max()-0.012*fc/250): continue   # edge guard
        bm=beam(ous,spw)
        sel=np.isfinite(s)&(np.abs((1-fr/fc)*C)<half+60)
        if sel.sum()<10: continue
        f=fr[sel]; y=toK(s[sel],f,bm); dvch=abs(np.median(np.diff(fr)))/fc*C
        # template: sum of unit-area (K km/s) Gaussians weighted ws, sum of weights normalised so that A = total integrated K km/s of the blend
        t=np.zeros_like(f)
        for fi,wi in zip(fr0s,ws/ws.sum()):
            fl=fi*(1-v0/C); sg=fw/2.3548/C*fl; t+=wi*np.exp(-0.5*((f-fl)/sg)**2)/(fw/2.3548*2.5066)
        x=(f-fc)/fc*C/100.
        M=np.vstack([t,np.ones_like(x),x]).T
        # mask other-line contamination: iterative 3-sigma clipping of the residuals outside +-1.2 FWHM
        use=np.ones(len(f),bool)
        for _ in range(4):
            p=np.linalg.lstsq(M[use],y[use],rcond=None)[0]; r=y-M@p
            far=np.abs((1-f/fc)*C-(v0-v0))>1.0*fw
            sg=1.4826*np.median(np.abs(r[far&use])) if (far&use).sum()>10 else np.std(r)
            use=np.abs(r)<4*sg
        cov=np.linalg.pinv(M[use].T@M[use])*sg**2
        A,eA=p[0],np.sqrt(cov[0,0])
        # convert K km/s -> Jy/beam km/s at line centre
        k2j=bm[0]*bm[1]*fc**2/1.222e6
        out.append(dict(ous=ous,region=reg,spw=spw,line=name,nu_sky=fc,v0=v0,fwhm=fw,I_Kkms=A,eI_Kkms=eA,I_Jykms=A*k2j,eI_Jykms=eA*k2j,snr=A/eA,rms_mK_ch=sg*1e3,
             rms_mJy_ch=sg*k2j/ (1.0) *1e3/ (1.0) / 1.0 * (1.0/ dvch) * dvch,dv=dvch,beam=f'{bm[0]:.2f}x{bm[1]:.2f}',lim3_Jykms=3*eI_ if False else 3*eA*k2j,nblend=len(fr0s)))
    return out
if __name__=='__main__':
    REG={'Xb55':['N'],'Xb47':[],'Xb4f':['N','S','Bridge','NEclump','Eclump','Soutflow'],'Xb4d':['N','S','Bridge','NEclump','Eclump','Soutflow'],'Xb4b':['N','S','Bridge']}
    rows=[]
    for ous,regs in REG.items():
        for reg in regs:
            for name,frs in LINES.items(): rows+=fitline(ous,reg,name,frs)
    df=pd.DataFrame(rows)
    # per-channel rms in mJy/beam: sg(K)*k2j
    df['rms_mJy']=df.rms_mK_ch*1e-3*df.I_Jykms/df.I_Kkms.replace(0,np.nan)*1e3
    df['status']=np.where(df.snr>=5,'detected',np.where(df.snr>=3,'marginal','limit'))
    df.drop(columns=['rms_mJy_ch']).round(4).to_csv('/workspace/out/ch3oh_table.csv',index=False)
    pd.set_option('display.width',250,'display.max_rows',500)
    print(df[['ous','region','spw','line','I_Jykms','eI_Jykms','snr','I_Kkms','v0','fwhm','rms_mJy','status']].round(3).to_string())
