# Q4g: shock/PDR/XDR-tracer ratios per region and beam, as ratios of integrated Tb (K km/s) within ONE OUS (same array, beams of the spws differ <=7%), from out/line_table_v2.csv (peak-pixel kind) and out/ch3oh_table.csv.
# value = ratio (detection >=3 sigma in the numerator) else 3-sigma upper limit; error = stat (numerator, denominator) + 5% same-OUS flux-scale systematic in quadrature.
# Output out/q4/chem_ratios.csv
import sys; sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
LT=pd.read_csv('/workspace/out/line_table_v2.csv'); LT=LT[LT.kind=='peakpix']
CH=pd.read_csv('/workspace/out/ch3oh_table.csv'); CH=CH.rename(columns={'line':'line_ch'}); CH['line']='CH3OH 2(0)-1(-1)'
PAIRS=[('CH3OH 2(0)-1(-1)','HCO+ 3-2',['Xb4d','Xb4f','Xb4b']),('HOC+ 3-2','HCO+ 3-2',['Xb4d','Xb4f','Xb4b']),
       ('HC3N 28-27','HCN 3-2',['Xb4d','Xb4f','Xb4b']),('34SO 6(6)-5(5)','HCO+ 3-2',['Xb4d','Xb4f','Xb4b']),('SiS 14-13','HCO+ 3-2',['Xb4d','Xb4f','Xb4b']),
       ('CS 7-6','HCO+ 4-3',['Xb47']),('H2CS 10(0,10)-9(0,9)','HCO+ 4-3',['Xb47']),('HC15N 4-3','HCN 4-3',['Xb47']),('SO 8_8-7_7','HCO+ 4-3',['Xb47'])]   # HCN v2=1f: see q4_vib.py (joint HCO+-wing fit; this peak-pixel window sum is not wing-corrected)
rows=[]
def get(ous,reg,line):
    if line=='CH3OH 2(0)-1(-1)':
        g=CH[(CH.ous==ous)&(CH.region==reg)]
        if g.empty: return None
        r=g.iloc[0]; K=r.I_Kkms/r.I_Jykms if abs(r.I_Jykms)>1e-6 else np.nan
        if abs(r.I_Jykms)<1e-9 or not np.isfinite(K): return None
        k=1.222e6/(r.nu_sky**2*float(r.beam.split('x')[0])*float(r.beam.split('x')[1]))   # K per Jy/beam; (I_Kkms in the table is the template-amplitude scale: use Jy -> K with the spw beam)
        return dict(I=r.I_Jykms*k,e=r.eI_Jykms*k,lim=(max(r.I_Jykms,0)+3*r.eI_Jykms)*k,snr=r.snr)
    g=LT[(LT.ous==ous)&(LT.region==reg)&(LT.line==line)]
    if g.empty: return None
    r=g.iloc[0]; return dict(I=r.I_Kkms,e=r.eI_Kkms,lim=max(r.I_Kkms,0)+3*r.eI_Kkms,snr=r.snr_I)   # limit convention: max(value,0)+3 sigma
for num,den,ouss in PAIRS:
    for ous in ouss:
        for reg in ('N','S','Bridge','Soutflow'):
            a=get(ous,reg,num); b=get(ous,reg,den)
            if a is None or b is None or not np.isfinite(a['I']): continue
            R=a['I']/b['I']; sd=np.hypot(a['e']/b['I'],R*b['e']/b['I']); sd=np.hypot(sd,0.05*abs(R))
            det=a['snr']>=3
            rows.append(dict(ous=ous,region=reg,numerator=num,denominator=den,I_num_Kkms=a['I'],eI_num=a['e'],I_den_Kkms=b['I'],snr_num=a['snr'],ratio=R if det else np.nan,eratio=sd if det else np.nan,ratio_3sig_upper=(max(a['lim'],0)/b['I']) if not det else np.nan,detected=det))
R=pd.DataFrame(rows); R.round(4).to_csv('/workspace/out/q4/chem_ratios.csv',index=False); pd.set_option('display.width',250,'display.max_rows',300)
print(R[['ous','region','numerator','denominator','I_num_Kkms','eI_num','I_den_Kkms','snr_num','ratio','eratio','ratio_3sig_upper']].round(4).to_string())
