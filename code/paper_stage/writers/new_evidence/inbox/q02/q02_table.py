"""q02: N and S HCN/HCO+ 3-2 separately, 12-m (Xb4d) vs 12-m+ACA (Xb4d+Xb4f), from team-a's uv-matched tclean images
(scripts/q64c_image.py, out/q64/ratios.csv). Run: python3 -I q02_table.py"""
import csv, math
R={r['tag']:r for r in csv.DictReader(open('/workspace/out/q64/ratios.csv'))}
f=lambda t,k:float(R[t][k])
def row(t,reg):
    h,o=f(t,f'{reg}_peakHCN'),f(t,f'{reg}_peakHCO'); eh,eo=f(t,'rmsHCN'),f(t,'rmsHCO')
    r=h/o; e=r*math.hypot(eh/h,eo/o); return h*1e3,o*1e3,r,e,f(t,f'{reg}_a0_0.25_R')
out=[]
print('tag set beam | reg HCN(mJy/bm) HCO+(mJy/bm) peak-ratio +-stat | r<0.25" mean ratio')
for b in ('1.0','1.5'):
    for lab,t in (('12m',f'd_u0_b{b}'),('12m+ACA',f'dc_b{b}')):
        for reg in 'NS':
            h,o,r,e,ra=row(t,reg); print(f'{lab:8s} {b}" {reg} {h:6.2f} {o:6.2f} {r:.3f} +-{e:.3f} | {ra:.3f}')
    for reg in 'NS':
        a=row(f'd_u0_b{b}',reg); c=row(f'dc_b{b}',reg)
        print(f'  {b}" {reg}: (12m+ACA)/12m ratio change {100*(c[2]/a[2]-1):+.1f}%  HCN {100*(c[0]/a[0]-1):+.1f}%  HCO+ {100*(c[1]/a[1]-1):+.1f}%')
    a=[row(f'd_u0_b{b}',x)[2] for x in 'NS']; c=[row(f'dc_b{b}',x)[2] for x in 'NS']
    print(f'  {b}" N/S: 12m {a[0]/a[1]:.3f}  12m+ACA {c[0]/c[1]:.3f}  change {100*((c[0]/c[1])/(a[0]/a[1])-1):+.1f}%')
print('\nvariants at 1.0" (peak ratio N, S, N/S): ')
for t in ('d_u0_b1.0','dc_b1.0','dB','dcB','d_w250','d_w500','eb_Xcfae6a','eb_Xda2c82','d_u25_b1.0','d_u50_b1.0'):
    n,s=row(t,'N')[2],row(t,'S')[2]; print(f'  {t:12s} N {n:.3f} S {s:.3f} N/S {n/s:.2f}')
