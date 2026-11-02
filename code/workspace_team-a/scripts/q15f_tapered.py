# Q15f: Xb4d tapered cubes (26.4/26.6 MHz = ~30 km/s channels, restricted range): ridge/lobe window fluxes, no line-free baseline possible -> compare with the Xb4d standard cube (1.05" beam) rebinned and
# measured with the same (zero) baseline offsets; the tapered range HCN 2350-3216, HCO+ 2240-3100 km/s. Window red +230..+470 truncated to the covered range.
import sys; sys.path.insert(0,'/workspace/scripts')
from q3b_lib import *
from q3_pos import POS
O='/workspace/out/q15/'; log=open(O+'q15f.log','w')
def pr(*a):
    s=' '.join(str(x) for x in a); print(s); log.write(s+'\n')
Sx,Sy=POS['S']; RIDGE=(0.3,Sy+3.3); LOBE=POS['Soutflow']; N=POS['N']
for sp in ('HCN','HCO+'):
    spw=SPW[sp]
    c=fcube('Xb4d',spw,None,9.,'cube.I.tapered.image.pbcor','_'); v=vel(c['fr'],F0[sp]); dv=abs(np.median(np.diff(v)))
    pr(f'## {sp} tapered: v {v.min():.0f}..{v.max():.0f}, dv {dv:.1f} km/s, native beam {c["beam"]}, pix {c["pix"]}')
    s={}
    for nm,(x,y) in (('ridge',RIDGE),('lobe',LOBE),('N',N),('S',POS['S']),('ctrl1',(3.5,-1.5)),('ctrl2',(-3,-6)),('ctrl3',(-3.5,3.0)),('ctrl4',(3,5))): s[nm]=aper(c,x,y,1.0)
    pr('v     '+' '.join(f'{x:6.0f}' for x in v))
    for nm in s: pr(f'{nm:6s}'+' '.join(f'{1e3*a:6.1f}' for a in s[nm]))
    # red window sum (use channels with centre in 2960-3200 or to the edge), blue mirror 2260-2500 (clipped)
    for nm in ('ridge','lobe','N','ctrl1','ctrl2','ctrl3','ctrl4'):
        r=s[nm]; red=(v>=2960)&(v<3200); blu=(v>=2375)&(v<2575)
        pr(f'  {nm}: red 2960-3200 sum {r[red].sum()*dv:+.3f} Jy/b km/s ({red.sum()} ch); blue 2375-2575 sum {r[blu].sum()*dv:+.3f}')
    # reference: Xb4d standard cube, same 1.05" (tapered beam ~0.96x0.83 ~ similar), zero baseline
    cd=vcube(sp,'Xb4d',1.05,9.); vd=cd['v']
    for nm,(x,y) in (('ridge',RIDGE),('lobe',LOBE),('ctrl1',(3.5,-1.5))):
        a=aper(cd,x,y,1.0); red=(vd>=2960)&(vd<3200); blu=(vd>=2375)&(vd<2575)
        # raw (no baseline) and with order-0 offset from the 1380-2100 mean
        off=a[(vd>=1380)&(vd<=2100)].mean()
        pr(f'  Xb4d standard cube raw {nm}: red {a[red].sum()*cd["dv"]:+.3f} (continuum-offset-removed {((a-off)[red].sum()*cd["dv"]):+.3f}); blue {a[blu].sum()*cd["dv"]:+.3f} ({((a-off)[blu].sum()*cd["dv"]):+.3f}); free-channel offset {off*1e3:+.2f} mJy/b')
