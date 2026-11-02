# Q22: blend integral with an empirical, continuum-scaled edge-ramp baseline (replaces the 5-baseline grid of s120b where the ramp, ch 40-110, falls inside the +-400 km/s window).
# ch = channel from the low-frequency spw22 edge (1.953 MHz = 3.58 km/s). Line centre ch~153 (163.544 GHz).
# Methods (per position, pix and ap):
#  R: ramp fitted ch 0-100 (linear, free amplitude) extrapolated, ramp clipped to <=0 beyond zero-crossing; + constant from blue side ch 230-420 -> integral +-W
#  L: straight line through anchors ch 55-100 and ch 230-420 (J=4,5 masked)
#  Q: quadratic through ch 0-100 and ch 230-420 (global)
#  Z: zero baseline (continuum-subtracted data taken at face value) -> bound
import numpy as np, json, warnings; warnings.filterwarnings('ignore')
z=np.load('/workspace/out/regspec/Xb55.spw22.cube.pbcor.npz'); nu=z['nu']/1e9 if z['nu'].max()>1e6 else z['nu']
J=json.load(open('/workspace/out/q22/s132.json')); cont=J['cont']
dv=3.5843  # km/s per channel: c*dnu/nu at 163.5 GHz
ch=np.arange(len(nu)); ok=(abs(ch-219)>36)&(abs(ch-309)>36)   # mask J=4 (165.190, ch~219) and J=5 (165.369, ch~309) +-130 km/s
VC={'N':2792,'S':2721,'Bridge':2757,'W':2757}; f0=165.061130
v=299792.458*(1-nu*1e9/(f0*1e9))  # v rel 165.061 (as s120b), increasing with ch
res={}
LINEF=[165.190475,165.369341,165.678649,166.169098]
for reg in VC:
  for st in ('pix','ap'):
    y=z[f'{reg}_{st}']*1e3  # mJy/b
    c0=int(np.argmin(abs(v-VC[reg])))
    out={}
    for W in (150,200,250,300):
        w=abs(v-VC[reg])<W; n=w.sum()
        def integ(b): return float(np.sum((y-b)[w])*dv/1e3)   # Jy/b km/s
        # R
        m=(ch<=100); p=np.polyfit(ch[m],y[m],1); ramp=np.minimum(np.polyval(p,ch),0) if p[0]>0 else np.minimum(np.polyval(p,ch),0)
        zc=-p[1]/p[0] if p[0]>0 else np.inf
        ramp=np.where(ch<zc,np.polyval(p,ch),0.) if p[0]>0 else np.zeros_like(y)
        mb=(ch>=230)&(ch<=420)&ok; cb=np.median((y-ramp)[mb])
        bR=ramp+cb
        # L
        ma=(((ch>=55)&(ch<=100))|((ch>=230)&(ch<=420)))&ok; pL=np.polyfit(ch[ma],y[ma],1); bL=np.polyval(pL,ch)
        # Q
        mq=((ch<=100)|((ch>=230)&(ch<=420)))&ok; pQ=np.polyfit(ch[mq],y[mq],2); bQ=np.polyval(pQ,ch)
        e=np.std((y-bL)[ma])*np.sqrt(n)*1.3*dv/1e3
        out[W]=dict(R=integ(bR),L=integ(bL),Q=integ(bQ),Z=integ(np.zeros_like(y)),stat=float(e),ramp_slope=float(p[0]),ramp_y0=float(p[1]),zero_ch=float(zc),blue_c=float(cb))
    res[f'{reg}_{st}']=out
    print(f'{reg:6s} {st:3s} ramp y(ch0)={out[200]["ramp_y0"]:.2f} mJy slope {out[200]["ramp_slope"]*1e3:.1f} uJy/ch zero@ch {out[200]["zero_ch"]:.0f} blue c={out[200]["blue_c"]:.2f}')
    for W in (150,200,250,300):
        o=out[W]; print(f'     +-{W}km/s: R {o["R"]:.2f} L {o["L"]:.2f} Q {o["Q"]:.2f} Z {o["Z"]:.2f}  stat+-{o["stat"]:.2f}')
json.dump(res,open('/workspace/out/q22/s134.json','w'),indent=1)
