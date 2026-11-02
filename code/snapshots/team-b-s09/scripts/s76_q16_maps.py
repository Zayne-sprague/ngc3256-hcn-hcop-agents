# Q16r2 (1): build and cache the moment-0 maps (|v-vsys|<W, line-free baseline) used by the deblend / injection tests.
#  Xb4f J3 native (beam B0=5.72x3.54 PA82), Xb47 J4 convolved to EXACTLY B0 (analytic Fourier kernel of the beam difference, flux rescaled by beam-area ratio),
#  Xb55 J2/J3/J4 convolved to the Xb55 native beam, Xb4d J3 native (0.94x0.74).  Also a convolution self-test: the Xb47 beam after convolution is checked on a point source.
#  usage: python3 -I s76_q16_maps.py -> out/q16b/maps.npz   (maps in Jy/beam km/s on each file's own cropped grid; wcs info stored as crpix-free arrays: ra,dec of every pixel)
import sys,numpy as np; sys.path.insert(0,'/workspace/scripts')
import s30_matchedJ_ratios as S
OUT='/workspace/out/q16b/'
def radec(wc,x0,y0,n):
    Y,X=np.mgrid[:n,:n]; ra,dec=wc.pixel_to_world_values(X+x0,Y+y0); return ra,dec
res={}
B0=None
for W in (350,):
  for mol,key in (('HCN','HCN'),('HCO+','HCOp')):
    for tag,(ous,J) in {'4f3':('Xb4f',3),'47_4':('Xb47',4),'4d3':('Xb4d',3),'55_2':('Xb55',2),'55_3':('Xb4f',3),'55_4':('Xb47',4)}.items():
        spw=S.SPW[(ous,J)][0 if mol=='HCN' else 1]
        cube,v,ps,wc,x0,y0,b0=S.load(ous,spw,S.F0[mol][J])
        if tag=='4f3': B0=b0 if mol=='HCN' else B0
        tgt=None
        if tag=='47_4': tgt=(5.72,3.54,82.0)   # Xb4f HCN native beam (HCO+ 3-2 beam is 5.73x3.54): common target
        if tag.startswith('55_') and tag!='55_2': tgt=(8.64,6.00,-90.0)
        if tag=='55_2': tgt=None
        cc=S.conv(cube,ps,b0,tgt) if tgt is not None else cube
        c,core,dv,ref=S.products(cc,v,W); m=abs(v)<W
        m0=c[m].sum(0)*dv; sig=c[ref].std(0)*np.sqrt(m.sum())*dv
        n=m0.shape[0]; ra,dec=radec(wc,x0,y0,n)
        if tag.startswith('55_') and tag in ('55_3','55_4'): pass
        res[f'{tag}|{mol}|m0']=m0; res[f'{tag}|{mol}|e']=sig; res[f'{tag}|{mol}|ra']=ra; res[f'{tag}|{mol}|dec']=dec; res[f'{tag}|{mol}|ps']=ps
        res[f'{tag}|{mol}|beam']=np.array(tgt if tgt else b0)
        print(tag,mol,m0.shape,ps,'beam',tgt if tgt else b0,'peak',m0.max(),flush=True)
# convolution self-test: delta function through the transfer function -> measure resulting beam on the coarse grid vs target
for ous,tg in (('Xb47',(5.72,3.54,82.0)),('Xb47',(8.64,6.0,-90.0))):
    cube,v,ps,wc,x0,y0,b0=S.load(ous,S.SPW[(ous,4)][0],S.F0['HCN'][4])
    n=cube.shape[1]; img=np.zeros((1,n,n),np.float32); img[0,n//2,n//2]=1.0
    # point source of native-beam-shaped Gaussian
    Y,X=np.mgrid[:n,:n]; th=np.radians(b0[2]); dx=-(X-n//2)*ps; dy=(Y-n//2)*ps   # dx = offset towards East (RA increases to decreasing array x), PA measured East of North
    u=dx*np.sin(th)+dy*np.cos(th); w_=dx*np.cos(th)-dy*np.sin(th); g=np.exp(-4*np.log(2)*((u/b0[0])**2+(w_/b0[1])**2))[None].astype(np.float32)
    out=S.conv(g,ps,b0,tg)[0]
    # second moments -> FWHM, PA
    tot=out.sum(); mx=(out*dx).sum()/tot; my=(out*dy).sum()/tot; cxx=(out*(dx-mx)**2).sum()/tot; cyy=(out*(dy-my)**2).sum()/tot; cxy=(out*(dx-mx)*(dy-my)).sum()/tot
    ev,evec=np.linalg.eigh(np.array([[cxx,cxy],[cxy,cyy]])); fw=2.3548*np.sqrt(ev[::-1]); pa=np.degrees(np.arctan2(evec[0,1],evec[1,1]))%180
    print(f'convolution self-test {ous}: native {b0[0]:.2f}x{b0[1]:.2f} PA{b0[2]:.0f} -> measured {fw[0]:.3f}x{fw[1]:.3f} PA{pa:.1f} (target {tg}); peak {out.max():.3f} (should be 1 for Jy/beam of a unit-peak source)')
np.savez(OUT+'maps.npz',**res)
