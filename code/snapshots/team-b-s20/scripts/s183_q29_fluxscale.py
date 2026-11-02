"""Q29: HCN/HCO+ 1-0 in three independent archive datasets (2016.1.00965.S X94a 0.77", 2015.1.00412.S 1.38x1.14", 2015.1.00993.S X961 2.1x1.7") at a common circular beam,
peak pixel and aperture, +-150 (Harada) and +-350 km/s windows; compares with Harada+2018 Table 13 (1.7", +-150 km/s about 2775) -> flux scale / ratio consistency.
usage: python3 -I s183_q29_fluxscale.py  -> out/q29/s183.json"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from s182_q29_measure import vel, kern_w, REG, F0, SPECIAL, xb4d
DS={'2016.1.00965 X94a':CUBES['X94a_spw27'],
    '2015.1.00412 X145':'/archive/files/2015.1.00412.S/member.uid___A001_X33e_X145.ari_l.NGC_3256_sci.spw1_88144MHz.12m.cube.I.pbcor.fits',
    '2015.1.00993 X961':'/archive/files/2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits'}
THETAS=[1.5,1.7,2.2]
HARADA={'N':(197.23,205.25,1111.22,1173.21),'S':(56.68,87.87,510.57,773.12)}   # HCN I, HCO+ I [K km/s], HCN Ipeak, HCO+ Ipeak [mK]
res={}
for name,path in DS.items():
    cut,nu,pix,bm,pxy=load(path,half=14)
    print(name,'beam %.2f x %.2f'%bm[:2],'dnu %.2f MHz'%(abs(nu[1]-nu[0])*1e3),'cut',cut.shape,flush=True)
    for th in THETAS:
        if th<max(bm[:2])*1.001: continue
        Tb=to_tb_conv(cut,nu,pix,bm,th,np.mean(nu))
        vH=vel(nu,F0['HCN']); vO=vel(nu,F0['HCOp']); free=(np.abs(vH)>600)&(np.abs(vO)>600)&(vO>-3000)&(vH<3000)
        free=(np.abs(vH)>600)&(np.abs(vO)>600)
        Tb=Tb-np.median(Tb[free],axis=0)
        ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]
        for rn in ('N','S','Bridge','Slobe','Ridge'):
            px,py=pxy(REG[rn]); dvO=abs(vO[1]-vO[0])
            k=kern_w(vO,dvO,-150,150); m0=np.tensordot(k,Tb,axes=(0,0))*dvO
            if rn in SPECIAL: mp=(px,py)
            else: mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm); mp=(float(X.flat[j]),float(Y.flat[j]))
            ix,iy=int(round(mp[0])),int(round(mp[1]))
            ap=np.hypot(X-mp[0],Y-mp[1])*pix<=(max(th,0.6) if rn in SPECIAL else th)
            o={}
            for ln,v in (('HCN',vH),('HCOp',vO)):
                dv=abs(v[1]-v[0]); r={}
                for wn,(lo,hi) in (('h150',(-150,150)),('w350',SPECIAL.get(rn,(-350,350)))):
                    k=kern_w(v,dv,lo,hi)
                    sp_p=Tb[:,iy,ix]; sp_a=Tb[:,ap].mean(axis=1)
                    for md,sp in (('pix',sp_p),('ap',sp_a)):
                        sig=np.std(sp[free]); r[wn+'_'+md]=[float((k*sp).sum()*dv),float(sig*np.sqrt((k**2).sum())*dv)]
                r['Tpk_pix']=float(np.max(np.convolve(Tb[:,iy,ix],np.ones(3)/3,'same')[np.abs(v)<300]))
                o[ln]=r
            res[f'{name}|{th}|{rn}']=o
            a,b=o['HCN'],o['HCOp']
            print(f'{name:18s} th={th} {rn:6s} pix150 HCN {a["h150_pix"][0]:6.1f}+-{a["h150_pix"][1]:.1f} HCO+ {b["h150_pix"][0]:6.1f}+-{b["h150_pix"][1]:.1f} R={a["h150_pix"][0]/b["h150_pix"][0]:.3f} | ap350 HCN {a["w350_ap"][0]:6.1f} HCO+ {b["w350_ap"][0]:6.1f} R={a["w350_ap"][0]/b["w350_ap"][0]:.3f} | Tpk HCN {a["Tpk_pix"]*1e3:.0f} HCO+ {b["Tpk_pix"]*1e3:.0f} mK',flush=True)
json.dump(res,open('/workspace/out/q29/s183.json','w'),indent=1)

res['xb4d_32']={'%.2f'%th:xb4d(th) for th in THETAS}
json.dump(res,open('/workspace/out/q29/s183.json','w'),indent=1,default=float)
