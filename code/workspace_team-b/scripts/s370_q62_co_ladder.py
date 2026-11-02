"""Q62: CO / 13CO / C18O line intensities at N, S and the N-S midpoint at matched circular beams (10" where the native beam allows, 14.5" for all).
RJ T_b (rest-frame, x(1+z) not applied: obs-frame K km/s, as s53 uses rest-frame via *(1+Z): here K km/s are multiplied by (1+Z) to be rest-frame RJ).
Integrate |v|<400 km/s around VSYS 2780 (radio), linear baseline from |v|>700 km/s (non-line-contaminated windows). usage: python3 -I s370_q62_co_ladder.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
import s180_q29_load as L
C=L.C; V0=2780.; Z=1/(1-V0/C)-1
P='/archive/files/'
CUBES={ # tag: (path, rest GHz, line)
 'CO10_X77':(P+'2016.2.00094.S/member.uid___A001_X1262_X77.NGC3256_sci.spw16.cube.I.pbcor.fits',115.27120,'CO1-0'),
 'CO10_Xd5':(P+'2016.2.00042.S/member.uid___A001_X1262_Xd5.ngc3256_sci.spw16.cube.I.pbcor.fits',115.27120,'CO1-0'),
 'CO10_X4f':(P+'2018.1.00223.S/member.uid___A001_X133f_X4f.NGC3256_sci.spw16.cube.I.pbcor.fits',115.27120,'CO1-0'),
 'CO21_7m':(P+'2015.1.00714.S/member.uid___A001_X2d6_Xb2.NGC_3256_sci.spw0_227531MHz.mos.7m.cube.I.pbcor.fits',230.53800,'CO2-1'),
 'CO32_Xd9':(P+'2016.2.00042.S/member.uid___A001_X1262_Xd9.ngc3256_sci.spw16.cube.I.pbcor.fits',345.79599,'CO3-2'),
 'CO32_X67':(P+'2016.2.00094.S/member.uid___A001_X1262_X67.NGC3256_sci.spw20.cube.I.pbcor.fits',345.79599,'CO3-2'),
 'CO32_X85a':(P+'2018.1.00223.S/member.uid___A001_X133d_X85a.NGC3256_sci.spw20.cube.I.pbcor.fits',345.79599,'CO3-2'),
 '13CO21_X6b':(P+'2016.2.00094.S/member.uid___A001_X1262_X6b.ngc_3256_sci.spw16.cube.I.pbcor.fits',220.39868,'13CO2-1'),
 'C18O21_X6b':(P+'2016.2.00094.S/member.uid___A001_X1262_X6b.ngc_3256_sci.spw16.cube.I.pbcor.fits',219.56035,'C18O2-1'),
}
POS={'N':(156.96346,-43.90390),'S':(156.96340,-43.90532),'mid':(156.96343,-43.90461)}
out={}
for tag,(path,f0,ln) in CUBES.items():
    try: cut,nu,pix,bm,pxy=L.load(path,half=45.)
    except Exception as e: print(tag,'load fail',e); continue
    v=C*(1-nu/(f0/(1+Z)*1.0)*0)  # placeholder
    v=C*(1-nu/f0)            # radio velocity, observed freq vs rest
    o=np.argsort(v); v=v[o]; cut=cut[o]
    fnu=f0/(1+Z)
    for th in (10.,14.5):
        if th<=bm[0]*1.02: out[f'{tag}|{th}']='beam too large'; continue
        Tb=L.to_tb_conv(cut,nu[o],pix,bm,th,fnu)
        dv=abs(np.gradient(v)); free=(abs(v-V0)>700)&(abs(v-V0)<1500)
        win=abs(v-V0)<400
        for rn,pos in POS.items():
            x,y=pxy(pos); sp=Tb[:,int(round(y)),int(round(x))].astype(float)
            ok=free if free.sum()>10 else None
            if ok is None: out[f'{tag}|{th}|{rn}']='no baseline'; continue
            p=np.polyfit(v[ok],sp[ok],1); s=sp-np.polyval(p,v); sig=s[ok].std()
            I=(s[win]*dv[win]).sum()*(1+Z); e=sig*np.sqrt(win.sum())*np.mean(dv[win])*(1+Z)
            out[f'{tag}|{th}|{rn}']=[float(I),float(e),float(s[win].max()),float(sig)]
            print(f'{tag:11s} beam {bm[0]:.1f}x{bm[1]:.1f} th={th:5.1f} {rn:3s} W={I:8.1f}+-{e:5.1f} K km/s peakTb={s[win].max():6.2f} rms/ch={sig:.3f} K nwin={win.sum()} nfree={free.sum()}',flush=True)
json.dump(out,open('/workspace/out/q62/s370.json','w'),indent=0)
