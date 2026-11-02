"""Q6: assemble table (sizes in pc at 194 pc/"), beam-dilution-corrected Tb, virial masses. python3 -I s154_q6_table.py"""
import json,numpy as np
PC=194.
d=json.load(open('/workspace/out/q6/s151.json')); out={}; L=[]
bm={ln:d['lines'][ln]['_beam'] for ln in d['lines']}
L.append('region line | Ddec(arcsec) [16-84] (pc) | m0 peak snr | Tb_pk(RJ) Tb_pl | Tb_src(RJ, Gaussian-source dilution corr) | FWHM_pk(halfmax/1G) W20 | sigma_1G | Mvir')
for ln in d['lines']:
    for n,x in d['lines'][ln].items():
        if n[0]=='_': continue
        s=x['size']; sp=x['spec']['peakpix']; g=sp['gauss1']; snr=s['peak_sm']/x['sig0']
        bg=np.sqrt(bm[ln][0]*bm[ln][1]); Dd=s['Ddec']; 
        corr=1+(bg/Dd)**2 if Dd>0.5*bg else np.nan
        Tsrc=sp['Tb_RJ']*corr
        sig=g['fwhm']/2.355 if g and g['chi2r']<2 else sp['fwhm_halfmax']/2.355
        R=1.91*(Dd/2.355)*PC
        Mv=1040*sig**2*R if Dd>0 else np.nan
        out[f'{ln}_{n}']=dict(Ddec=Dd,Ddec_lo=s['Ddec_p16_84'][0],Ddec_hi=s['Ddec_p16_84'][1],Ddec_pc=Dd*PC,snr=snr,Tb=sp['Tb_RJ'],eTb=sp['eTb'],Tpl=sp['Tb_planck'],Tsrc=Tsrc,corr=corr,fwhm_hm=sp['fwhm_halfmax'],fwhm_1G=g and g['fwhm'],w20=sp['w20'],sigma=sig,Reff_pc=R,Mvir=Mv,Tb_ap=x['spec']['ap0.25']['Tb_RJ'])
        L.append(f"{n:7s}{ln:8s}| {Dd:.2f} [{s['Ddec_p16_84'][0]:.2f},{s['Ddec_p16_84'][1]:.2f}] ({Dd*PC:.0f} pc) | snr {snr:.0f} | {sp['Tb_RJ']:.2f}+-{sp['eTb']:.2f} {sp['Tb_planck']:.1f} | {Tsrc:.2f} | {sp['fwhm_halfmax']:.0f}/{g['fwhm'] if g else float('nan'):.0f} {sp['w20']:.0f} | {sig:.0f} | {Mv:.2e}")
open('/workspace/out/q6/s154_table.txt','w').write('\n'.join(L)); json.dump(out,open('/workspace/out/q6/s154.json','w'),default=float,indent=1)
print('\n'.join(L))
