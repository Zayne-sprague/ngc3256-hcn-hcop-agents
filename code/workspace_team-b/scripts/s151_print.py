import json
d=json.load(open('/workspace/out/q6/s151.json'))
for n,s in d['cont'].items(): print('cont',n,'pk %.3g'%s['peak_sm'],'Deq %.2f Ddec %.2f [%.2f,%.2f] p95 %.2f'%(s['Deq_obs'],s['Ddec'],*s['Ddec_p16_84'],s['Ddec_p95']), s['g2'] and 'g2 core %.2fx%.2f dec %.2fx%.2f ext %.2fx%.2f A %.3g/%.3g chi %.1f'%(s['g2']['a'],s['g2']['b'],s['g2']['core_dec_maj'],s['g2']['core_dec_min'],s['g2']['aext'],s['g2']['bext'],s['g2']['Acomp'],s['g2']['Aext'],s['g2']['chi']))
for ln,r in d['lines'].items():
    for n,x in r.items():
        if n[0]=='_':continue
        s=x['size']; sp=x['spec']['peakpix']; sa=x['spec']['ap0.25']
        print(ln,n,'pk %.3g snr %.0f'%(s['peak_sm'],s['peak_sm']/x['sig0']),'Deq %.2f Ddec %.2f [%.2f,%.2f] p95 %.2f'%(s['Deq_obs'],s['Ddec'],*s['Ddec_p16_84'],s['Ddec_p95']),
          s['g2'] and 'g2 core %.2fx%.2f dec %.2fx%.2f ext %.2fx%.2f A %.3g/%.3g chi %.1f'%(s['g2']['a'],s['g2']['b'],s['g2']['core_dec_maj'],s['g2']['core_dec_min'],s['g2']['aext'],s['g2']['bext'],s['g2']['Acomp'],s['g2']['Aext'],s['g2']['chi']))
        print('     pk: peak %.1f mJy snr %.0f FWHM %.0f sig2 %.0f sigint %s Tb %.2f(+-%.2f) planck %.1f | ap: FWHM %.0f sig2 %.0f Tb %.2f | grad %s'%(sp['peak_mJy'],sp['snr'],sp['fwhm_halfmax'],sp['sigma_mom2'],sp.get('sigma_int_mom2'),sp['Tb_RJ'],sp['eTb'],sp['Tb_planck'],sa['fwhm_halfmax'],sa['sigma_mom2'],sa['Tb_RJ'],x['grad'] and '%.0f km/s/" th %.0f'%(x['grad']['G_kms_per_arcsec'],x['grad']['theta_deg'])))
