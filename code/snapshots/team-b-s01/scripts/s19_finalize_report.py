# Round-3 report clean-up: (1) purge/relabel stale rows, (2) explicit-field line rows + one compact N/S table, (3) CH3OH treatment (s18), (4) literature status.
import json,re,pickle,numpy as np,pandas as pd,shutil
R='/output/report.json'; shutil.copy(R,'/workspace/out/report_before_s19.json')
d=json.load(open(R)); res=d['results']
SYS="Line centres are per region (out/s13_region_vsys.csv): N 2792+-19, S 2721+-21 km/s radio-LSRK (N-S offset +71+-28, 2.5 sigma); N-S midpoint 2756+-15 ~ nominal 2750. The 2780+-5 km/s of round 1 is the flux-weighted centre of the N+S-blended low-resolution lines, NOT a systemic velocity of either nucleus. Rest frequencies assigned at 2780 shift by <=0.9e-4 x nu (<=0.03 GHz at 265 GHz, <=25 km/s) if 2756 is used: no identification changes (blend widths 150-300 km/s)."
# (1) purge stale
res=[r for r in res if not r['quantity'].startswith('Continuum morphology: S compact')]
for r in res:
    for a in r.get('assumptions',[]):
        if a['name']=='systemic_velocity': a['value']=SYS
    if r['quantity'].startswith('Targeted lines present'):
        for a in r['assumptions']:
            if a['name']=='rest_frequencies': a['value']="Splatalogue/CDMS/JPL via astroquery (scripts/s0_linelist.py); rest = nu_obs/(1-v/c) with v=2780 (blended centre); see systemic_velocity assumption for the effect of 2756 (<=0.03 GHz)"
    if r['quantity'].startswith('Rest-frame spectral coverage'):
        r['quantity']='Rest-frame spectral coverage (computed at v=2780 km/s; shifts by <=+0.03 GHz at 2756 km/s)'
    if r['quantity'].startswith('Continuum peak (line-free'):
        r['value']+=" [Morphology note (replaces the round-1 'S compact, N extended' row, withdrawn): the native-beam S/N peak ratio 3.3 (0.24\") -> 1.2 (0.94\") -> 0.7 (5.8\") is beam dilution + resolved-out flux; at matched beams N/S = 1.4-1.6 (see matched-beam continuum entry). Other >5 sigma sources: E1, SE, bridge (see catalogue entries).]"
    if 'Michiyama+2018 abstract' in r['value']:
        r['value']=r['value'].replace("the published outflow axis/position angle was NOT verified against the papers (only abstracts available), so the geometric consistency is qualitative","Michiyama+2018 (arXiv:1810.04821, ApJ 868, 95; abstract read via arXiv API) states that the HCN/HCO+ outflow lines up with a diffuse molecular outflow from the southern nucleus and that N shows no detection of it; the outflow position angle was not found in the abstract and not read from the paper, so geometric consistency with the N-S elongation is qualitative")
    if r['quantity'].startswith('HOC+ 2-1'):
        r['value']=r['value'].split(" Comparison (literature")[0]+" Literature comparison REMOVED: the earlier numbers (HOC+/HCO+ ~0.04 'in NGC 3256, Harada+2018'; HCO+/HOC+ 30-100 in NGC 1068, Usero+2004) came from search snippets, could not be verified against the papers (ADS/arXiv searches of 2018ApJ...855...49H returned nothing usable; the HOC+ papers found by Harada+ concern NGC 253, arXiv:2109.06476), and are NOT used. Region: N+S blend (Xb55 cannot separate them; 'N'/'S' labels are only aperture centres)."
        r['target']='N+S blend (Xb55)'; r['confidence']='medium'
# (2) explicit-field rows
for r in res:
    if r['quantity']=='Line properties (see value)':
        m=re.match(r'(.+?) (\S+) \((Xb\w+), beam ([\d.x]+)", (.+)\)',r['target'])
        if m:
            r['line'],r['region'],r['ous'],r['beam_arcsec'],r['aperture']=m.groups()
            r['quantity']=f"Line properties: {r['line']} at {r['region']} ({r['ous']}, {r['beam_arcsec']}\", {r['aperture']})"
            r['convention']="peak=raw channel; 1G/2G Gaussians fitted to spectrum (1G values INVALID where chi2r>2); I=int over |v-vsys|<350 km/s, Jy/beam km/s; rms per channel in mJy/beam (native cube units)"
f=pd.read_csv('/workspace/out/s7_line_fits.csv'); t=pd.read_csv('/workspace/out/s12_line_table.csv')
want=[('Xb55','HCN 2-1'),('Xb55','HCO+ 2-1'),('Xb47','HCN 4-3'),('Xb47','HCO+ 4-3'),('Xb4d','HCN 3-2'),('Xb4d','HCO+ 3-2'),('Xb47','CS 7-6'),('Xb55','HOC+ 2-1'),('Xb55','CH3OH')]
rows=[]
for ous,ln in want:
    for reg in ('N','S'):
        tt=t[(t.ous==ous)&t.line.str.startswith(ln)&(t.region==reg)&(t.spec=='pix')]
        if not len(tt): continue
        tt=tt.iloc[0]; ff=f[(f.ous==ous)&f.line.str.startswith(ln)&(f.region==reg)&(f.spec=='pix')]
        row=dict(line=ln if ln!='CH3OH' else 'CH3OH J_-1-J_0 E (J=1-3 blend)',ous=ous,beam_arcsec=tt.beam,region=reg,aperture='peak pixel, native beam',rms_chan_mJy_per_beam=round(tt.rms_chan_mJy,2),peak_smoothed3ch_mJy_per_beam=round(tt.peak_sm3_mJy,1),I_template_Jy_per_beam_kms=round(tt.I_T,3),I_err_total=round(tt.sig_tot,3),snr=round(tt.snr,1))
        if len(ff):
            ff=ff.iloc[0]; row['peak_raw_mJy_per_beam']=round(ff.peak_raw*1e3,1)
            if ff.g1_chi2r<=2: row.update(G1_peak_mJy=round(ff.g1_A*1e3,1),G1_v0=round(ff.g1_v0,1),G1_fwhm=round(ff.g1_fwhm,1))
            else: row['G1']=f'blank (chi2r={ff.g1_chi2r:.1f}>2)'
            row.update(G2_blue=f'v={ff.g2_vblue:.0f},fwhm={ff.g2_wblue:.0f},A={ff.g2_Ablue*1e3:.1f} mJy',G2_red=f'v={ff.g2_vred:.0f},fwhm={ff.g2_wred:.0f},A={ff.g2_Ared*1e3:.1f} mJy',G2_chi2r=round(ff.g2_chi2r,2))
        rows.append(row)
pd.DataFrame(rows).to_csv('/workspace/out/s19_compact_line_table.csv',index=False)
res.append(dict(quantity="Compact line table: strongest lines at the N and S peak pixels (native beam)",target="N, S",value=json.dumps(rows),uncertainty="rms per channel and I_err_total (stat + baseline, lag-1 corrected) in the columns; add flux scale 5% (B6) / 10% (B5,B7) for absolute fluxes",unit="mJy/beam; km/s; Jy/beam km/s",
 convention="One convention: spectrum at the peak pixel of each region in the native-resolution pbcor cube (beam in column); rest frequency at the region's velocity; I_template = integral of a template (HCO+ profile of the same region/array) fitted to the line; peak_smoothed3ch is a 3-channel boxcar peak; G1 blank where chi2r>2; G2 = two-Gaussian fit (v = absolute radio-LSRK km/s); full N/S/beam-aperture/pixel table in out/s12_line_table.csv and s7_line_fits.csv; N and S are blended at Xb55 (response 0.7).",
 assumptions=[{"name":"systemic_velocity","value":SYS}],confidence="high",evidence="scripts/s12_line_table.py, s7b_fit_lines.py, s19_finalize_report.py -> out/s19_compact_line_table.csv"))
# (3) CH3OH
A=pickle.load(open('/workspace/out/s18_ch3oh_radex.pkl','rb'))['adopt']
for i,r in enumerate(res):
    if r['quantity'].startswith('CH3OH J_-1-J_0 E series blend'):
        v=r['value']; v=v.split(" Ladder check:")[0]
        sentences=("Blend integral J=1-3 (HCO+ template, +-900 km/s; scripts/s18): N %.2f+-%.2f (stat), S %.2f+-%.2f; the integral depends on template width/window: N %.2f-%.2f, S %.2f-%.2f (HCO+ and 150/220/300 km/s Gaussian templates) => ~30%% systematic (sigma/mean %.0f%%/%.0f%%), NOT included in the statistical error. J=4,5 not detected: 3-sigma J=4 <%.2f (N), <%.2f (S), J=5 <%.2f, <%.2f Jy/beam km/s. The individual J=1,2,3 amplitudes are degenerate (J=1,2 are 20 km/s apart; the unconstrained fit gives J=1<0, i.e. unphysical) and the previous Boltzmann fit (Trot=5-6 K, chi2 20-29 for 3 dof) is WITHDRAWN. "
          "Replacement: one-scale fit of LTE and RADEX (e-ch3oh.dat, LVG, dV=200, n=1e3-1e7, N=1e14-1e16, Tkin 4-150 K) amplitude patterns to the spectrum (all 5 lines, common v/shape, chi2 ~ 285-300 for ~397 channels): LTE Tex <~10-12 K (2 sigma, best 8-10 K, range over templates 6-12 K, S 6-12 K); non-LTE: n<=1e4 cm-3 excluded at any Tkin (J=4,5 would be too strong), n=1e5: Tkin 15-40 K, n>=1e6: Tkin 6-15 K. Tkin and n are degenerate; robust statement: J>=4 E-methanol lines are weak relative to J=1-3 => low excitation (Tex<~12 K) in the 9\" beam, model dependent; Tkin<~40 K for all allowed densities. Competing species at 165.05: CH3OCH3 27(1,26)-27(2,25) Eu=363 K and HC5N J=62-61 Eu=250 K are implausible.")%(A['N']['blend'],A['N']['blend_stat'],A['S']['blend'],A['S']['blend_stat'],min(A['N']['sums']),max(A['N']['sums']),min(A['S']['sums']),max(A['S']['sums']),100*A['N']['sysfrac'],100*A['S']['sysfrac'],A['N']['ul4'],A['S']['ul4'],A['N']['ul5'],A['S']['ul5'])
        r['value']=v+" "+sentences; r['evidence']+=", s18_ch3oh_radex.py -> out/s18_ch3oh_radex.txt"
        r['uncertainty']="statistical+baseline; blend-integral template/window systematic ~30% (not in stat error); cal 10%; N/S blended (beam response 0.7-0.76): apertures not independent"
        r['confidence']='medium'
# literature caveats
d['caveats']=[c for c in d['caveats'] if 'literature comparison values' not in c.lower()]
d['caveats'].append("Literature: Michiyama+2018 (arXiv:1810.04821) abstract verified via arXiv API (outflow from S nucleus; N no detection); distance 40 Mpc (Pereira-Santaella+2024) and Nishimura+2024 were taken from earlier notes and not re-verified; the HOC+ literature ratios (Harada+2018 2018ApJ...855...49H, Usero+2004) could NOT be verified and were removed. CH3OH excitation is model dependent (blend degenerate; 30% window/template systematic).")
d['results']=res; json.dump(d,open(R,'w'),indent=1)
print(len(res))
