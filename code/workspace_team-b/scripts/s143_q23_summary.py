# Q23: ladder-ratio shifts (T43/T32 HCN, T21/T32 HCN, HCO+ T43/T32, T32/T21) and pooled blue-vs-red line-free level test.
import json, numpy as np
R=json.load(open('/workspace/out/q23/s142.json')); L=json.load(open('/workspace/out/q23/s141.json'))
print('Fractional change of the integral when the ramp is removed (corr/raw-1 = -bias/I), median [16,84], % ; narrow = ramp ends <=~120ch, wide = extends to 160-400 ch')
for kind in ('pix','ap'):
  for scen in ('narrow','wide'):
    print(f'--- {kind} {scen}')
    for r in ('N','S','Bridge','W'):
        row=[]
        for ous,J in (('Xb55',2),('Xb4f',3),('Xb47',4)):
            o=R[f'{kind}.{scen}.{ous}.J{J}']['out'][r]
            row.append(f'J{J}: HCN {-100*o["bias_HCN"][0]/o["I_HCN"]:+.1f}% HCO+ {-100*o["bias_HCO"][0]/o["I_HCO"]:+.1f}%')
        # ladder ratio change, HCN: (I43c/I32c)/(I43/I32)-1 ≈ f43-f32 (use medians, not MC-paired: conservative upper envelope = f43 - f32 with opposite sign extremes)
        f=lambda ous,J,m: -R[f'{kind}.{scen}.{ous}.J{J}']['out'][r][f'bias_{m}'][0]/R[f'{kind}.{scen}.{ous}.J{J}']['out'][r][f'I_{m}']
        d43_32=(1+f('Xb47',4,'HCN'))/(1+f('Xb4f',3,'HCN'))-1; d21_32=(1+f('Xb55',2,'HCN'))/(1+f('Xb4f',3,'HCN'))-1
        h43_32=(1+f('Xb47',4,'HCO'))/(1+f('Xb4f',3,'HCO'))-1; h21_32=(1+f('Xb55',2,'HCO'))/(1+f('Xb4f',3,'HCO'))-1
        print(f'  {r:6s}',' | '.join(row),f'|| HCN 43/32 {100*d43_32:+.1f}%  HCN 21/32 {100*d21_32:+.1f}%  HCO+ 43/32 {100*h43_32:+.1f}%  HCO+ 21/32 {100*h21_32:+.1f}%')
print('\nHCN/HCO+ ratio shifts and their size vs total 1-sigma (stat + 2x7% cal):')
for kind in ('pix',):
  for scen in ('narrow','wide'):
    for ous,J in (('Xb55',2),('Xb4f',3),('Xb47',4)):
        o=R[f'{kind}.{scen}.{ous}.J{J}']['out']
        print(f'  {scen} J={J}: '+'  '.join(f'{r} {100*o[r]["dR_med"]:+.1f}% [{100*o[r]["dR_16"]:+.1f},{100*o[r]["dR_84"]:+.1f}] ({o[r]["dR_over_sig"]:+.2f}sig)' for r in ('N','S','Bridge')))
print('\nDirect test: mean level blue(-600..-420 km/s) minus red(420..600) in mJy/ch, pooled N,S,Bridge (peak pixel, raw spectra); expected from the ramp model ~ f*cont*shape at blue side')
for k,o in L.items():
    d=[o[f'{r}_W350']['blue']-o[f'{r}_W350']['red'] for r in ('N','S','Bridge')]; rms=np.mean([o[f'{r}_W350']['rms'] for r in ('N','S','Bridge')])
    nb=int(round(180/o['dv']))
    e=rms*1.3*np.sqrt(2/nb)
    print(f'  {k:22s} blue-red: N {d[0]:+.2f} S {d[1]:+.2f} Bridge {d[2]:+.2f}; mean {np.mean(d):+.2f} +- {e/np.sqrt(3):.2f} (rms/ch {rms:.2f}, {nb} ch each side)')
