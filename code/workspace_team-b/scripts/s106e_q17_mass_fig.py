# Q17 round 2, step e: revised ridge mass / Mdot / momentum (lobe-scaled, core-null-corrected flux ratio, RADEX excitation factor from s106d) + figure of wing-vs-core.
# usage: python3 -I s106e_q17_mass_fig.py   (needs s106a, s106d outputs)
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q17r2'; txt=[]; P=lambda s:(txt.append(s),print(s))
a=json.load(open(OUT+'/s106a.json')); d=json.load(open(OUT+'/s106d.json'))
# flux ratio: core-null excess ridge / lobe excess (baseline A, excl ridge+lobe model); spread over baselines & exclusions as a systematic
rat=[]
for mol in ('HCN','HCO+'):
    for k in 'ABC':
        for ex in ('excl ridge+lobe (r<1.5")','excl ridge+lobe+N+S cores (r<1.5")'):
            r=a[f'{mol}|{k}|{ex}|red']; l=a[f'{mol}|{k}|{ex}|blue']; rat.append((mol,k,ex,r['excess']/l['excess'],r['obs']/l['obs']))
for mol in ('HCN','HCO+'):
    x=[t[3] for t in rat if t[0]==mol]; y=[t[4] for t in rat if t[0]==mol]; P(f'{mol}: ridge/lobe excess-over-core-null flux ratio {np.mean(x):.2f} (range {min(x):.2f}-{max(x):.2f} over baselines/exclusions); raw ratio {np.mean(y):.2f}')
rm=np.mean([t[3] for t in rat]); rlo=min(t[3] for t in rat); rhi=max(t[3] for t in rat)
P(f'adopted ridge/lobe flux ratio (null-corrected): {rm:.2f} (range {rlo:.2f}-{rhi:.2f}; statistical +-0.12 from matched-core scatter)')
ex=(0.8,1.45,2.4)   # mass per unit 3-2 flux ridge/lobe from RADEX r31 (16%,median,84%) s106d
Mw=(2.1e6,7.0e6)    # lobe window -360..-180 dense mass, HCO+ X=1e-8 .. 3e-9 (Q2r2)
for lab,m in (('X(HCO+)=1e-8',Mw[0]),('X(HCO+)=3e-9',Mw[1])):
    P(f'ridge mass in window ({lab}): flux-scaled only {m*rm:.2e}; with RADEX excitation factor (0.8/1.45/2.4): {m*rm*ex[0]:.2e} / {m*rm*ex[1]:.2e} / {m*rm*ex[2]:.2e} Msun')
P(f'full range of ridge window mass: {Mw[0]*rlo*ex[0]:.1e} .. {Mw[1]*rhi*ex[2]:.1e} Msun (lobe-scaled, tentative); central {Mw[0]*rm:.1e}..{Mw[1]*rm:.1e}')
D=194.; Rl=1.88*D; Rr=3.3*D; vl=253.; vr=257.; f=1.0227e-6
for i_ in (70,80,85):
    t=np.tan(np.deg2rad(i_))
    for lab,m in (('X=1e-8',Mw[0]),('X=3e-9',Mw[1])):
        Mr=m*rm; Mdl=m*vl*f*t/Rl; Mdr=Mr*vr*f*t/Rr
        P(f'i={i_} {lab}: Mdot lobe {Mdl:.1f}, ridge (flux-scaled, excitation x1) {Mdr:.1f} (x{ex[0]}..{ex[2]}: {Mdr*ex[0]:.1f}..{Mdr*ex[2]:.1f}) Msun/yr; momentum ridge {Mr*vr:.2e} (lobe {m*vl:.2e}) Msun km/s (los)')
P('Not reported: bipolar total (1+ratio) as a result; ridge counted only if counter-lobe hypothesis holds. M18 southern (red+blue) HCN-based dense mass (2.8-120)e5 Msun, Mdot 1-48 Msun/yr (v=1600 km/s, l=400 pc: different convention).')
open(OUT+'/s106e.txt','w').write('\n'.join(txt))
fig,ax=plt.subplots(1,2,figsize=(12,5))
for j,mol in enumerate(('HCN','HCO+')):
    m=np.load(f'{OUT}/s106a_maps_{mol}.npz'); sel=m['sel']&m['grid']; ok=sel&~((m['dR']<1.5)|(m['dL']<1.5))
    ax[j].scatter(m['core'][ok],m['red'][ok],s=12,c='r',alpha=.6,label='red 180-320 (control apertures)'); ax[j].scatter(m['core'][ok],m['blue'][ok],s=12,c='b',alpha=.6,label='blue -320..-180 (controls)')
    ir=m['ir']; il=m['il']; ax[j].scatter([m['core'][ir]],[m['red'][ir]],marker='*',s=250,c='r',ec='k',label='ridge'); ax[j].scatter([m['core'][il]],[m['blue'][il]],marker='*',s=250,c='b',ec='k',label='lobe')
    ax[j].set_xscale('symlog',linthresh=0.3); ax[j].set_xlabel('core flux -150..150 (Jy km/s, r=1" aperture mean)'); ax[j].set_ylabel('wing flux (Jy km/s)'); ax[j].set_title(f'{mol} 3-2 wing vs core, 1.05" beam'); ax[j].legend(fontsize=7); ax[j].axhline(0,color='k',lw=.5)
plt.tight_layout(); plt.savefig(OUT+'/s106e_wing_vs_core.png',dpi=80)
