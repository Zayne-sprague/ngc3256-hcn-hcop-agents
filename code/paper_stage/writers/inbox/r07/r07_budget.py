"""R07: N/S HCN/HCO+ 3-2 contrast per beam with error budget. Inputs: out/q26/s170_*.json (0.26-1.29", 6.5"), out/q26/s400_*.json (1.0, 1.5, 2.2"; scripts/s400_r07_extra_beams.py),
out/q34/s192.json (12m+7m uv-matched 1.5", 2.2"). All: T_b ratio (RJ), aperture r = beam FWHM, window +-350 km/s about 2780 unless stated.
usage: python3 -I r07_budget.py"""
import json, numpy as np
O='/workspace/out/'
b=json.load(open(O+'q26/s170_Xb4b.json')); d=json.load(open(O+'q26/s170_Xb4d.json')); f=json.load(open(O+'q26/s170_Xb4f.json'))
b2=json.load(open(O+'q26/s400_Xb4b.json')); d2=json.load(open(O+'q26/s400_Xb4d.json')); s192=json.load(open(O+'q34/s192.json'))['res']
def con(v,k):
    n,s=v['N'][k],v['S'][k]; c=n[0]/s[0]; return c,c*np.hypot(n[1]/n[0],s[1]/s[0])
def est(v): return {k:con(v,k) for k in ('R_ap','R_pix','R_chanpeak')}
X={0.26:est(b['0.26"(50pc)']),0.2835:est(b['55pc(0.2835")']),0.6186:est(b['120pc(0.6186")']),1.0:est(b2['1.0"']),1.2887:est(b['250pc(1.2887")']),1.5:est(b2['1.5"']),2.2:est(b2['2.2"'])}
D={1.0:est(d2['1.0"']),1.2887:est(d['250pc(1.2887")']),2.2:est(d2['2.2"']),6.5:est(d['6.5"'])}
def comb(tag,w):
    r={}
    for p in 'NS':
        m=s192[f'{tag}|comb|{p}']; r[p]=(m[f'HCN32|ap|{w}'][0]/m[f'HCOp32|ap|{w}'][0], 0)
    return r['N'][0]/r['S'][0]
COMB={1.5:{'w350':comb('1p5','w350'),'h150':comb('1p5','h150')},2.2:{'w350':comb('2p2','w350'),'h150':comb('2p2','h150')}}
out=['| beam (pc) | aperture (Xb4b) | stat | pixel | peak-chan | estimator range | Xb4d aperture | 12m+7m aperture (+-350 / +-150) | uv-recovery term | total (adopted) |','|---|---|---|---|---|---|---|---|---|---|']
for th,e in X.items():
    ap,st=e['R_ap']; allv=[e[k][0] for k in e]; rng=(min(allv),max(allv))
    dx=D.get(th,{}).get('R_ap',(np.nan,))[0]; cb=COMB.get(th)
    if cb: uv=cb['w350']-ap; uvs='%.2f (to 12m+7m +-350)'%uv; uvt=abs(uv)
    elif not np.isnan(dx): uv=dx-ap; uvs='%.2f (Xb4d-Xb4b)'%uv; uvt=abs(uv)
    else: uvs='not measured; <=0.06 assumed from 1.0-1.29" Xb4d-Xb4b'; uvt=0.06
    tot=np.hypot(np.hypot(st,(rng[1]-rng[0])/2),uvt)
    out.append('| %.2f" (%d) | %.2f | %.2f | %.2f | %.2f | %.2f-%.2f | %s | %s | %s | %.2f +- %.2f (stat+est.) / -%.2f (uv, one-sided) |'%(th,round(th*194),ap,st,e['R_pix'][0],e['R_chanpeak'][0],*rng,'%.2f'%dx if not np.isnan(dx) else '-', '%.2f / %.2f'%(cb['w350'],cb['h150']) if cb else '-',uvs,ap,np.hypot(st,(rng[1]-rng[0])/2),uvt))
out.append('\n6.5" blends: Xb4d aperture %.2f, pixel %.2f, peak-chan %.2f; Xb4f aperture %.2f'%(D[6.5]['R_ap'][0],D[6.5]['R_pix'][0],D[6.5]['R_chanpeak'][0],con(f['6.5"'],'R_ap')[0]))
open('/workspace/paper_answers/r07/budget_table.md','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
