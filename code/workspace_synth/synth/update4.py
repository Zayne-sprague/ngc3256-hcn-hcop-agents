import json
o=json.load(open('/output/report.json'))
t=json.load(open('/snapshots/team-b-s12/TEAM_REPORT.json'))
for r in t['results'][-3:]:
    r=dict(r); r['evidence']='p0018 (team-b), /snapshots/team-b-s12/ (TEAM_REPORT.json); '+r['evidence'].split('Q20:',1)[-1].strip().join(['Q20: ','']) if False else 'p0018 (team-b), /snapshots/team-b-s12/ (TEAM_REPORT.json); '+r['evidence']
    o['results'].append(r)
o['results'][-1]['value']+=' Additional 4.5x3.0" Band 7 limits (s110e): SO 6_5-5_4 ridge <0.16 Jy/beam km/s (74x the core SO/HCN ratio), lobe <0.23: uninformative; the 1" Xb4d SO limit above is the useful one.'
o['interpretation']+='\n\n(Update, p0018.) Independent-line check of the red ridge: HCN/HCO+ 4-3 (Band 7, matched 4.5x3.0" beam, 3-2 excess as template) show only a marginal counterpart (stack 3.4 sigma, individual lines 2.3/2.9 sigma, template-dependent), with no wrong-sign signal in mirror controls; the ridge cannot be resolved from N (1.9") and S (3.3") at this beam, so Band 7 adds corroboration, not a precision measurement. The matched-beam r43/32 and HCN/HCO+ of the ridge and lobe excess are indistinguishable within 1.3 sigma, so line ratios do not separate an S counter-lobe from Bridge streaming; the ridge origin stays undecided (S counter-lobe ~55-65%). CS 7-6, SO, HC3N and H13CN limits are 4-74x above the core ratios and test nothing. The ridge mass/Mdot remain conditional on the counter-lobe hypothesis.'
o['caveats'].append('p0018: the 4-3 ridge detection is ~3 sigma and depends on the 3-2 template and taper; Band 7 beam (4.1x2.7") blends N, S and ridge; H13CN 4-3 sits on the CO 3-2 tail (negative formal excess); weak-line limits are too shallow (would need >10x lower noise). Resolving the origin needs a velocity-sign map of CO/JWST beyond 250 pc or deeper data.')
o['recipe'].append({'order':16,'stage':'verification','description':'Synthesis update 4: added team-b p0018 (Q20) results from TEAM_REPORT.json of team-b-s12; spot-checked stacks and upper limits against out/q20/s110e.txt and s110d.txt.','software':'python json','parameters':[],'scripts':['/workspace/synth/update4.py']})
json.dump(o,open('/output/report.json','w'),indent=1)
