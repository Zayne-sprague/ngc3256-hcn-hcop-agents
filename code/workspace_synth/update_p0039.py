import json
r=json.load(open('/output/report.json'))
T=json.load(open('/snapshots/team-a-s18/TEAM_REPORT.json'))['results']
sel=[x for x in T if x['quantity'].startswith(('Q23g','Q23h','Q30'))]
assert len(sel)==6
for x in sel:
    q=x['quantity']
    if q.startswith('Q23'):
        x['evidence']="team-a p0037 (/snapshots/team-a-s17 and s18 TEAM_REPORT.json): "+x['evidence']
    else:
        x['evidence']="team-a p0039, /snapshots/team-a-s18: "+x['evidence']+" (post: 13CO(2-1) detected at both nuclei; no published 13CO/12C/13C for NGC 3256 found by the team)"
    x['assumptions']=x.get('assumptions',[])+[{"name":"distance/scale","value":"D=40 Mpc, 1\"=213 pc (Q23); R=12C/13C 25-90 assumed, not measured"}] if q.startswith('Q30') else x.get('assumptions',[])
    r['results'].append(x)
r['interpretation']+=("\n\nCO opacity and temperature bound (team-a p0037, p0039, /snapshots/team-a-s18): CO(2-1) at 0.26\" has T_b(S)~31 K peak (20 K matched to the HCN velocities) vs N ~13 K (reproducing the 37 K maximum of Brunetti+2021 within ~10%). 13CO(2-1) is detected at both nuclei (new; no published 13CO for NGC 3256 found); matched 1.0\" 12CO/13CO is 11-15 (S 10.8-13.3, N 13.6-15.2, +-20%), nearly uniform across the field, giving tau(12CO 2-1)~3-9 for assumed 12C/13C=40-90 (S 3.8-8.8; N 2.9-6.9; only ~1.5-2 for R=25). Hence the Tk lower limits with ff=1: S >=26-37 K, N >=17-19 K; ff=0.5 raises S to 46-68 K, ff=0.3 to 73-110 K. These are LOWER LIMITS applying a 1.0\" ratio to the 0.26\" component, not Tk measurements; R and the filling factor are open. The S nucleus is at least ~2x hotter than N in CO T_b at matched resolution irrespective of tau. At 2.5\" CO(3-2)/(2-1)~0.8 and (1-0)/(2-1)~1.16 imply thermalised thick CO with Tex 12-16 K, i.e. beam-averaged cool gas that does not probe the compact component (cf. Michiyama+2018). This is consistent with, but does not discriminate, the abundance-vs-excitation question for HCN/HCO+ above.")
r['caveats'].append("CO opacity/Tk (p0037/p0039): 12C/13C is assumed (25-90), the 12CO/13CO ratio is measured at 1.0\" and applied to the 0.26\" T_b, 13CO is assumed optically thin with equal Tex and filling for both isotopologues, missing flux differs between the CO projects, a CH3OH 220.078 GHz blend in the 13CO red wing is not excluded, and C18O was not analysed; Tk(S)>=26-37 K is a conditional lower limit. Early co13.log ratios (S 10.0/12.4, N 12.7/14.7) come from an earlier sampling and are superseded by the final table (10.8/13.3, 13.6/15.2); the 'native S core' 12CO noise in co13.log is a known sampling bug and unused.")
json.dump(r,open('/output/report.json','w'),indent=1)
