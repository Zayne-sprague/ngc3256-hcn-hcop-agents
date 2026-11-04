import json
p='/output/report.json'
d=json.load(open(p))
R=d['results']
A=json.load(open('/snapshots/team-a-s08/TEAM_REPORT.json'))['results']
# fix char-exploded caveats
cav=[];buf=''
for c in d['caveats']:
    if len(c)<=1:
        buf+=c
    else:
        if buf.strip(): cav.append(buf.strip()); 
        buf=''; cav.append(c)
if buf.strip(): cav.append(buf.strip())
d['caveats']=cav
snap='/snapshots/team-a-s08/'
def tag(i,txt):
    R[i]['value']=txt+' || '+R[i]['value']
    R[i]['evidence']+=' || Revised in light of team-a p0021 (q12c) '+snap+'out/q12c/'
tag(94,'CONDITIONAL (p0011-r02, p0021): Tk > ~56 K is the fixed-geometry one-zone result (team-a fixed-geometry N fit 62 K [39,96], P(T>=80 K)=0.28, compatible); with beam geometry measured from the 3-2 peak-Tb growth curve the one-zone N posterior moves to 15 K [12,25] (prior-edge), so the lower bound depends on the beam-filling assumption and is NOT a robust measurement.')
R[94]['confidence']='low'
tag(103,'SUPERSEDED for N temperature by p0021: N compact Tkin is unconstrained (12-96 K; 62 [39,96] fixed-geometry one-zone, 15 [12,25] measured-geometry one-zone, 27 [14,70] two-comp); the pooled median 80 K below depends on fixed geometry. The one-zone N misfit (chi2 9.1/3, p~0.03) is absorbed by beam geometry.')
R[103]['confidence']='low'
tag(104,'p0021: S one-zone adequate (chi2 5.2/3; 3.1/3 measured geometry), Tkin 19-23 K (16-84%: 14-35 K) in team-a default deblend, P(T>=80 K)<=0.06; but p0015-r01 confirms the S upper limit is deblend/aperture dependent (not tabulated over variants; Tbg>5 K and IR pumping untested). Standing value: S Tkin ~20-80 K, lower edge >~15-20 K; "cool, >=80 K excluded" NOT established.')
tag(105,'Floor test stands; S temperature statements as in the S entry (p0021).')
tag(108,'p0021 joint fit over 16 deblend variants: Delta T = T_N - T_S = +3.7 K [-10,+36], P(T_N>T_S)=0.60 (per-variant 0.42-0.81, measured geometry; 0.72-0.93 at fixed geometry); copula rho +-0.5 shifts P by <=0.1. The N-S temperature contrast is UNCONSTRAINED (size effect u=1.6 at N can absorb a T43/T32 excess). Confirms p0015 / team-b.')
R[108]['confidence']='low'
tag(75,'Bound valid only for equal Tex/filling; p0007-r02/p0021 use H13CN as a likelihood term (alpha<0.10 N, <0.15 S; 12C/13C 30-100 marginalised): N two-comp posterior survival 0.20-0.50, compact tau32 median 23-37 -> 5-14, P(tau32>10) 0.71-0.73 -> 0.38-0.54.')
tag(76,'p0007-r02/p0021: very thick, abundance-matched compact cores are disfavoured: P(both tau32>10 and abundance-matched) 0.03-0.11 at N after the H13CN term; opacity vs abundance still undecided.')
tag(102,'p0021 supersedes the "no constraint" statement for two-component models: folding the H13CN limit into the likelihood removes 50-80% of the N two-comp posterior mass (S 20-30%); it does not fix opacity vs abundance.')
tag(111,'See p0021: H13CN folded into likelihood removes 50-80% of N two-comp mass.')
new=[A[i] for i in range(109,118)]
for r in new:
    r['evidence']=r['evidence']+' || forum p0021 (team-a), replies p0007-r02, p0011-r02, p0015-r01; snapshot '+snap
R.extend(new)
d['interpretation']+="""

(Update, p0021 and replies p0007-r02, p0011-r02, p0015-r01; team-a q12c.) A second dense-gas (compact-core + envelope) component at N is neither required nor excluded. The one-zone N misfit (chi2 9.1/3, p=0.03) is acceptable (3.2/3) once beam geometry is taken from the 3-2 peak-Tb growth curve (N compact FWHM ~0.44", ~85 pc), and the Bayesian preference for two components (dlnZ +3.0 at fixed geometry) occurs in 14.5% of data-consistent one-zone simulations with wrong sizes and reverses with measured geometry. Consequently every physical scale rests on the geometry assumption: the N Tkin is unconstrained (12-96 K), team-b's Tk>~56 K (10") holds only for fixed geometry, S Tkin is ~20 K in the default deblend but deblend dependent (20-80 K), and Delta T(N-S) = +4 K [-10,+36] with P(T_N>T_S)~0.6 is unconstrained. Including the H13CN 4-3 limit as a likelihood term disfavours very thick (tau32>10), abundance-matched compact cores at N (0.03-0.11 of posterior mass) but opacity vs abundance stays open. Standing community statement: N and S differ robustly in HCN/HCO+ (0.83 vs 0.40 at 1.05"), not demonstrably in temperature or density. Conflicts resolved: where the older entries above quote N Tk >=56 K or pooled 80 K, S 'cool' or '>=80 K excluded', the q12c/p0015 values stand."""
d['caveats'].append('q12c (p0021): the 1-0 point is a literature value at 1.7" (Harada+2018) and N and S blend there; the 2-1/4-3 compact sizes are unmeasurable (Band 5/7 beams 4.5-9") and the relative size scatter (0.15, 0.30) is assumed; models are one- and two-zone LVG spheres with shared Tbg 2.73 K (5 K as test); Tbg>5 K, IR pumping and the range of S P(T>=80 K) over deblend variants were not tested.')
json.dump(d,open(p,'w'),indent=1)
