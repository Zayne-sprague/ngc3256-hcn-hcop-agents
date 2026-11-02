"""Q18: reconcile /output/report.json with Q2r2 (idempotent). Usage: python3 -I q18_reconcile.py report.json"""
import json, re, sys
p = sys.argv[1]
d = json.load(open(p))
R = d['results']
def rep(s, old, new):
    assert old in s, old[:60]
    return s.replace(old, new)

FLUX = ("FLUX-SCALE CONVENTION (Q2r2, single convention for all wing fluxes): residual-included pbcor image fluxes (wing channels lie outside the clean mask, so Xb4d wing = residual in dirty-beam units), "
        "with a measured imaging spread (Xb4d standard vs tapered +-5-8%; Xb4b gives 0.5 (HCN) / 0.8 (HCO+) of Xb4d for extended structure) and a synthetic-uv dirty/clean bound epsilon=0.73-0.95; "
        "adopted systematic +25%/-10%. They are NOT lower limits (earlier 'lower limit' claim dropped; fluxes are more likely low by 5-25% than high).")
# entries with stale 'lower limit' wording
R[54]['value'] = rep(R[54]['value'], "numbers are lower limits (uncleaned)", "wing channels are uncleaned (outside the clean mask); flux scale per Q2r2: image fluxes with +25%/-10% systematic, not lower limits (see Q2r2 flux-scale entry)")
R[54]['assumptions'][0]['value'] = rep(R[54]['assumptions'][0]['value'], "uncleaned, dirty-beam-affected", "residual-included (uncleaned) image fluxes; flux-scale systematic +25%/-10% (Q2r2), not a lower limit")
R[58]['value'] = rep(R[58]['value'], "flux is a lower limit)", "flux scale: image flux with +25%/-10% systematic per Q2r2, not a lower limit)")
R[58]['uncertainty'] = rep(R[58]['uncertainty'], "uncleaned => missing-flux/sidelobe systematic estimated <=0.01-0.02", "uncleaned => flux-scale +25%/-10% (Q2r2: cross-imaging spread + epsilon 0.73-0.95); the earlier 'missing-flux/sidelobe <=0.01-0.02' term is superseded")
R[138]['value'] = rep(R[138]['value'], "Residual-image flux has dirty-beam units, so its flux scale carries an unknown efficiency factor (JvM epsilon)", "Residual-image flux has dirty-beam units, so its flux scale carries an efficiency factor (JvM epsilon) that Q2r2 bounded at 0.73-0.95 (synthetic uv; adopted +25%/-10%)")
R[101]['value'] += " [Q2r2: this R_wing=1.40 is for the fixed -405..-150 window; the HCN and HCO+ broad components have different kinematics (shared-kinematics rejected, dchi2=109/4 par), so the ratio is window dependent: 1.46 fixed window (Xb4d), 0.83 (Xb4b), 0.6-0.9 for the whole broad component; HDO vs HCN-rich gas undecided. Do not quote 1.40 as a gas property.]"
# mass / Mdot
R[150]['value'] = ("[MODEL-FAILURE FLAG (Q2r2 review): the single-zone LVG fit is HCO+-driven and fails the HCN 3-2/1-0 ratio (0.07 vs 0.38+-0.12) and needs N(HCN)/N(HCO+)~10; it is not a validated description. "
    "Quote the dense mass as the scaled quantity M(HCO+) = 7.8e6 x (1e-8/X(HCO+)) Msun (68% 5.2e6-1.3e7 at X=1e-8; 2.6e7 at 3e-9, 7.8e7 at 1e-9), i.e. an abundance (x10) times excitation (x3-5) estimate, 'low' confidence; HCN-based ~10x higher. "
    "The 'central' 1.4e7 used downstream is a geometric mean chosen with the CO cap and is not independent of it.] ") + R[150]['value']
R[152]['value'] = ("[RANGE FIRST (Q2r2 review): Mdot = M v_los tan i / R_proj (lobe), a scaling with an assumed inclination i=70-85 deg from the line of sight (flow near the plane of the sky; Sakamoto+2014, Pereira-Santaella+2024). "
    "With the CO-capped HCO+ mass bracket (0.8-2.6e7 Msun) Mdot ~15-60 Msun/yr for i=70-85; over the whole assumption grid 2.8-2100 Msun/yr (1.5 dex, 'low' confidence). "
    "The single 'central' case below (M=1.4e7 = geometric mean of the HCO+ mass and the CO cap, hence not independent of the CO cap; Mdot=30) must not be quoted without this range. "
    "The mass includes gas at |v|<150 km/s (mass-weighted |v|=134 km/s) that is not unambiguously outflow, whereas Sakamoto's HV mass is |v|>224 km/s, so the comparison with 50 X20 Msun/yr is order-of-magnitude only.] ") + R[152]['value']
# interpretation
t = d['interpretation']
t = rep(t, "so wing fluxes are lower limits/dirty.", "so wing fluxes are residual-included image fluxes (flux-scale systematic +25%/-10%, Q2r2; not lower limits).")
t = rep(t, "this is the only dense-gas outflow candidate in the data and supports regional comparisons for S vs its wing, but with uncleaned (lower-limit) fluxes.",
        "this is the main dense-gas outflow candidate in the data (the S blue lobe; the red feature between S and N is a second, unassigned candidate, see Q2r2) and supports regional comparisons for S vs its lobe, with residual-included image fluxes (+25%/-10% flux-scale systematic, not lower limits).")
t = rep(t, "shows the dense mass is ~1e7 Msun, comparable with the CO outflow mass of the lobe (dense fraction 0.1-1, unconstrained), and Mdot_dense ~10-60 Msun/yr for a flow within 70-85 deg of the line of sight (any value from 3 to 2000 if all assumptions are varied)",
        "gives a dense mass M = 7.8e6 (1e-8/X(HCO+)) Msun (bracket 0.8-2.6e7 with the CO cap, up to 8e7 uncapped; a factor ~10 abundance x 3-5 excitation estimate, low confidence; the single-zone model fails the HCN 3-2/1-0 ratio, so it is not a validated description, and the often-used central 1.4e7 is only a geometric mean chosen with the CO cap), comparable with the CO outflow mass of the lobe (dense fraction 0.1-1, unconstrained), and Mdot_dense ~15-60 Msun/yr for a flow 70-85 deg from the line of sight with the CO-capped mass (3-2000 over the whole grid; scaling M v tan i/R, 1.5 dex, low confidence)")
t = rep(t, "(1.5 in -405..-150, 0.6-0.9 over the whole broad component, 0.42 in the narrow component)", "(1.5 in -405..-150, 0.6-0.9 over the whole broad component, 0.42 in the narrow component; hence the Q1 'wing R=1.4' and the Q11 outer-blue excess are window values, not gas properties)")
# remove verbatim duplicated sentences in the interpretation (Q1/Q4 paragraphs were repeated)
parts = re.split(r'(?<=[.)]) (?=[A-Z])', t); seen=set(); out=[]
for s_ in parts:
    k = s_.strip()
    if len(k) > 120 and k in seen: continue
    seen.add(k); out.append(s_)
d['interpretation'] = ' '.join(out)
# caveats
c = d['caveats']
c[3] = ("Wings and the Sblue structure are real in two independent arrays (Xb4d, Xb4b convolved) but lie outside the clean mask. Single flux-scale convention (Q2r2): residual-included image fluxes with a measured imaging spread (+-5-8% Xb4d; Xb4b 0.5 of Xb4d for HCN, 0.8 for HCO+) and a synthetic-uv bias bound epsilon=0.73-0.95, adopted systematic +25%/-10%; "
        "they are NOT lower limits and the earlier sidelobe/missing-flux term of <=0.01-0.02 Jy/beam km/s is superseded. The Xb4b HCN test has no power against sidelobes on its own.")
c[25] = c[25] + " Mass and Mdot must be quoted as brackets: M=7.8e6(1e-8/X) Msun (0.8-2.6e7 CO-capped; model fails HCN 3-2/1-0; 1.4e7 is not an independent point value) and Mdot ~15-60 Msun/yr for i=70-85 (3-2000 whole grid; 30 Msun/yr is a central case of a scaling M v tan i/R). Mass includes |v|<150 km/s gas of uncertain outflow nature."
json.dump(d, open(p, 'w'), indent=1)
