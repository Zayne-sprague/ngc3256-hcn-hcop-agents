import json
p='/output/report.json'; r=json.load(open(p))
A=[{"name":"data","value":"as Q52 (s320): 2015.1.00993 X961 (HCN, HCO+, HNC 1-0), 2018.1.00223 X4d (CN 1-0 G2/G1), Xb4d (HCN, HCO+ 3-2), 2015.1.00412 X169 (HNC 3-2); 13CO(2-1) from 2015.1.00412 X155 spw0 (0.85\")"},
 {"name":"method","value":"s320 recipe (pbcor, RJ Tb, circular-beam kernel, fixed N/S pixel, linear baseline from line-free channels), integrated in bins about v_sys (N 2785, S 2750 km/s): core |v|<100, blue -350..-100, red 100..350, wings=blue+red. Scripts s330, s331."},
 {"name":"caveat","value":"wing channels are 8% (N) and 21% (S) of HCO+ 1-0 flux, so wing ratios at N have 10-15% errors (1-0) to 50% (3-2); 13CO(2-1) is itself of moderate opacity and has no wings above noise; image-plane only (no uv matching)"}]
new=[
{"quantity":"Q53: velocity-resolved N/S of X/HCO+ (core |v|<100 vs wings 100-350 km/s), 2.5 arcsec beam, J=1-0",
 "target":"N and S, 2.5 arcsec circular beam",
 "value":"Core (|v|<100): HCN/HCO+ N 0.846+-0.005, S 0.627+-0.005, N/S 1.35+-0.01; HNC/HCO+ 0.484/0.271, N/S 1.78+-0.05; CN(G2)/HCO+ 0.892/0.474, N/S 1.88+-0.02. Full +-350 (Q52): N/S 1.56, 1.90, 2.07. Wings (blue+red): HCN/HCO+ N 1.72+-0.14, S 0.425+-0.028 (N/S 4.1+-0.4); CN(G2)/HCO+ N/S 4.4+-0.4; HNC/HCO+ N/S 2.1+-1.0 (noisy). HCO+ wings carry 7.6% of the 1-0 flux at N (11.2 of 148 K km/s) and 21% at S (22.1 of 105; red 17.5, blue 4.6). Removing the wings lowers the N/S contrast of HCN/HCO+ by 13% (1.56->1.35), HNC/HCO+ by 6%, CN/HCO+ by 9%; the contrast stays 1.35-1.9 (>30 sigma) in the core alone.",
 "unit":"ratio of K km/s","uncertainty":"stat only; core ratios +-1-3%, wings +-10% (HCN, CN) to 50% (HNC); systematics as Q52 (window, project scale cancel in N/S)","confidence":"high for core, medium for wings",
 "evidence":"scripts/s330_q53_velbins.py; out/q53/s330.txt","assumptions":A},
{"quantity":"Q53: velocity-resolved N/S of X/HCO+, J=3-2 at 1.4 and 2.5 arcsec",
 "target":"N and S, J=3-2",
 "value":"Core |v|<100, 1.4\": HCN/HCO+ N 0.784, S 0.374, N/S 2.10+-0.02 (all-velocity 2.14); HNC/HCO+ N/S 1.91+-0.08 (1.84). 2.5\": HCN/HCO+ core N/S 1.90+-0.02 (all 1.93), HNC 2.12+-0.11 (2.26). Wings: HCO+ 3-2 wing flux 2.6% of line at N and 15% at S (1.4\"); wing HCN/HCO+ N 2.6+-0.5 vs S 0.48+-0.03 (N/S ~5, 1.4\"). The 3-2 contrast is therefore unchanged (<=2%) by excluding the wings.",
 "unit":"ratio of K km/s","uncertainty":"core +-1-5%; wings 20-100% (N wings very weak in HCO+ 3-2: 0.7-2 K km/s)","confidence":"high for core, low for wings",
 "evidence":"scripts/s330_q53_velbins.py; out/q53/s330.txt","assumptions":A},
{"quantity":"Q53: X/13CO(2-1) at N and S as a column reference independent of the 1.1 mm continuum (2.5 arcsec for 1-0 lines, 1.4 and 2.5 arcsec for 3-2 lines)",
 "target":"N and S",
 "value":"13CO(2-1) integrated intensity is nearly equal at the nuclei: N 60.3+-0.7, S 57.8+-0.6 K km/s (N/S 1.04) at 2.5\"; 88.2 vs 100.7 (N/S 0.88) at 1.4\". Line/13CO N/S (2.5\", J=1-0): HCO+ 1.36+-0.03, HCN 2.12+-0.04, HNC 2.58+-0.13, CN 2.80+-0.06; (core only: 1.37, 1.84, 2.44, 2.57). J=3-2 lines/13CO(2-1) N/S at 1.4\": HCO+ 1.66, HCN 3.55+-0.09, HNC 3.05+-0.19; at 2.5\": 1.43, 2.75, 3.22. So against a column proxy HCO+ is only 1.4-1.7x higher at N, while HCN, HNC, CN are 2-3.5x higher at N (S deficit of the N-bearing group), same pattern as against the 1.1 mm continuum (Q52, 1.25-1.54 N/S).",
 "unit":"ratio","uncertainty":"stat 2-6%; 13CO(2-1) opacity/excitation (Tk higher at S would raise its S/N ratio, making N-species/13CO N/S larger still... if 13CO is thick at N the N/S of X/13CO is underestimated); 13CO(2-1) has no detectable wings, so no velocity-resolved version","confidence":"medium",
 "evidence":"scripts/s331_q53_13co_ref.py; out/q53/s331.txt","assumptions":A},
{"quantity":"Q53 verdict: depletion of N-species at S versus HCO+ boost from outflow/shock gas",
 "target":"NGC 3256 N and S",
 "value":"Mostly a deficit of HCN/HNC/CN at S (or enhancement at N), not an HCO+ boost. (1) The contrast is present in the line core (|v|<100 km/s): HCN/HCO+ N/S 1.35 (J=1-0) and 2.10 (3-2, 1.4\"), HNC 1.78/1.91, CN 1.88. (2) HCO+ wings are a larger share of the S line (21% at 1-0, 15% at 3-2) than the N line (8%, 3%), a real outflow contribution at S, but excising them removes only 13% of the HCN/HCO+ contrast at 1-0 and <=2% at 3-2 (HNC 6%, CN 9% at 1-0); an HCO+ boost accounts for at most ~15% of the observed N/S, the remaining 85-100% is in the core. (3) Against two independent references, 1.1 mm continuum (Q52) and 13CO(2-1) (this item), HCO+ is within a factor 1.4-1.7 of equal at N and S while HCN, HNC, CN are 2-3.5x higher at N. (4) In the wings the N-bearing species are not depleted at N but are at S: wing HCN/HCO+ is 1.7 at N vs 0.43 at S (N/S ~4, 1-0), i.e. the wing gas at S is also HCO+-rich/N-poor, so the S outflow contributes the same sense of chemistry as the S core, not a distinct HCO+-only component. Unresolved: composition (abundance) versus excitation/opacity within the N-bearing group; no matched-beam optically thin H13CN/H13CO+; 13CO(2-1) opacity.",
 "unit":"qualitative","uncertainty":"see components","confidence":"medium-high","evidence":"s330, s331; out/q53/","assumptions":A}]
r['results']+=new
json.dump(r,open(p,'w'),indent=1)
