import json
p='/output/report.json'; r=json.load(open(p)); R=r['results']
def find(s): 
    for x in R:
        if x['quantity'].startswith(s): return x
    raise KeyError(s)
TAG='[Q56 REVISED: the 0.306/0.217 HCN R32/10 and the "30% lower HCN excitation at S" below used the Xb4d 12m-only image-plane 3-2, which loses 25% (HCN) / 15% (HCO+) of the S flux at 2.5" (comb/Xb4d = 1.25/1.16, N 1.00/0.98). With uv-complete 12m+7m 3-2 (s342) the S/N of R32/10 is 0.89 (HCN) and 1.05 (HCO+) at 2.5", |v|<100; the project statement of equal HCN 3-2/1-0 at N and S (Q34, 2.2": 1.007+-0.020) stands to within ~10%.] '
e=find('Q53: opacity/excitation checks')
e['value']=TAG+e['value']; e['confidence']='low'
e=find('Q53 verdict')
e['value']=e['value'].replace("HCN 3-2/1-0 is 30% lower at S (HCO+ 11%): a lower excitation of HCN at S contributes, so 'depletion' (abundance) is not demonstrated","[Q56: the claimed 30% lower HCN 3-2/1-0 at S was an Xb4d uv-incompleteness artefact; with uv-complete 3-2, HCN R32/10 S/N = 0.87-0.96 and HCO+ 0.98-1.05, i.e. at most a 5-13% HCN excitation difference]; abundance vs opacity is therefore not separated by the 3-2/1-0 ratio, which is nearly equal at N and S")
e['assumptions']=e['assumptions'] if isinstance(e['assumptions'],list) else e['assumptions']
e=find('Q35: is the N-vs-S CS excitation')
e['value']='[Q56 note: the HCN and HCO+ 3-2/1-0 N/S values quoted here (1.16, 0.90) are image-plane from Harada 1-0 at 1.7" with uv-incomplete Xb4d 3-2; uv-complete matched-beam values are HCN 1.007+-0.020, HCO+ 1.025+-0.010 (Q34) and 0.89-0.96 / 1.0 at 2.5" (Q56); the CS-vs-HCN/HCO+ excitation contrast is therefore at most weakly present.] '+e['value']
A=[{"name":"data","value":"3-2: 2017.1.00379.S Xb4d 12m EB + Xb4f 7m EBs restored from visibilities, tclean mosaic, uvtaper+restoringbeam circular 2.5/3.5\" (s341 -> out/q34/*_2p5_comb, *_3p5_comb), 20 km/s channels; Xb4d archive image convolved for comparison; 1-0: 2015.1.00993 (X961, 12m+ACA)"},
 {"name":"method","value":"RJ Tb, fixed N/S pixel (s320 POS), linear/median baseline from line-free channels, windows |v|<100/150/350 about N 2785, S 2750 km/s; comb cube baseline from |v|>380 (cube range -600..+550); HCN 1-0 hyperfine satellites (-7,+5 km/s) lie inside the windows; 3-2 lines fully in band"}]
new=[
{"quantity":"Q56: R32/10 = Tb(3-2)/Tb(1-0) for HCN and HCO+ at N and S, uv-complete 3-2 (12m+7m) vs Xb4d image-plane, matched circular 2.5\" and 3.5\" beams",
 "target":"N and S, 2.5\" and 3.5\", core |v|<100 (|v|<150, <350 in brackets)",
 "value":"uv-complete 2.5\": HCN N 0.304 (0.303; 0.293), S 0.271 (0.284; 0.283), S/N 0.89 (0.94; 0.96); HCO+ N 0.396 (0.395; 0.391), S 0.415 (0.413; 0.406), S/N 1.05 (1.04; 1.04). 3.5\": HCN N 0.283, S 0.246, S/N 0.87 (0.90; 0.92); HCO+ N 0.372, S 0.368, S/N 0.99 (0.98; 0.98). Xb4d image-plane (Q53) gave HCN S/N 0.71 / 0.69 and HCO+ 0.89 / 0.84: the Xb4d 3-2 at S is 25-32% (HCN) and 14-21% (HCO+) below the uv-complete flux, at N only 0-5%. Reconciliation: the Q53 '30% lower HCN 3-2/1-0 at S' is a missing-flux artefact (more extended 3-2 emission at S; same effect as Q34 recovery); the earlier 'equal HCN R32/10 at N and S' (Q34: HCN N 0.263, S 0.261 at 2.2\", all velocities) is right to ~10%. Residual HCN S/N 0.87-0.96 (pixel vs aperture, window dependent) is 1-2 sigma of systematics; HCN-HCO+ difference of S/N (0.89-0.96 vs 1.0-1.05) is a <=10-15% hint of lower HCN 3-2 excitation relative to HCO+ at S, not a 30% effect.",
 "unit":"ratio of K km/s","uncertainty":"stat +-0.002-0.005; window/aperture systematic on S/N ~+-0.05; 3-2 uv-complete imaging with auto-multithresh, flux scale cancels in S/N","confidence":"medium-high that the 30% deficit is an artefact; medium for the residual 5-13%",
 "evidence":"scripts/s340_q56_sweep.py, s342_q56_comb.py, s341_q56_casa_image.py; out/q56/s340.txt, s342.txt","assumptions":A},
{"quantity":"Q56: N/S HCN/HCO+ contrast at 1-0 and at uv-complete 3-2, matched 2.5\" and 3.5\" beams (D = contrast32/contrast10)",
 "target":"N/S, core windows",
 "value":"2.5\": 1-0 1.48 (|v|<100), 1.55 (150), 1.56 (350); 3-2 uv-complete 1.74, 1.73, 1.68 (Xb4d image-plane 1.85, 1.83, 1.78); D = 1.18, 1.12, 1.08. 3.5\": 1-0 1.41, 1.47, 1.48; 3-2 1.61, 1.60, 1.57 (Xb4d 1.72, 1.70, 1.66); D = 1.14, 1.09, 1.06. So the 3-2 contrast is 1.6-1.7 and the 1-0 contrast 1.4-1.55: D = 1.06-1.18 depending on window (Q34 aperture-mean 2.2\": 0.98-1.04). The image-plane 1.85 is inflated by ~6-8% by the missing S flux. Excitation alone (HCN R32/10 down by 5-13% at S relative to HCO+) can explain at most ~10% of the contrast; ~90% needs abundance/opacity/different columns.",
 "unit":"ratio","uncertainty":"+-0.03 stat+window","confidence":"medium","evidence":"out/q56/s342.txt","assumptions":A}]
R+=new; json.dump(r,open(p,'w'),indent=1)
