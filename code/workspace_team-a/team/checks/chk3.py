import csv,sys
rows=list(csv.DictReader(open('/workspace/out/splatalogue_raw.csv')))
print(len(rows),list(rows[0].keys())[:20])
def find(lo,hi):
    for r in rows:
        try:
            f=float(r.get('orderedfreq') or r.get('Freq-GHz(rest frame,redshifted)') or r.get('freq') or 0)
        except: continue
        if lo<f<hi: print(f,{k:v for k,v in r.items() if k in('name','chemical_name','resolved_QNs','upper_state_energy_K','linelist','intensity')})
find(164.9,165.2)
print('--')
find(253.85,254.15)
