# Unit/regression tests: the pipeline must recover the known HCO+/HCN/CS/SO lines and the CH3OH lines (identification + measurement).
import sys,subprocess,pandas as pd
ok=subprocess.run([sys.executable,'-I','/workspace/scripts/ident.py','--test']).returncode==0
t=pd.read_csv('/workspace/out/line_table_v2.csv'); c=pd.read_csv('/workspace/out/ch3oh_table.csv'); b=pd.read_csv('/workspace/out/blind_detections_reassigned.csv')
def chk(name,cond):
    global ok; print(('PASS ' if cond else 'FAIL ')+name); ok&=bool(cond)
for ln in ('HCO+ 3-2','HCN 3-2','HCO+ 4-3','HCN 4-3','HCO+ 2-1','HCN 2-1','CS 7-6'):
    chk(ln+' detected >=5 sigma at N',((t.line==ln)&(t.region=='N')&(t.snr_I>=5)).any())
chk('SO 4_5-3_4 blind feature assigned to SO',b[(b.ous=='Xb55')&b.fsky.between(176.9,177.0)].id_local.str.contains('SO 3').any())
chk('163.54 GHz assigned to CH3OH',b[(b.ous=='Xb55')&b.fsky.between(163.5,163.6)].id_local.str.startswith('CH3OH').all())
chk('251.66 GHz feature has CH3OH candidate',b[b.fsky.between(251.64,251.70)].id_local.str.contains('CH3OH').all())
chk('CH3OH Q-branch J=1-3 >=20 sigma in Xb55',((c.ous=='Xb55')&c.line.str.contains('J=1-3')&(c.snr>=20)).any())
chk('CH3OH 2(0)-1(-1) >=8 sigma at Xb4f S',((c.ous=='Xb4f')&(c.region=='S')&c.line.str.startswith('CH3OH 2(0')&(c.snr>=8)).any())
print('ALL PASS' if ok else 'FAILURES'); sys.exit(0 if ok else 1)
