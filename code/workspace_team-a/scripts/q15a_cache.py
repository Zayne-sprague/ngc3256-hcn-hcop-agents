# Q15a: cache baseline-subtracted HCN/HCO+ 3-2 cubes (1.05" circular beam, +-9" crop) for Xb4d and Xb4b; Q13a machinery (linear baseline in line-free channels).
import sys; sys.path.insert(0,'/workspace/scripts')
from q13a_xb4_compare import *
for ous in ('Xb4d','Xb4b'):
    M=wmaps(ous)
    np.savez_compressed(f'/workspace/out/q15/cube_{ous}.npz',dx=M['dx'],dy=M['dy'],dv=M['dv'],pix=M['pix'],
        HCN=M['HCN_spec'].astype(np.float32),HCO=M['HCO+_spec'].astype(np.float32),vH=M['v_HCN'],vO=M['v_HCO+'],rmsH=M['HCN_rms'],rmsO=M['HCO+_rms'])
    print(ous,M['HCN_spec'].shape,M['pix'],M['dv'])
