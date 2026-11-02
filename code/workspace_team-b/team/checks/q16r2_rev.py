import sys;sys.path.insert(0,'/workspace/scripts')
from q16b_lib import *
for mol in ['HCN','HCO+']:
  for r in ['N','S','Bridge']:
    e,n=POS[r]; print(mol,r,round(ratio43(float(sample('47_4',mol,e,n)),float(sample('4f3',mol,e,n)),mol),3))
