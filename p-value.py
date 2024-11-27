import math
from os import close
from scipy import stats
import sys


scr_path = sys.argv[1]
out_path = sys.argv[2]
#data = []
with open(scr_path, 'r') as f1:
    line = f1.readlines()
    for i in range(0,len(line)) :
        ls = line[i].split('\t')
        p = stats.poisson.sf(k = int(ls[0]), mu = int(ls[1]) )
        line[i] = line[i].strip() + '\t' + str(float(p)) + '\n'
     

with open(out_path, 'w') as f2:
    f2.writelines(line)
f2.close()
