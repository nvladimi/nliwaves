import sys
import twomode as tm
import numpy as np


fbase  = sys.argv[1]
bs     = float(sys.argv[2])
b1     = int(sys.argv[3])
b2     = int(sys.argv[4])

tm.Probability3D(fbase, istart=1, iend=-1, binsize=bs, numbins=(b1, b2, 32), fbaseout='POSTMI/3/'+fbase)
tm.MI('POSTMI/3/'+fbase)
 
