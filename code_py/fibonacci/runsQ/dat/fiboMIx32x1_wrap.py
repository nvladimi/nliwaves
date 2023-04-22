import sys
import fiboMI

fbase  = sys.argv[1]
imode  = sys.argv[2]

i = int(imode)

fbaseout = "MIx32x1/" + fbase + "_i" + imode

fiboMI.probability4D(i, fbase, istart=1, iend=-1, binsize=1, numbins=32, numbinsPhi=32, fbaseout=fbaseout)

fiboMI.MI(fbaseout)



