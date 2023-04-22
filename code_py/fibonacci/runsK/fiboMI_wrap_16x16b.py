import sys
import fiboMI

fbase  = sys.argv[1]
imode  = sys.argv[2]

i = int(imode)

fbaseout = "MI_16x16b/" + fbase + "_i" + imode

fiboMI.probability4D(i, fbase, istart=1, iend=-1, binsize=0.0625, numbins=256, numbinsPhi=32, fbaseout=fbaseout)

fiboMI.MI(fbaseout)



