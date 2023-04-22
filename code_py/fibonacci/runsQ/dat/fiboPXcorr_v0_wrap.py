import sys
import fiboPX


imode  = sys.argv[1]
i = int(imode)

fbase = "q7_v0dt1_s1"

fbaseout = "PX/" + fbase + "_i" + imode

# tmaxM = 14 for q7_v0dt1_s?.1500.ak
# tmaxM =  5 for q7_v0dt2_s?.0540.ak

fiboPX.corrSigma(fbase, i, jmax=7, tmaxM=12, ntaumax = 1024, istart=1, iend=-1,  fbaseout=fbaseout)

