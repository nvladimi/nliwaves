import sys
import fiboPX



imode  = sys.argv[1]
i = int(imode)

fbase = "q7b_g3500dt5_s4"

fbaseout = "PX/" + fbase + "_i" + imode

# tmaxM =  4  for q7b_g3500dt6_s4.0421.ak
# tmaxM = 23  for q7b_g3500dt5_s?.2500.ak

fiboPX.corrSigma(fbase, i, jmax=7, tmaxM=12, ntaumax = 1024, istart=1, iend=-1,  fbaseout=fbaseout)

