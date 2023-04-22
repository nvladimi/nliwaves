import sys
import fiboPX

fbase  = sys.argv[1]

#fbase = "q7_v0dt1"
#fbase = "q7b_g3500dt5"

fbaseout = "PX/" + fbase
    
fiboPX.probX(fbase, istart=1, iend=-1, binsize=0.1, numbins=200, fbaseout=fbaseout)
