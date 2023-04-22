import sys
import fiboPost

fbase  = sys.argv[1]
fbaseout  = sys.argv[2]

#fiboPost.Moments(fbase, nmom=12, istart=1, iend=-1, fbaseout=fbaseout)

fiboPost.modeCorr(fbase, jj=[-2, -1], istart=1, iend=-1, convert=True, fbaseout=fbaseout)

#fiboPost.avgSpectra(fbase, fnameout=fbaseout)


#Both moments and flux take  time; separate to run in parallel.
