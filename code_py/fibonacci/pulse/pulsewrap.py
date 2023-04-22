
import sys
import fiboPost

outdir="Pulse"
f="s5m90"
dj = int(sys.argv[1])

fbaseout=outdir+'/'+f+'_aa'
print(f, dj, fbaseout)


j = [20, 40, 60, 80, 100, 120, 140, 160, 180]
#frange=(1,5001)
frange=(1000,1400)
#frange=(1000,1100)
#frange=(1000,1010)

fiboPost.pulse2D(f, j, dj=dj, ntime=200, n=10000, frange=frange, fbaseout=fbaseout);



