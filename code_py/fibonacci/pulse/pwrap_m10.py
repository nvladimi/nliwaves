
import sys
import fiboPost

outdir="Pulse"
f="s5m10"
dj = int(sys.argv[1])

fbaseout=outdir+'/'+f+'_aa2500_'
print(f, dj, fbaseout)


j = [100, 120, 140, 160, 180]
frange=(1,2501)
#frange=(1000,1400)
#frange=(1000,1100)
#frange=(1000,1010)

fiboPost.pulse2D(f, j, dj=dj, ntime=200, n=10000, frange=frange, fbaseout=fbaseout);



