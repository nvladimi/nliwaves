import sys
import gfProbFlat

fdir =   sys.argv[1]
fbase =  sys.argv[2]
iz     = sys.argv[3]

iz = int(iz)


#fdir="set03"
#fbase="flat_dx2n08_dz1nz100"

fbasein = fdir + "/" + fbase 
fbaseout = 'DATA/' + fbase + '_iz' + str(iz).zfill(4)

gfProbFlat.ProbFromCount(fbasein, iz,
              nseeds=200, dI=0.125, Imax=512.0,
              dlnI=0.01, lnImin=-10, lnImax=6,
              cfold=100, fbaseout=fbaseout)

 


