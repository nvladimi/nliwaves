import sys
import gfProbFlat
import numpy as np

fdir = "set21"
fbase = "fr51_dx1n08_dz1nz200"
outdirRebin = "DATA_rebin21"
outdirTrim  = "DATA_trim21"
seeds = np.arange(40)

izall = (5, 8, 10, 12,   14, 16, 18, 20,   24, 28, 30, 34,
    40, 50, 60, 70,    80, 90, 100, 120,   140, 160, 180)

Iranges = ( (0.125, 4), (0.25, 8), (0.5, 16), (1, 24), (2, 32), (4, 1000) )


fbasein = fdir + "/" + fbase 

for iz in izall:

    fbaseout = outdirRebin + '/' + fbase + '_iz' + str(iz).zfill(4)

    gfProbFlat.ProbFromCount(fbasein, iz, seeds, 
              dI=1e-3, Imax=100.0,
              dlnI=0.01, lnImin=-10, lnImax=6,
              ranges=Iranges, fbaseout=fbaseout)


    fbaseout = outdirTrim + '/' + fbase + '_iz' + str(iz).zfill(4)

    gfProbFlat.ProbFromCount(fbasein, iz, seeds,
              dI=1.0e-3, Imax=100.0,
              dlnI=0.01, lnImin=-10, lnImax=6,
              fbaseout=fbaseout)



