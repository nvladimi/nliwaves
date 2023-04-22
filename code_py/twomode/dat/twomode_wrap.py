import sys
import twomode as tm
import numpy as np


fbase  = sys.argv[1]
#istart = int(sys.argv[2])
#iend   = int(sys.argv[3])

#print(fbase, istart, iend)


[p1, p2,  g1, g2, dt, isave, nsave] = tm.ReadParam(fbase)


bursts = False
quad   = False

#-- pure direct cascade --

if ((p2 == 0) & (g1 == 0)):
    
    nu2  = p1 / (4*g2)
    chi  = g2*g2*g2 / p1

    if (chi < 1):
        dphi = 64
        dn   = nu2
    else:
        dphi = 64
        dn   = nu2
        quad = True

    
#-- pure inverse cascade --
    
if ((p1 == 0) & (g2 == 0)):
    
    nu1 = p2/g1
    chi = g1*g1*g1 / (2*p2)
    
    if (chi < 1):
        dphi = 64
        dn   = nu1/2
    else:
        dphi = 360
        dn   = 20*nu1
        bursts = True

#-- post-processing --
            
dat = tm.ReadData(fbase)  #, istart, iend)

tm.EvolutionNavg(dat, tavg=nsave*dt, fbaseout = "POST/" + fbase,  showplot=False)

tm.ProbabilityN(dat, numbins=100, binsize=dn, fbaseout = "POST/" + fbase,  showplot=False)

tm.ProbabilityPhi(dat, numbins=dphi, fbaseout = "POST/" + fbase,  showplot=False)

if bursts:
     tm.Bursts(dat, fbaseout="POST/" + fbase + "_factor2",  factor=2, showplot=False)
     tm.Bursts(dat, fbaseout="POST/" + fbase + "_factor10", factor=10, showplot=False)

if quad:    
    tm.ProbabilityN(dat, numbins=100, binsize=dn*4, fbaseout = "POST/" + fbase + "x4",  showplot=False)

dr = np.sqrt(dn)

if quad:
    tm.ProbabilityRho(dat, numbins=(100,100), binsize=(4*dr,dr), fbaseout = "POST/" + fbase,  showplot=False)
else:
    tm.ProbabilityRho(dat, numbins=(100,100), binsize=(dr,dr), fbaseout = "POST/" + fbase,  showplot=False)
