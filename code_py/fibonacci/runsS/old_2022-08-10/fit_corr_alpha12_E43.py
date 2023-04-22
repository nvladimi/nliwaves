import numpy as np
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import sys
import math

c1 = sys.argv[1]
ishift = int(sys.argv[2])

#c2="11-30";  c3="21-0";
#c2="21-0";   c3="11-30";

#c2="322-0";  c3="21-0";
#c2="21-0";   c3="322-0";

#c2="431-0";  c3="21-0";
#c2="21-0";   c3="431-0";

toff = 4
rdir = "Post04"
runs = ("s5m10", )
nmodes=200; i1=20; i2=180;



#---------------------------------------------------------------------------------



combfactor = 1
for c in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A"):
    combfactor = combfactor * math.factorial(c1.count(c))

    
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,4))
    
#-- common parameters --

MS = 1
LW = 1


#-- load and plot data --

for run in runs:

    if run in ('s5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
          
    iax = ax
 
    rname = rdir + '/' + run  + "_" + c1 + "_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_" + c2 + "_mcr.txt"
    dat2 = np.loadtxt(rname)
    
    rname = "Post02/" + run  + "_" + c3 + "_mcr.txt"
    dat3 = np.loadtxt(rname)

    x  = dat1[i1:i2,0]
    sc = 1/dat1[i1:i2,2] * np.abs(x-d) / np.sqrt(combfactor)
    y1 = dat1[i1:i2,3] * sc
    
    y2 = dat2[i1:i2,3] * dat3[i1-ishift:i2-ishift,3] * sc
    y = y1-y2
    
    y1avg = np.average(y1)
    y2avg = np.average(y2)
    yavg  = np.average(y)
    yrms  = np.sqrt(np.average( (y - yavg)**2 ) )

    #print(run + " " + c1 + ":  {:>+7.4f} +/-  {:>+7.4f} ".format(yavg, yrms)) 
     
    print(c1 + ":  ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-{:>7.4f} ".format(y1avg, y2avg, yavg, yrms)) 

    label =  c1 + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'or',   mfc='none', ms=MS, lw=LW)
    iax.plot(x, y1, ':g',  mfc='none', ms=MS, lw=LW, label="y1")
    iax.plot(x, y2, ':b',  mfc='none', ms=MS, lw=LW, label="y2")
    iax.plot(x, x*0 + yavg, 'k-',  lw=LW/2)

 

#-----------------------

#plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.09, hspace=0.15, wspace=0.15)

#plt.savefig(outfile, pad_inches=0)
plt.ion()
plt.show()
plt.pause(toff)
plt.ioff()
plt.close()

   



#===========================================



#=========================================================
