import numpy as np
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import sys
import math

c0 = sys.argv[1]     # main correlator 
c1 = sys.argv[2]     # subtracted correlator

ishift1 = int(sys.argv[3])  #shift for subtracted correlator
ishift2 = int(sys.argv[4])  #shift for nk

multfactor = int(sys.argv[5])

toff = 4
rdir = "Post05"
run =  "s5m90"
nmodes=200; i1=20; i2=180;
outfile = "fit_tmp_E53.pdf"

#---------------------------------------------------------------------------------

    
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,4))
    
MS = 1
LW = 1


if run in ('s5m10', 's4v34',  's4v58'):
    d = nmodes
else:
    d = 1

rname = rdir + '/' + run  + "_" + c0 + "_mcr.txt"
dat0 = np.loadtxt(rname)

rname = rdir + '/' + run  + "_" + c1 + "_mcr.txt"
dat1 = np.loadtxt(rname)


# 0.i  1.nk  2.sqrt|nj|   3.Re(Q)  4.Im(Q)   5.R

x  = dat0[i1:i2,0]
xd = np.abs(x-d)

corr  = dat0[i1:i2, 3]
norm  = dat0[i1:i2, 5]
nk    = dat0[i1-ishift2 : i2-ishift2, 1]
csub  = dat1[i1-ishift1 : i2-ishift1, 3]

y1 = xd*corr/norm
y2 = xd*csub*nk/norm * multfactor

y  = y1-y2

y1avg = np.average(y1)
y2avg = np.average(y2)
yavg  = np.average(y)
yrms  = np.sqrt(np.average( (y - yavg)**2 ) )


print(c0 + "  ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-{:>7.4f} ".format(y1avg, y2avg, yavg, yrms)) 

#label =  c0 + "{:>+7.4f}".format(yavg) + "$"

ax.plot(x, y, 'or',   mfc='none', ms=MS, lw=LW, label="Diff")
ax.plot(x, y1, ':g',  mfc='none', ms=MS, lw=LW, label="Corr")
ax.plot(x, y2, ':b',  mfc='none', ms=MS, lw=LW, label="Csub")
ax.plot(x, x*0 + yavg, 'r-',  lw=LW/2)
ax.set_title(c0)
ax.set_xlabel("|i-d|")
ax.axhline(y=0, color='k', lw=LW/2, ls='--')
ax.legend(frameon=False)
 

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
