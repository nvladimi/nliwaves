import numpy as np
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import sys
import math

c0 = sys.argv[1]     # main correlator 

rdir = "Post05"
run = "s5m90"
nmodes=200; i1=20; i2=180;
#outfile="fit_corr_1.pdf"
toff = 4

#---------------------------------------------------------------------------------


    
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6,4))
    
#-- common parameters --

MS = 1
LW = 1


#-- load and plot data --


if run in ('s5m10', 's4v34',  's4v58'):
    d = nmodes
else:
    d = 1

rname = rdir + '/' + run  + "_" + c0 + "_mcr.txt"
dat0 = np.loadtxt(rname)

# 0.i  1.nk  2.sqrt|nj|   3.Re(Q)  4.Im(Q)   5.R

x  = dat0[i1:i2,0]
xd = np.abs(x-d)

corr  = dat0[i1:i2,3]
norm  = dat0[i1:i2,5] 
y = xd*corr/norm

yavg  = np.average(y)
yrms  = np.sqrt(np.average( (y - yavg)**2 ) )

print(c0 + "  {:>+7.4f} +/- {:>6.4f} ".format(yavg, yrms)) 
     

ax.plot(x, y, 'or',   mfc='none', ms=MS, lw=LW)
#ax.plot(x, y1, ':g',  mfc='none', ms=MS, lw=LW, label="y1")
#ax.plot(x, y2, ':b',  mfc='none', ms=MS, lw=LW, label="y2")
ax.plot(x, x*0 + yavg, 'k-',  lw=LW/2)


#-----------------------

#plt.subplots_adjust(left=0.20, right=0.98, top=0.96, bottom=0.05, hspace=0.15, wspace=0.15)
#plt.savefig(outfile, pad_inches=0)


plt.ion()
plt.show()
plt.pause(toff)
plt.ioff()

#plt.close()

   



#===========================================



#=========================================================
