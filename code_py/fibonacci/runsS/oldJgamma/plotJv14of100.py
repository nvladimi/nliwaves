
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
    
outfile = "plotJv14of100.pdf"
alpha = 1/4


MS = 0.5
LW = 1
m = 100



#-- load data --

g015  = np.loadtxt("Post/s2_v14_G015_dt2_mcr.txt")
g029  = np.loadtxt("Post/s2_v14_G029_dt2_mcr.txt")
g030  = np.loadtxt("Post/s2_v14_G030_dt2_mcr.txt")
g031  = np.loadtxt("Post/s2_v14_G031_dt2_mcr.txt")
g060  = np.loadtxt("Post/s2_v14_G060_dt2_mcr.txt")


#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)


iax = ax[0]
iax.plot(x, g015[:,1], 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.15$");
iax.plot(x, g029[:,1], 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.29$");
iax.plot(x, g030[:,1], 'r--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.30$");
iax.plot(x, g031[:,1], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.31$");
iax.plot(x, g060[:,1], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.60$");


iax = ax[1]
iax.plot(x, g015[:,4]*c, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.15$");
iax.plot(x, g029[:,4]*c, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.29$");
iax.plot(x, g030[:,4]*c, 'r--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.30$");
iax.plot(x, g031[:,4]*c, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.31$");
iax.plot(x, g060[:,4]*c, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.60$");
 

#-- canvas options --




iax = ax[0]  
iax.set_xlim(0,100)
iax.set_ylim(0,9)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(70, 1, "$\\alpha = 1/4$")
iax.axhline(y=4.0, color='gray', lw = 0.5)


iax = ax[1]  
iax.set_xlim(0,100)
#iax.set_ylim(-0.02,0.02)
iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")
iax.axhline(y=-1, color='gray', lw = 0.5)


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
