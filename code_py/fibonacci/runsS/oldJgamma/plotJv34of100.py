
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
    
outfile = "plotJv34of100.pdf"
alpha = 3/4


MS = 0.5
LW = 1
m = 100



#-- load data --

g250  = np.loadtxt("Post/s2_v34_G250_dt3_mcr.txt")
g300  = np.loadtxt("Post/s2_v34_G300_dt3_mcr.txt")
g310  = np.loadtxt("Post/s2_v34_G310_dt3_mcr.txt")
g320  = np.loadtxt("Post/s2_v34_G320_dt3_mcr.txt")
g350  = np.loadtxt("Post/s2_v34_G350_dt3_mcr.txt")


#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)


iax = ax[0]
iax.plot(x, g250[:,1], 'b--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 250$");
iax.plot(x, g300[:,1], 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 300$");
iax.plot(x, g310[:,1], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 310$");
iax.plot(x, g320[:,1], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 320$");
iax.plot(x, g350[:,1], 'g--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 350$");

iax = ax[1]
iax.plot(x, g250[:,4]*c, 'b--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 250$");
iax.plot(x, g300[:,4]*c, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 300$");
iax.plot(x, g310[:,4]*c, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 310$");
iax.plot(x, g320[:,4]*c, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 320$");
iax.plot(x, g350[:,4]*c, 'g--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 350$");
 

#-- canvas options --


iax = ax[0]  
iax.set_xlim(0,100)
iax.set_ylim(3,5)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(60, 4.75, "$\\alpha = 3/4$")
iax.axhline(y=4.3, color='gray', lw = 0.5)


iax = ax[1]  
iax.set_xlim(0,100)
#iax.set_ylim(-0.02,0.02)
iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")
iax.axhline(y=1, color='gray', lw = 0.5)



#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
