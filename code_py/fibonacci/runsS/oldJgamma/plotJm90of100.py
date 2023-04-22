
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
    
outfile = "plotJm90of100.pdf"
alpha = 1/2


MS = 0.5
LW = 1
m = 100



#-- load data --

g050  = np.loadtxt("Post/s2_m90_G050_dt2_mcr.txt")
g060  = np.loadtxt("Post/s2_m90_G060_dt2_mcr.txt")
g067  = np.loadtxt("Post/s2_m90_G067_dt2_mcr.txt")
g075  = np.loadtxt("Post/s2_m90_G075_dt2_mcr.txt")
g100  = np.loadtxt("Post/s2_m90_G100_dt2_mcr.txt")


#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)


iax = ax[0]
iax.plot(x, g050[:,1], 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.50$");
iax.plot(x, g060[:,1], 'y-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.60$");
iax.plot(x, g067[:,1], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.67$");
iax.plot(x, g075[:,1], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.75$");
iax.plot(x, g100[:,1], 'm-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 1.00$");


iax = ax[1]
iax.plot(x, g050[:,4]*c, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.50$");
iax.plot(x, g060[:,4]*c, 'y-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.60$");
iax.plot(x, g067[:,4]*c, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.67$");
iax.plot(x, g075[:,4]*c, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.75$");
iax.plot(x, g100[:,4]*c, 'm-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 1.00$");
 

#-- canvas options --




iax = ax[0]  
iax.set_xlim(0,100)
iax.set_ylim(0,25)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(10, 23, "$\\alpha = 1/2$")
#iax.axhline(y=7, color='gray', lw = 0.5)


iax = ax[1]  
iax.set_xlim(0,100)
#iax.set_ylim(-0.02,0.02)
iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")
iax.axhline(y=-1, color='gray', lw = 0.5)
iax.axhline(y=-0.3732, color='gray', lw = 0.5)


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
