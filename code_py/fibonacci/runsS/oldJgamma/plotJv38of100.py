
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
    
outfile = "plotJv38of100.pdf"
alpha = 3/8


MS = 0.5
LW = 1
m = 100



#-- load data --

g010  = np.loadtxt("Post/s2_v38_G010_dt2_mcr.txt")
g019  = np.loadtxt("Post/s2_v38_G019_dt2_mcr.txt")
g020  = np.loadtxt("Post/s2_v38_G020_dt2_mcr.txt")
g021  = np.loadtxt("Post/s2_v38_G021_dt2_mcr.txt")
g040  = np.loadtxt("Post/s2_v38_G040_dt2_mcr.txt")


#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)


iax = ax[0]
iax.plot(x, g010[:,1], 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.10$");
iax.plot(x, g019[:,1], 'r--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.19$");
iax.plot(x, g020[:,1], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.20$");
iax.plot(x, g021[:,1], 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.21$");
iax.plot(x, g040[:,1], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.40$");


iax = ax[1]
iax.plot(x, g010[:,4]*c, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.10$");
iax.plot(x, g019[:,4]*c, 'r--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.19$");
iax.plot(x, g020[:,4]*c, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.20$");
iax.plot(x, g021[:,4]*c, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.21$");
iax.plot(x, g040[:,4]*c, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 0.40$");
 

#-- canvas options --




iax = ax[0]  
iax.set_xlim(0,100)
iax.set_ylim(0,14)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(60, 2, "$\\alpha = 3/8$")
iax.axhline(y=7, color='gray', lw = 0.5)


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
