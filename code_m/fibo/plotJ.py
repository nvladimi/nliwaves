
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  

    
outfile = "plotJ.pdf"

MS = 0.5
LW = 1
m = 60

x = np.arange(1,m+1)

#-- load data --

dat38  = np.loadtxt("postA_mcr.txt")


#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))


phi = (1 + np.sqrt(5))/2

f38 = phi**(1 + 3/8)


iax = ax[0]
y = dat38[:,1]
iax.plot(x,y, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

iax = ax[1]
y = dat38[:,4]*f38
iax.plot(x,y, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");


#-- canvas options --




iax = ax[0]  
#iax.set_xlim(0,60)
#iax.set_ylim(0,8)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")


iax = ax[1]  
#iax.set_xlim(0,60)
#iax.set_ylim(-0.02,0.02)
iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
