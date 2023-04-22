
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  

    
outfile = "plotJJ.pdf"

MS = 0.5
LW = 1



#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(7,6.75))

phi = (1 + np.sqrt(5))/2
#f38 = phi**(1 + 3/8)


dat  = np.loadtxt("testout/s5m10_21-0_mcr.txt")

x=dat[:,0]
ax[0][0].plot(x, dat[:,1], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$C_{0-0}");
ax[0][1].plot(x, dat[:,3], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$C_{21-0}");

dat  = np.loadtxt("testout/s5m10_2211-00_mcr.txt")

x=dat[:,0]
ax[1][0].plot(x, dat[:,3], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$C_{2211-00}");
ax[1][0].plot(x, dat[:,4], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$C_{2211-00}");

dat  = np.loadtxt("testout/s5m10_210-210_mcr.txt")

x=dat[:,0]
ax[1][1].plot(x, dat[:,3], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$C_{2211-00}");
ax[1][1].plot(x, dat[:,4], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$C_{2211-00}");
ax[1][1].plot(x, (dat[:,1])**3, 'gray',  mfc='none', ms=MS,  lw=LW/2,  label="$C_{2211-00}");




#-- canvas options --




iax = ax[0][0]  
iax.set_xlim(0,200)
iax.set_ylim(0,40)
iax.set(xlabel='$i$', ylabel='$n=C_{0\\bar{0}}$')
#iax.grid(axis="y")


iax = ax[0][1]  
iax.set_xlim(0,200)
iax.set_ylim(-0.8,0)
iax.set(xlabel='$i$', ylabel='$J=C_{21\\bar{0}}$')
#iax.grid(axis="y")
iax.axhline(y=-phi**(-3/2), color='grey', lw = LW/2, linestyle='-')

iax = ax[1][0]  
iax.set_xlim(0,200)
iax.set_ylim(-400,400)
iax.set(xlabel='$i$', ylabel='$JJ=C_{2211\\bar{00}}$')
#iax.grid(axis="y")
iax.axhline(y=0, color='grey', lw = LW/2, linestyle='-')


iax = ax[1][1]  
iax.set_xlim(0,200)
iax.set_ylim(-10000,50000)
iax.set(xlabel='$i$', ylabel='$JJ^*=C_{2\\bar{2}1\\bar{1}0\\bar{0}}$')
#iax.grid(axis="y")




#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.1, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
