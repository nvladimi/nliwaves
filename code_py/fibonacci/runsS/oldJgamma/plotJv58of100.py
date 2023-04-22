
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
    
outfile = "plotJv58of100.pdf"
alpha = 5/8


MS = 0.5
LW = 1
m = 100



#-- load data --

g350  = np.loadtxt("Post/s2_v58_G350_dt2_mcr.txt")
g380  = np.loadtxt("Post/s2_v58_G380_dt2_mcr.txt")
g390  = np.loadtxt("Post/s2_v58_G390_dt2_mcr.txt")
g390a = np.loadtxt("Post/s2_v58_G390_dt3_mcr.txt")
g400  = np.loadtxt("Post/s2_v58_G400_dt2_mcr.txt")
g450  = np.loadtxt("Post/s2_v58_G450_dt2_mcr.txt")



#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)


iax = ax[0]
iax.plot(x, g350[:,1], 'g--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.50$");
iax.plot(x, g380[:,1], 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.80$");
iax.plot(x, g390[:,1], 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.90$");
iax.plot(x, g390a[:,1], 'y-', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.90$");
iax.plot(x, g400[:,1], 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 4.00$");
iax.plot(x, g450[:,1], 'b--', mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 4.50$");

iax = ax[1]
#iax.plot(x, g380[:,4]*c, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.80$");
iax.plot(x, g390[:,4]*c, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.90, \\Delta t = 0.01$");
iax.plot(x, g390a[:,4]*c, 'y-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 3.90, \\Delta t = 0.001$");
#iax.plot(x, g400[:,4]*c,  'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\gamma = 4.00$");

#-- canvas options --

iax = ax[0]  
iax.set_xlim(0,100)
iax.set_ylim(6,8)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(60, 7.75, "$\\alpha = 5/8$")
iax.axhline(y=7, color='gray', lw = 0.5)


iax = ax[1]  
iax.set_xlim(0,100)
#iax.set_ylim(-0.02,0.02)
iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")
iax.axhline(y=1, color='gray', lw = 0.5)
iax.legend(frameon=False)


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
