
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
    
outfile = "plotJm10.pdf"
alpha = 1/2


MS = 0.5
LW = 1
m = 60



#-- load data --

datP1dt1  = np.loadtxt("Post/s1_m10_P1_dt1_mcr.txt")
datP1dt2  = np.loadtxt("Post/s1_m10_P1_dt2_mcr.txt")
datP1dt3  = np.loadtxt("Post/s1_m10_P1_dt3_mcr.txt")
datPtdt2  = np.loadtxt("Post/s1_m10_Pt_dt2_mcr.txt")


#-- plot data --


fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)


iax = ax[0]
y = datP1dt1[:,1]
iax.plot(x,y, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\Delta t = 0.1$");
y = datP1dt2[:,1]
iax.plot(x,y, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\Delta t = 0.01$");
y = datP1dt3[:,1]
iax.plot(x,y, 'y-',  mfc='none', ms=MS,  lw=LW,  label="$\\Delta t = 0.001$");
y = datPtdt2[:,1]
iax.plot(x,y, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$c_P = 0.3732$");



iax = ax[1]
y = datP1dt1[:,4]*c
iax.plot(x,y, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\Delta t = 0.1$");
y = datP1dt2[:,4]*c
iax.plot(x,y, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\delta t = 0.01$");
y = datP1dt3[:,4]*c
iax.plot(x,y, 'y-',  mfc='none', ms=MS,  lw=LW,  label="$\\delta t = 0.001$");
y = datPtdt2[:,4]*c
iax.plot(x,y, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\delta t = 0.01$");


#-- canvas options --




iax = ax[0]  
iax.set_xlim(0,60)
iax.set_ylim(0,17.5)
iax.set(xlabel='$i$', ylabel='$n=C_0$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(40, 1, "$\\alpha = 1/2$")
#iax.axhline(y=7, color='gray', lw = 0.5)


iax = ax[1]  
iax.set_xlim(0,60)
#iax.set_ylim(-0.02,0.02)
iax.set(xlabel='$i$', ylabel='$J$')
#iax.grid(axis="y")
iax.axhline(y=1, color='gray', lw = 0.5)
iax.axhline(y=0.3732, color='gray', lw = 0.5)


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.10, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
