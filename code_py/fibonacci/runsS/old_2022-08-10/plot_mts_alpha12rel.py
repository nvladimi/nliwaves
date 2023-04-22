
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=13)
  
    
outfile = "plot_mts_alpha12rel.pdf"

MS = 0.5
LW = 1
phi = (1 + np.sqrt(5))/2


#-- load data --


#nrange = range(12,0,-1)
nrange = (1,2,3,6,9,12)

#-- plot data --

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,2.9))


inax0 = ax[0].inset_axes([0.4, 0.39, 0.55, 0.55])
inax1 = ax[1].inset_axes([0.2, 0.39, 0.55, 0.55])


   
#---------------------------

datM  = np.loadtxt("Post/s5m90_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-1

inax0.set_prop_cycle(None)
for n in nrange:
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   inax0.plot(x,y, mfc='none', ms=MS,  lw=LW) #  label="n = {}".format(n));

inax0.plot(x,1.13*x, ":k", label="$1.13|i-d|$")

ax[0].set_prop_cycle(None)
for n in nrange:
   y = datM[:,n+1] / gamma(n/2+1) / (datM[:,2+1])**(n/2)
   ax[0].plot(x,y, mfc='none', ms=MS,  lw=LW,  label="{}".format(n));


#---------------------------

datM  = np.loadtxt("Post/s5m10_mts.txt")
M = datM[:,1:];  m = len(M)
x = m-np.arange(1,m+1)

inax1.set_prop_cycle(None)
for n in nrange:
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   inax1.plot(x,y, mfc='none', ms=MS,  lw=LW); #  label="n = {}".format(n));
inax1.plot(x,1.13*x, ':k', label="$1.13|i-d|$")

ax[1].set_prop_cycle(None)
for n in nrange:
   y = datM[:,n+1] / gamma(n/2+1) / (datM[:,2+1])**(n/2)
   ax[1].plot(x,y, mfc='none', ms=MS,  lw=LW,  label="{}".format(n));
   
         
#---------------------------



   
#-- canvas options --


for i in (0,1):
   iax = ax[i]  
   iax.set_xlim(0,200)
   iax.set_ylim(0.75,3)
   iax.set(xlabel='$|i-d|$')


for iax in (inax0,inax1):
   iax.set_xlim(0,200)
   iax.set_ylim(0,250)
   iax.set(xlabel='$|i-d|$')
   iax.legend(frameon=False, handlelength=1)


ax[0].set(ylabel='$A_m  A_2^{-m/2} / \Gamma(m/2+1)$')
ax[0].text(20, 0.84, "inverse", fontsize=12)
ax[1].text(20, 0.84, "direct", fontsize=12)
ax[1].legend(frameon=False, labelspacing = 0.35, handlelength=1,
             loc="upper right", bbox_to_anchor=(1.02, 0.90))
ax[1].text(170, 2.7, "$m$", fontsize=12)

inax0.set(ylabel='$[A_m / \Gamma(m/2+1)]^{3/m}$')


#iax.grid(axis="y")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')
#iax.set_yticks(np.arange(0, 2.5, step=0.5)) 





#-----------------------

plt.subplots_adjust(left=0.10, right=0.98, top=0.97, bottom=0.15, wspace=0.2)
#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.10, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
