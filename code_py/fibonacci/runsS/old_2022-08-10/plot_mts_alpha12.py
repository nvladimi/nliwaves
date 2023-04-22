
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 
  
    
outfile = "plot_mts_alpha12.pdf"

MS = 0.5
LW = 1
phi = (1 + np.sqrt(5))/2


#-- load data --


#nrange = range(12,0,-1)
nrange = (1,2,3,6,9,12)

#-- plot data --

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6,5.8))



   
#---------------------------

datM  = np.loadtxt("Post/s5m90_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-1
ax[1][0].set_prop_cycle(None)
for n in nrange:
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   ax[1][0].plot(x,y, mfc='none', ms=MS,  lw=LW) #  label="n = {}".format(n));

ax[1][0].plot(x,1.13*x, ":k", label="$1.13|i-d|$")
   
x = np.arange(1,m+1)-1
ax[0][0].set_prop_cycle(None)
for n in nrange:
   y = ( np.log(datM[:,n+1]) - np.log(gamma(n/2+1)) ) / n
   ax[0][0].plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));



#---------------------------

datM  = np.loadtxt("Post/s5m10_mts.txt")
M = datM[:,1:];  m = len(M)
x = m-np.arange(1,m+1)
ax[1][1].set_prop_cycle(None)
for n in nrange:
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   ax[1][1].plot(x,y, mfc='none', ms=MS,  lw=LW); #  label="n = {}".format(n));
ax[1][1].plot(x,1.13*x, ':k', label="$1.13|i-d|$")

   
x = np.arange(1,m+1)-1
ax[0][1].set_prop_cycle(None)
for n in nrange:
   y = ( np.log(datM[:,n+1]) - np.log(gamma(n/2+1)) ) / n
   ax[0][1].plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

   
         
#---------------------------



   
#-- canvas options --


for i in (0,1):
   iax = ax[0][i]  
   iax.set_xlim(0,200)
   iax.set_ylim(-0.5,2)
   iax.set(xlabel='$i$')
   iax.set(ylabel='$[\ln A_n - \ln \Gamma(n/2+1)]/n$')
   iax.legend(frameon=False, labelspacing = 0.35)
   iax.text(10, 210, "inverse")


for i in (0,1):
   iax = ax[1][i]  
   iax.set_xlim(0,200)
   iax.set_ylim(0,250)
   iax.set(xlabel='$|i-d|$')
   iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')
   iax.legend(frameon=False, loc="lower right")

   
   ax[0][0].text(10, 1.8, "inverse")
   ax[0][1].text(160, 1.8, "direct")

   ax[1][0].text(10, 230, "inverse")
   ax[1][1].text(10, 230, "direct")

   
#iax.grid(axis="y")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')
#iax.set_yticks(np.arange(0, 2.5, step=0.5)) 





#-----------------------

plt.subplots_adjust(left=0.10, right=0.98, top=0.98, bottom=0.08, hspace=0.25, wspace=0.3)
#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.10, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
