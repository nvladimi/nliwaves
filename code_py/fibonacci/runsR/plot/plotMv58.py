
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 
  
    
outfile = "plotMv58.pdf"
alpha = 5/8


MS = 0.5
LW = 1
m = 60



#-- load data --

datM  = np.loadtxt("../POST58/r1_v58R19_dt3T_mts.txt")
M = datM[:,1:]


#-- plot data --

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,3))

x = np.arange(1,m+1)

phi = (1 + np.sqrt(5))/2

c = phi**(1 + alpha)

iax = ax[0]
for n in range(12,0,-1):
   c = np.log(gamma(n/2+1))
   #c = np.log(factorial(n))
   y = np.log(datM[:,n+1]) - c
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));
   
iax = ax[1]
for n in range(12,0,-1):
   c = np.log(gamma(n/2+1))
   y = (np.log(datM[:,n+1]) - c)/n
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));
   

#-- canvas options --

iax = ax[0]  
iax.set_xlim(0,60)
iax.set_ylim(0,13)
iax.set(xlabel='$i$', ylabel='$\ln A_n - \ln \Gamma(n/2+1)$')
#iax.grid(axis="y")
iax.legend(frameon=False)
#iax.text(10, 0.5, "$\\alpha = 3/8$")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')

iax = ax[1]  
iax.set_xlim(0,60)
iax.set_ylim(0.96,1.10)
iax.set(xlabel='$i$', ylabel='$[\ln A_n - \ln \Gamma(n/2+1)]/n$')
#iax.grid(axis="y")
#iax.legend(frameon=False)
iax.text(10, 1.08, "$\\alpha = 5/8$")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')



#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
