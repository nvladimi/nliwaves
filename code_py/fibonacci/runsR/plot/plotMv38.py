
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 
    
outfile = "plotMv38.pdf"
alpha = 3/8


MS = 0.5
LW = 1
m = 60



#-- load data --

datM  = np.loadtxt("../POST38/r1_v38L06_dt3T_mts.txt")
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
#iax.set_ylim(0.1,5)
iax.set(xlabel='$i$', ylabel='$\ln A_n - \ln \Gamma(n/2+1)$')
#iax.grid(axis="y")
iax.legend(frameon=False)
#iax.text(10, 0.5, "$\\alpha = 3/8$")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')

iax = ax[1]  
iax.set_xlim(0,60)
iax.set_ylim(0.95,1.05)
iax.set(xlabel='$i$', ylabel='$[\ln A_n - \ln \Gamma(n/2+1)]/n$')
#iax.grid(axis="y")
#iax.legend(frameon=False)
iax.text(45, 1.04, "$\\alpha = 3/8$")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
