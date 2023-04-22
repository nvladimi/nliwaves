
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 
  
    
outfile = "plotMv34.pdf"
alpha = 3/4


MS = 0.5
LW = 1
m = 60



#-- load data --

datM  = np.loadtxt("../POST34/q7_v34R28T_mts.txt")
M = datM[:,1:]


#alpha=3/4
oldfit = np.array([
   [0,   0,                  0],
   [1,   1.94595550094121,  0.584848661840946],
   [2,   4.12692626390787,  1.16757480235658],
   [3,   6.46492218106627,  1.74823849249817],
   [4,   8.92002195451477,  2.3269178882561],
   [5,   11.4683748858707,  2.90370517856197],
   [6,   14.0944870189613,  3.47870803404926],
   [7,   16.787897777618,   4.05205527099105],
   [8,   19.5415965198986,  4.62390609581986],
   [9,   22.3512535811863,  5.19446207352164],
   [10,  25.2148032752631,  5.76397702496796],
   [11,  28.1319390583751,  6.3327511685066],
   [12,  31.1030453386246,  6.90109015224401]
])

q=oldfit[:,0]
dzeta = oldfit[:,2]
delta = q*dzeta[3]/3 - dzeta
#delta = dzeta[3] - dzeta



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
   


iax.set_prop_cycle(None)

for n in range(12,0,-1):
   y = 0.76 + (x-10) * delta[n] * np.log(phi)/n
   iax.plot(x,y, '--', mfc='none', ms=MS,  lw=LW/2,  label="n = {}".format(n));


   
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
iax.set_ylim(0.7,1.15)
iax.set(xlabel='$i$', ylabel='$[\ln A_n - \ln \Gamma(n/2+1)]/n$')
#iax.grid(axis="y")
#iax.legend(frameon=False)
iax.text(10, 1.1, "$\\alpha = 3/4$")
iax.axvline(x=15, ymax=0.4, color='gray', lw = 0.5)
iax.axvline(x=30, ymax=0.4, color='gray', lw = 0.5)
#iax.set_yscale('log')


#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
