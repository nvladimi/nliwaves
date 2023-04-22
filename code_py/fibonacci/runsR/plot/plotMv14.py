
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 
  
    
outfile = "plotMv14.pdf"
alpha = 1/4


MS = 0.5
LW = 1
m = 60



#-- load data --

datM  = np.loadtxt("../POST14/q7_v14L10T_mts.txt")
M = datM[:,1:]


oldfit = np.array([
   [0,   0,                 0],
   [1,   1.7178958450344,   0.414681602162036],
   [2,   3.73527871383263,  0.831762089133451],
   [3,   5.97233759711943,  1.25113641083553],
   [4,   8.38697948606695,  1.67271674803373],
   [5,  10.9529360301602,   2.09642649203429],
   [6,  13.652039760714,    2.52219550638833],
   [7,  16.470826667468,    2.94995774557976],
   [8,  19.398841272157,    3.37965210815355],
   [9,  22.4277250651091,   3.81122627793124],
   [10, 25.5506696168943,   4.24464090325529],
   [11, 28.7619743282338,   4.67986854465576],
   [12, 32.0565277706819,   5.11688098293115]
])

q=oldfit[:,0]
dzeta = oldfit[:,2]
delta = q*dzeta[3]/3 - dzeta
#delta = dzeta[3] - dzeta

print(delta)



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
   y = 0.745 + (x-50) * delta[n] * np.log(phi)/n
   iax.plot(x,y, '--', mfc='none', ms=MS,  lw=LW/2,  label="n = {}".format(n));



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
iax.set_ylim(0.6,1.05)
iax.set(xlabel='$i$', ylabel='$[\ln A_n - \ln \Gamma(n/2+1)]/n$')
#iax.grid(axis="y")
#iax.legend(frameon=False)
iax.text(45, 1.00, "$\\alpha = 1/4$")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')
iax.axvline(x=45, ymin=0, ymax=0.7, color='gray', lw = 0.5)
iax.axvline(x=30, ymin=0, ymax=0.7, color='gray', lw = 0.5)



#-----------------------

#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)
plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.15, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
