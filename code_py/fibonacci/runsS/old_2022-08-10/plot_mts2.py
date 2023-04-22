
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 
  
    
outfile = "plot_mts2.pdf"

MS = 0.5
LW = 1
phi = (1 + np.sqrt(5))/2


#-- load data --


#alpha = 1/2; c = phi**(1 + alpha)



#-- plot data --

fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(10,6.5))


#---------------------------

iax = ax[0][0]
datM  = np.loadtxt("Post/s5m10_mts.txt")
M = datM[:,1:];  m = len(M)
x = m-np.arange(1,m+1)
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

datM  = np.loadtxt("Post/s4m10_mts.txt")
M = datM[:,1:];  m = len(M)
x = m-np.arange(1,m+1)
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW);


   
#---------------------------

iax = ax[1][0]
datM  = np.loadtxt("Post/s5m90_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-1
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

datM  = np.loadtxt("Post/s4m90_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-1
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW);

#---------------------------
#---------------------------

iax = ax[0][1]
datM  = np.loadtxt("Post/s3v58_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-10
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

datM  = np.loadtxt("Post/s4v58_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-10
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

#---------------------------

iax = ax[1][1]
datM  = np.loadtxt("Post/s3v38_mts.txt")
M = datM[:,1:];  m = len(M)
x = m-10-np.arange(1,m+1)
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

datM  = np.loadtxt("Post/s4v38_mts.txt")
M = datM[:,1:];  m = len(M)
iax.set_prop_cycle(None)
x = m-10-np.arange(1,m+1)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

#---------------------------
#---------------------------

iax = ax[0][2]
datM  = np.loadtxt("Post/s3v34_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-10
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

datM  = np.loadtxt("Post/s4v34_mts.txt")
M = datM[:,1:];  m = len(M)
x = np.arange(1,m+1)-10
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

#---------------------------

iax = ax[1][2]
datM  = np.loadtxt("Post/s3v14_mts.txt")
M = datM[:,1:];  m = len(M)
x = m-10-np.arange(1,m+1)
iax.set_prop_cycle(None)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

datM  = np.loadtxt("Post/s4v14_mts.txt")
M = datM[:,1:];  m = len(M)
iax.set_prop_cycle(None)
x = m-10-np.arange(1,m+1)
for n in range(12,0,-1):
   y = (datM[:,n+1] / gamma(n/2+1) )**(3/n)
   iax.plot(x,y, mfc='none', ms=MS,  lw=LW,  label="n = {}".format(n));

#---------------------------



   
#-- canvas options --

iax = ax[0][0]  
iax.set_xlim(0,200)
iax.set_ylim(0,250)
iax.set(xlabel='$|i-d|$')
iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')
#iax.grid(axis="y")
#iax.legend(frameon=False)
iax.text(10, 230, "$\\alpha = 1/2$")
iax.text(10, 210, "direct")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')
#iax.set_yticks(np.arange(0, 2.1, step=0.5)) 

iax = ax[1][0]  
iax.set_xlim(0,200)
iax.set_ylim(0,250)
iax.set(xlabel='$|i-d|$')
iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')
#iax.grid(axis="y")
iax.legend(frameon=False)
iax.text(10, 230, "$\\alpha = 1/2$")
iax.text(10, 210, "inverse")
#iax.axhline(y=7, color='gray', lw = 0.5)
#iax.set_yscale('log')
#iax.set_yticks(np.arange(0, 2.5, step=0.5)) 




iax = ax[0][1]  
iax.set_xlim(-10,90)
#iax.set_ylim(2.6,3.1)
iax.set_ylim(17,27)
iax.set(xlabel='$|i-p|$')
#iax.set_yticks(np.arange(0.95, 1.11, step=0.05)) 
iax.text(0, 26, "$\\alpha = 5/8$")
iax.text(0, 25.3, "direct")
iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')

iax = ax[1][1]  
iax.set_xlim(-10,90)
#iax.set_ylim(2.6,3.1)
iax.set_ylim(17,27)
iax.set(xlabel='$|i-p|$')
#iax.set_yticks(np.arange(0.95, 1.11, step=0.05)) 
iax.text(0, 26, "$\\alpha = 3/8$")
iax.text(0, 25.3, "inverse")
iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')



iax = ax[0][2]  
iax.set_xlim(-10,90)
#iax.set_ylim(1.9,3.0)
iax.set_ylim(6,24)
iax.set(xlabel='$|i-p|$')
#iax.set_yticks(np.arange(0.95, 1.11, step=0.05)) 
iax.text(0, 22, "$\\alpha = 3/4$")
iax.text(0, 20.5, "direct")
iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')

iax = ax[1][2]  
iax.set_xlim(-10,90)
#iax.set_ylim(1.9,3.0)
iax.set_ylim(6,24)
iax.set(xlabel='$|i-p|$')
#iax.set_yticks(np.arange(0.95, 1.11, step=0.05)) 
iax.text(0, 22, "$\\alpha = 1/4$")
iax.text(0, 20.5, "inverse")
iax.set(ylabel='$[A_n / \Gamma(n/2+1)]^{3/n}$')



#-----------------------

plt.subplots_adjust(left=0.06, right=0.99, top=0.98, bottom=0.08, hspace=0.25, wspace=0.23)
#plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.10, wspace=0.3)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
