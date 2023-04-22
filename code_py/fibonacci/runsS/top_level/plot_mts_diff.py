
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
from scipy.special import gamma, factorial 

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=13)



    
outfile = "plot_mts_diff.pdf"

direct=True

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4,4))

LW = 1.2
mycolors=['#d00000',  '#40a000', '#0040ff', '#e07000', '#00a0f0',  '#c000f0']


#---------------------------------

  
   
#---------------------------
iax=ax

#MS=0.5

#nrange = (1,2,3,6,9,12)
nrange = (1,4,8,12)

m=200

if direct:
   datM  = np.loadtxt("Post/s5m10_mts.txt"); d=m
else:
   datM  = np.loadtxt("Post/s5m90_mts.txt"); d=1


x = np.abs(d-np.arange(1,m+1))

ic=3
iax.set_prop_cycle(None)
for n in nrange:
   y = np.abs( datM[:,n+1] / gamma(n/2+1) / (datM[:,2+1])**(n/2) - 1 )
   iax.plot(x,y, mfc='none',  lw=LW,  label="{}".format(n), color=mycolors[ic]);
   ic = ic-1

iax.plot(x,1/x, '--k', mfc='none',  lw=LW/2,  label="$1/|i-d|$");
   
iax.set_xlim(1,200)
iax.set_ylim(1e-4,10)
iax.set(xlabel='$|i-d|$')
iax.set_yscale('log')
iax.set_xscale('log')

iax.set(ylabel='$| A_m  A_2^{-m/2} / \Gamma(\\frac{m}{2}+1) - 1 |$')
iax.yaxis.set_label_coords(-0.13,0.5)

iax.legend(frameon=False) # labelspacing = 0.32, handlelength=1,
           # loc="upper left", bbox_to_anchor=(0.02, 0.93))

#---------------------------

#-----------------------

plt.subplots_adjust(left=0.20, right=0.98, top=0.97, bottom=0.16, wspace=0.27)

plt.savefig(outfile) #pad_inches=0

plt.close()





#===========================================



#=========================================================
