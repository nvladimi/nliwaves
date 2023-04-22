import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=12)


outfile = "plot_corr_alpha12_E33.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"

runs = ("s5m10", )

nmodes=200; i1=20; i2=180;

#c1="2211-00"; ishift=0;  combfactor=8
#c1="210-210"; ishift=0;  combfactor=8
c1="3221-10"; ishift=1;  combfactor=4
#c1="211-320"; ishift=1;  combfactor=4
#c1="4321-20"; ishift=2;  combfactor=2
#c1="221-430"; ishift=2;  combfactor=2
#c1="5421-30";  ishift=3;  combfactor=1
#c1="321-540";  ishift=3;  combfactor=1
#c1="6521-40";  ishift=4;  combfactor=1
#c1="421-650";  ishift=4;  combfactor=1
#c1="7621-50";  ishift=5;  combfactor=1
#c1="521-760";  ishift=5;  combfactor=1
#c1="8721-60";  ishift=6;  combfactor=1
#c1="621-870";  ishift=6;  combfactor=1
#c1="9821-70";  ishift=6;  combfactor=1
#c1="721-980";  ishift=6;  combfactor=1



c2="21-0";

#---------------------------------------------------------------------------------


#fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.0,4.5))
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4,3))

        
#-- common parameters --

MS = 1
LW = 1


#-- load and plot data --

for run in runs:

    if run in ('s5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
    
      
    iax = ax
 
    rname = "Post02/" + run  + "_" + c1 + "_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_" + c2 + "_mcr.txt"
    dat2 = np.loadtxt(rname)

    x  = dat1[i1:i2,0]
    sc = 1/dat1[i1:i2,2] * np.abs(x-d) / np.sqrt(combfactor)
    y1 = dat1[i1:i2,3] * sc
    y2 = dat2[i1:i2,3]*dat2[i1-ishift:i2-ishift,3] * sc
    y = y1-y2

    y1avg = np.average(y1)
    y2avg = np.average(y2)
    yavg  = np.average(y)
    yrms  = np.sqrt(np.average( (y - yavg)**2 ) )

  
    print(c1 + ":  ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  c1 + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'or',   mfc='none', ms=MS, lw=LW, label="y1-y2")
    iax.plot(x, y1, ':g',  mfc='none', ms=MS, lw=LW, label="y1")
    iax.plot(x, y2, ':b',  mfc='none', ms=MS, lw=LW, label="y2")
    iax.plot(x, x*0 + yavg, 'k-',  lw=LW/2)

    iax.legend(frameon=False)

    
#-- canvas options --

"""

for j in range(0,3):
    for i in range(0,2):
        ax[j][i].set_ylim(-ymax[j][i],ymax[j][i])
        ax[1][i].set(xlabel='$i$') # ylabel="$D_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
        ax[j][i].set_xlim(0,nmodes)
        ax[j][i].legend(frameon=False, handlelength=0.1, loc="upper left", bbox_to_anchor=(1.00,0.95))
        ax[j][i].text(1.05*nmodes, 0.9*ymax[j][i], leglbl[j][i], fontsize=12) 

"""

#-----------------------

plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.09, hspace=0.15, wspace=0.15)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
