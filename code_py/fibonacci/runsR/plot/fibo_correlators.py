
#import importlib ; import fibo_correlators
#importlib.reload(fibo_correlators) ; fibo_correlators.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo_correlators.pdf"


    #-- load data --

    corr="21"
    ylabel = '$C_{' + corr + '}$' 
    
    Tv0   = np.loadtxt("../POST0/q7_v0dt1T"       + corr + "_mcr.txt")
    Tv1   = np.loadtxt("../POST1/q7b_g3500dt5T"   + corr + "_mcr.txt")
    Tv12L = np.loadtxt("../POST12L/r1_m50_L180T"  + corr + "_mcr.txt")
    Tv12R = np.loadtxt("../POST12R/r1_m10_R180T"  + corr + "_mcr.txt")
    Tv14  = np.loadtxt("../POST14/q7_v14L10T"     + corr + "_mcr.txt")
    Tv34  = np.loadtxt("../POST34/q7_v34R28T"     + corr + "_mcr.txt")
    Tv38  = np.loadtxt("../POST38/r1_v38L06_dt3T" + corr + "_mcr.txt")
    Tv58  = np.loadtxt("../POST58/r1_v58R19_dt3T" + corr + "_mcr.txt")


    
    
    #-- common parameters --
 

    fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(9,3))
    
    MS = 0.5
    LW = 1
    m = 60
    
    x = np.arange(1,m+1)


     #-- plot data --


    iax = ax[0]
    y1 =  Tv0[:,4]
    y2 = Tv14[:,4]
    y3 = Tv38[:,4]
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");
    y1 =  Tv0[:,3]
    y2 = Tv14[:,3]
    y3 = Tv38[:,3]
    iax.plot(x,y1, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");


    
    
    iax = ax[1]
    y1 = Tv12R[:,4]
    y2 = Tv12L[:,4]
    iax.plot(60-x,y1, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x-1,y2, 'b-',  mfc='none', ms=MS,  lw=LW);     
    y1 = Tv12R[:,3]
    y2 = Tv12L[:,3]
    iax.plot(60-x,y1, 'r:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x-1,y2, 'b:',  mfc='none', ms=MS,  lw=LW);     

    
    
    iax = ax[2]
    y1 =  Tv1[:,4] 
    y2 = Tv34[:,4] 
    y3 = Tv58[:,4] 
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");
    y1 =  Tv1[:,3] 
    y2 = Tv34[:,3] 
    y3 = Tv58[:,3] 
    iax.plot(x,y1, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    
    
    
    #-- canvas options --

 
    iax = ax[0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-1,1)
    iax.set(xlabel='$i$', ylabel=ylabel)
    #iax.legend(frameon=False)
    #ticks = (-0.65, -0.62, -0.55, -0.51, -0.45)
    #iax.set_yticks(ticks)
    #iax.grid(axis="y")
    #iax.text(10,-0.48, "scaled data")
    #iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))
    iax.yaxis.labelpad = -2
    iax.xaxis.labelpad = 2


    iax = ax[1]
    iax.set_xlim(0,60)
    iax.set_ylim(-1,1)
    #iax.set_ylim(0,0.04)
    iax.set(xlabel='$|i - d|$', ylabel=ylabel)
    #ticks = (0.44, 0.46, 0.48, 0.50, 0.52)
    #iax.set_yticks(ticks)
    #iax.set_yticklabels(["$%g$"%y for y in ticks])
    #iax.grid(axis="y")
    #iax.text(10,0.51, "scaled data")
    #iax.legend(frameon=False)
    iax.yaxis.labelpad = -2  
    iax.xaxis.labelpad = 2
  

    iax = ax[2]
    iax.set_xlim(0,60)
    iax.set_ylim(-1,1)
    iax.set(xlabel='$i$', ylabel=ylabel)
    #iax.legend(frameon=False)
    #iax.text(10,0.48, "scaled data")
    #ticks = (0.34, 0.38, 0.43, 0.46, 0.5)
    #iax.set_yticks(ticks)
    #iax.grid(axis="y")
    iax.yaxis.labelpad = -2
    iax.xaxis.labelpad = 2

 



    
  
    #-----------------------

    plt.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.15, hspace=0.25, wspace=0.28)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 
   



#===========================================



#=========================================================
