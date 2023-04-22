
#import importlib ; import fibo_newksi
#importlib.reload(fibo_newksi) ; fibo_newksi.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo_newksi.pdf"


    #-- load data --

    Fv0   = np.loadtxt("../POST0/q7_v0dt1F21_mcr.txt")
    Fv1   = np.loadtxt("../POST1/q7b_g3500dt5F21_mcr.txt")
    Fv12L = np.loadtxt("../POST12L/r1_m50_L180F21_mcr.txt")
    Fv12R = np.loadtxt("../POST12R/r1_m10_R180F21_mcr.txt")
    Fv14  = np.loadtxt("../POST14/q7_v14L10F21_mcr.txt")
    Fv34  = np.loadtxt("../POST34/q7_v34R28F21_mcr.txt")
    Fv38  = np.loadtxt("../POST38/r1_v38L06_dt3F21_mcr.txt")
    Fv58  = np.loadtxt("../POST58/r1_v58R19_dt3F21_mcr.txt")

    Tv0   = np.loadtxt("../POST0/q7_v0dt1T21_mcr.txt")
    Tv1   = np.loadtxt("../POST1/q7b_g3500dt5T21_mcr.txt")
    Tv12L = np.loadtxt("../POST12L/r1_m50_L180T21_mcr.txt")
    Tv12R = np.loadtxt("../POST12R/r1_m10_R180T21_mcr.txt")
    Tv14  = np.loadtxt("../POST14/q7_v14L10T21_mcr.txt")
    Tv34  = np.loadtxt("../POST34/q7_v34R28T21_mcr.txt")
    Tv38  = np.loadtxt("../POST38/r1_v38L06_dt3T21_mcr.txt")
    Tv58  = np.loadtxt("../POST58/r1_v58R19_dt3T21_mcr.txt")


    
    
    #-- common parameters --
 

    fig, ax = plt.subplots(ncols=3, nrows=3, figsize=(9,9))
    #inax1 = ax[0,0].inset_axes([0.55, 0.08, 0.40, 0.30])
    #inax2 = ax[0,1].inset_axes([0.10, 0.08, 0.4, 0.30])
    
    MS = 0.5
    LW = 1
    m = 60
    
    x = np.arange(1,m+1)


     #-- plot data --

    #-- row 1: unscaled data

    iax = ax[0,0]
    y1 =  Fv0[:,4] /  Fv0[:,1]**(3/2);
    y2 = Fv14[:,4] / Fv14[:,1]**(3/2);
    y3 = Fv38[:,4] / Fv38[:,1]**(3/2);
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    
    iax = ax[0,1]
    y1a = Fv12R[:,4] / Fv12R[:,2];
    y1b = Fv12R[:,4] / Fv12R[:,1]**(3/2);
    y2a = Fv12L[:,4] / Fv12L[:,2];
    y2b = Fv12L[:,4] / Fv12L[:,1]**(3/2);

    iax.plot(60-x, y1a, 'r-',  mfc='none', ms=MS,  lw=LW, label="new ksi");
    iax.plot(60-x, y1b, 'r--', mfc='none', ms=MS,  lw=LW, label="old ksi");
    iax.plot(x-1, -y2a, 'b-',  mfc='none', ms=MS,  lw=LW, label="new ksi") ;     
    iax.plot(x-1, -y2b, 'b--', mfc='none', ms=MS,  lw=LW, label="old ksi");     
    iax.plot(x,0.87/x, 'k:',  mfc='none', ms=MS,  lw=LW);
 
    
    iax = ax[0,2]
    y1 =  Fv1[:,4] /  Fv1[:,1]**(3/2);
    y2 = Fv34[:,4] / Fv34[:,1]**(3/2);
    y3 = Fv58[:,4] / Fv58[:,1]**(3/2);   
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    #-- row 2: rescaled data

    iax = ax[1,0]
    y1 =  Tv0[:,4] /  Tv0[:,2];
    y2 = Tv14[:,4] / Tv14[:,2];
    y3 = Tv38[:,4] / Tv38[:,2];
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    
    iax = ax[1,1]
    y1a = Tv12R[:,4] / Tv12R[:,2];
    y1b = Tv12R[:,4] / Tv12R[:,1]**(3/2);
    y2a = Tv12L[:,4] / Tv12L[:,2];
    y2b = Tv12L[:,4] / Tv12L[:,1]**(3/2);
    iax.plot(60-x,y1a, 'r-',  mfc='none', ms=MS,  lw=LW,  label="new ksi");
    iax.plot(60-x,y1b, 'r--', mfc='none', ms=MS,  lw=LW,  label="old ksi");
    iax.plot(x-1,-y2a, 'b-',  mfc='none', ms=MS,  lw=LW,  label="new ksi");     
    iax.plot(x-1,-y2b, 'b--', mfc='none', ms=MS,  lw=LW,  label="old ksi");     
    iax.plot(x,0.87/x, 'k:',  mfc='none', ms=MS,  lw=LW);
 
    
    iax = ax[1,2]
    y1 =  Tv1[:,4] /  Tv1[:,2];
    y2 = Tv34[:,4] / Tv34[:,2];
    y3 = Tv58[:,4] / Tv58[:,2];   
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    #----- row 3

    iax = ax[2,0]
    y1 =  Tv0[:,4]
    y2 = Tv14[:,4]
    y3 = Tv38[:,4]
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    
    iax = ax[2,1]
    y1 = Tv12R[:,4]
    y2 = Tv12L[:,4]
    iax.plot(60-x,y1, 'r-',  mfc='none', ms=MS,  lw=LW, label="$p = 10$");
    iax.plot(x-1,-y2, 'b-',  mfc='none', ms=MS,  lw=LW, label="$p = 50$");     
 
    
    iax = ax[2,2]
    y1 =  Tv1[:,4] 
    y2 = Tv34[:,4] 
    y3 = Tv58[:,4] 
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    
    
    
    #-- canvas options --

    #-- row 1:

    iax = ax[0,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.3,0.05)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    ticks = (-0.22, -0.11, -0.055, 0)
    iax.set_yticks(ticks)
    iax.grid(axis="y")
    #iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))
    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0
    iax.text(20,0.02, "unscaled data, old ksi")
    

    iax = ax[0,1]
    iax.set_xscale('log')
    iax.set_yscale('log')
    iax.set_xlim(1,100)
    iax.set_ylim(0.01,1)
    #iax.set_ylim(0,0.04)
    iax.set(xlabel='$|i - d|$', ylabel='$|\\xi|$')
    ticks = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
    iax.legend(frameon=False)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])
    iax.text(2,0.7, "unscaled data") 
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0
   
    

    iax = ax[0,2]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.05,0.25)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    ticks = (0, 0.053, 0.105, 0.180)
    iax.set_yticks(ticks)
    iax.grid(axis="y")
    iax.text(20,0.23, "unscaled data, old ksi")
 
    iax.yaxis.labelpad = -4
    iax.xaxis.labelpad = 0



    #--row 2

    iax = ax[1,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.18,0.03)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    ticks = (-0.145, -0.06, -0.028, 0)
    iax.set_yticks(ticks)
    iax.grid(axis="y")
    iax.text(5,0.01, "scaled or unscaled data, new ksi")
    #iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))
    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0


    iax = ax[1,1]
    iax.set_xscale('log')
    iax.set_yscale('log')
    iax.set_xlim(1,100)
    iax.set_ylim(0.01,1)
    #iax.set_ylim(0,0.04)
    iax.set(xlabel='$|i - d|$', ylabel='$|\\xi|$')
    ticks = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])
    iax.text(2,0.7, "scaled data")
    iax.legend(frameon=False)
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0
  

    iax = ax[1,2]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.02,0.095)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    iax.text(5,0.085, "scaled or unscaled data, new ksi")
    ticks = (0, 0.024, 0.045, 0.070)
    iax.set_yticks(ticks)
    iax.grid(axis="y")
    iax.yaxis.labelpad = -4
    iax.xaxis.labelpad = 0


    
    #--row 3

    iax = ax[2,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.65,-0.45)
    iax.set(xlabel='$i$', ylabel='$J_i$')
    #iax.legend(frameon=False)
    ticks = (-0.65, -0.62, -0.55, -0.51, -0.45)
    iax.set_yticks(ticks)
    iax.grid(axis="y")
    iax.text(10,-0.48, "scaled data")
    #iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))
    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0


    iax = ax[2,1]
    iax.set_xlim(0,60)
    iax.set_ylim(0.44,0.52)
    #iax.set_ylim(0,0.04)
    iax.set(xlabel='$|i - d|$', ylabel='$|J_i|$')
    ticks = (0.44, 0.46, 0.48, 0.50, 0.52)
    iax.set_yticks(ticks)
    #iax.set_yticklabels(["$%g$"%y for y in ticks])
    iax.grid(axis="y")
    iax.text(10,0.51, "scaled data")
    iax.legend(frameon=False)
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0
  

    iax = ax[2,2]
    iax.set_xlim(0,60)
    iax.set_ylim(0.34,0.5)
    iax.set(xlabel='$i$', ylabel='$J_i$')
    #iax.legend(frameon=False)
    iax.text(10,0.48, "scaled data")
    ticks = (0.34, 0.38, 0.43, 0.46, 0.5)
    iax.set_yticks(ticks)
    iax.grid(axis="y")
    iax.yaxis.labelpad = -4
    iax.xaxis.labelpad = 0







    
  
    #-----------------------

    plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.06, hspace=0.25, wspace=0.28)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 
   



#===========================================



#=========================================================
