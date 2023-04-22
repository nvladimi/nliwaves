
#import importlib ; import fibo_compensated_ksi
#importlib.reload(fibo_compensated_ksi) ; fibo_compensated_ksi.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo_compensated_ksi.pdf"


    #-- load data --

    Tv0   = np.loadtxt("../POST0/q7_v0dt1T21_mcr.txt")
    Tv1   = np.loadtxt("../POST1/q7b_g3500dt5T21_mcr.txt")
    Tv14  = np.loadtxt("../POST14/q7_v14L10T21_mcr.txt")
    Tv34  = np.loadtxt("../POST34/q7_v34R28T21_mcr.txt")

    Tv12L = np.loadtxt("../POSTcorr/r1_m50L18T21_mcr.txt")
    Tv12R = np.loadtxt("../POSTcorr/r1_m10R18T21_mcr.txt")
    Tv38  = np.loadtxt("../POSTcorr/r1_v38L06T21_mcr.txt")
    Tv58  = np.loadtxt("../POSTcorr/r1_v58R19T21_mcr.txt")
    
    
    
    #-- common parameters --
 

    fig, ax = plt.subplots(ncols=3, nrows=3, figsize=(9,9))
    
    MS = 0.5
    LW = 1
    m = 60
    
    x = np.arange(1,m+1)

    phi = (1 + np.sqrt(5))/2

    f0  = phi**(1 + 0)
    f1  = phi**(1 + 1)
    f12 = phi**(1 + 1/2)
    f14 = phi**(1 + 1/4)
    f34 = phi**(1 + 3/4)
    f38 = phi**(1 + 3/8)
    f58 = phi**(1 + 5/8)
    
     #-- plot data --

    #-- row 1: rescaled ksi

    iax = ax[0,0]
    y1 =  Tv0[:,4] /  Tv0[:,2] *f0;
    y2 = Tv14[:,4] / Tv14[:,2] *f14;
    y3 = Tv38[:,4] / Tv38[:,2] *f38;
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    
    iax = ax[0,1]
    y1 = Tv12R[:,4] / Tv12R[:,2] *f12;
    y2 = Tv12L[:,4] / Tv12L[:,2] *f12;
    iax.plot(60-x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="p = 10");
    iax.plot(x-1,-y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="p = 50");     
    iax.plot(x,0.87/x, 'k:',  mfc='none', ms=MS,  lw=LW);
 
    
    iax = ax[0,2]
    y1 =  Tv1[:,4] /  Tv1[:,2] *f1;
    y2 = Tv34[:,4] / Tv34[:,2] *f34;
    y3 = Tv58[:,4] / Tv58[:,2] *f58;   
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    #----- row 2: rescaled J

    iax = ax[1,0]
    y1 =  Tv0[:,4] *f0 
    y2 = Tv14[:,4] *f14 
    y3 = Tv38[:,4] *f38
    iax.plot(x,y1+1 , 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2+1, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3+1, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    
    iax = ax[1,1]
    y1 = Tv12R[:,4] *f12
    y2 = Tv12L[:,4] *f12
    iax.plot(60-x,y1-1, 'r-',  mfc='none', ms=MS,  lw=LW, label="$p = 10$");
    iax.plot(x-1,y2+1, 'b-',  mfc='none', ms=MS,  lw=LW, label="$p = 50$");     
 
    
    iax = ax[1,2]
    y1 =  Tv1[:,4] *f1
    y2 = Tv34[:,4] *f34
    y3 = Tv58[:,4] *f58
    iax.plot(x,y1-1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2-1, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3-1, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    #----- row 3: rescaled n_i = C_0

    iax = ax[2,0]
    y1 =  Tv0[:,1]
    y2 = Tv14[:,1]
    y3 = Tv38[:,1]
    iax.plot(x,y1 , 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    
    iax = ax[2,1]
    y1 = Tv12R[:,1]
    y2 = Tv12L[:,1]
    iax.plot(60-x,y1, 'r-',  mfc='none', ms=MS,  lw=LW, label="$p = 10$");
    iax.plot(x-1,y2, 'b-',  mfc='none', ms=MS,  lw=LW, label="$p = 50$");     
 
    
    iax = ax[2,2]
    y1 =  Tv1[:,1]
    y2 = Tv34[:,1]
    y3 = Tv58[:,1]
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    
    
    
    #-- canvas options --



    #--row 1

    iax = ax[0,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.3,0.05)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    iax.legend(frameon=False, loc='upper left', labelspacing=0.3)
    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0


    iax = ax[0,1]
    iax.set_xscale('log')
    iax.set_yscale('log')
    iax.set_xlim(1,100)
    iax.set_ylim(0.01,1)
    iax.set(xlabel='$|i - d|$', ylabel='$|\\xi|$')
    ticks = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])
    iax.legend(frameon=False)
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0
  

    iax = ax[0,2]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.05,0.3)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    iax.legend(frameon=False, loc='upper left', labelspacing=0.3)
    iax.yaxis.labelpad = -4
    iax.xaxis.labelpad = 0


    
    #--row 2


    iax = ax[1,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.02,0.02)
    iax.set(xlabel='$i$', ylabel='$J_i + 1$')
    iax.grid(axis="y")
    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0


    iax = ax[1,1]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.02,0.02)
    iax.set(xlabel='$|i - d|$', ylabel='$J_i \pm 1$')
    iax.grid(axis="y")
    iax.legend(frameon=False)
    iax.yaxis.labelpad = -4  
    iax.xaxis.labelpad = 0
  

    iax = ax[1,2]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.02,0.02)
    iax.set(xlabel='$i$', ylabel='$J_i - 1$')
    iax.grid(axis="y")
    iax.yaxis.labelpad = -4
    iax.xaxis.labelpad = 0


    
    #--row 3


    iax = ax[2,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(0,8)
    iax.set(xlabel='$i$', ylabel='$n=C_0$')
    iax.grid(axis="y")
    iax.yaxis.labelpad = 4
    iax.xaxis.labelpad = 0


    iax = ax[2,1]
    iax.set_xlim(0,60)
    iax.set_ylim(0,16)
    iax.set(xlabel='$|i - d|$', ylabel='$n=C_0$')
    iax.grid(axis="y")
    iax.legend(frameon=False)
    iax.yaxis.labelpad = 4  
    iax.xaxis.labelpad = 0
  

    iax = ax[2,2]
    iax.set_xlim(0,60)
    iax.set_ylim(0,8)
    iax.set(xlabel='$i$', ylabel='$n=C_0$')
    iax.grid(axis="y")
    iax.yaxis.labelpad = 4
    iax.xaxis.labelpad = 0





    
  
    #-----------------------

    plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.06, hspace=0.35, wspace=0.28)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 
   



#===========================================



#=========================================================
