
#import importlib ; import fibo_compensated_corr
#importlib.reload(fibo_compensated_corr) ; fibo_compensated_corr.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo_compensated_corr.pdf"


    #-- load data --

    #corr="21";   n=(2+1)/3;  ylim=(-1.5,1.5)
    #corr="322";  n=(3+2+2)/3; ylim=(-8,8)
    #corr="431";  n=(4+3+1)/3; ylim=(-8,8)
    corr="6531"; n=(6+5+3+1)/3; ylim=(-20,20)
    
    ylabel = '$C_{' + corr + '}$' 
    
    Tv0   = np.loadtxt("../POST0/q7_v0dt1T"       + corr + "_mcr.txt")
    Tv1   = np.loadtxt("../POST1/q7b_g3500dt5T"   + corr + "_mcr.txt")
    Tv14  = np.loadtxt("../POST14/q7_v14L10T"     + corr + "_mcr.txt")
    Tv34  = np.loadtxt("../POST34/q7_v34R28T"     + corr + "_mcr.txt")

    Tv12L = np.loadtxt("../POSTcorr/r1_m50L18T"  + corr + "_mcr.txt")
    Tv12R = np.loadtxt("../POSTcorr/r1_m10R18T"  + corr + "_mcr.txt")
    Tv38  = np.loadtxt("../POSTcorr/r1_v38L06T"  + corr + "_mcr.txt")
    Tv58  = np.loadtxt("../POSTcorr/r1_v58R19T"  + corr + "_mcr.txt")


       
    #-- common parameters --
 

    fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(9,3))
    
    MS = 0.5
    LW = 1
    m = 60
    
    x = np.arange(1,m+1)


    phi = (1 + np.sqrt(5))/2

    f0  = phi**((1 + 0)*n)
    f1  = phi**((1 + 1)*n)
    f12 = phi**((1 + 1/2)*n)
    f14 = phi**((1 + 1/4)*n)
    f34 = phi**((1 + 3/4)*n)
    f38 = phi**((1 + 3/8)*n)
    f58 = phi**((1 + 5/8)*n)


    
     #-- plot data --


    iax = ax[0]
    y1 =  Tv0[:,4]*f0
    y2 = Tv14[:,4]*f14
    y3 = Tv38[:,4]*f38
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");
    y1 =  Tv0[:,3]*f0
    y2 = Tv14[:,3]*f14
    y3 = Tv38[:,3]*f38
    iax.plot(x,y1, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x,y2, 'b:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x,y3, 'g:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");


    
    
    iax = ax[1]
    y1 = Tv12R[:,4]*f12
    y2 = Tv12L[:,4]*f12
    iax.plot(60-x,y1, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x-1,y2, 'b-',  mfc='none', ms=MS,  lw=LW);     
    y1 = Tv12R[:,3]*f12
    y2 = Tv12L[:,3]*f12
    iax.plot(60-x,y1, 'r:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x-1,y2, 'b:',  mfc='none', ms=MS,  lw=LW);     

    
    
    iax = ax[2]
    y1 =  Tv1[:,4]*f1 
    y2 = Tv34[:,4]*f34 
    y3 = Tv58[:,4]*f58 
    iax.plot(x,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");
    y1 =  Tv1[:,3]*f1 
    y2 = Tv34[:,3]*f34 
    y3 = Tv58[:,3]*f58 
    iax.plot(x,y1, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x,y2, 'b:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x,y3, 'g:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    
    
    
    #-- canvas options --

    
    iax = ax[0]  
    iax.set_xlim(0,60)
    iax.set_ylim(ylim)
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
    iax.set_ylim(ylim)
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
    iax.set_ylim(ylim)
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
