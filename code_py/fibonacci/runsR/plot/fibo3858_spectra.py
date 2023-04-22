
#import importlib ; import fibo3858_spectra
#importlib.reload(fibo3858_spectra) ; fibo3858_spectra.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo3858_spectra.pdf"


    #-- load data --

    datv0   = np.loadtxt("../POSTFrisch/q7_v0dt1_flx.txt")
    datv1   = np.loadtxt("../POSTFrisch/q7b_g3500dt5_flx.txt")
    datv14  = np.loadtxt("../POSTFrisch/q7_v14L10_flx.txt")
    datv34  = np.loadtxt("../POSTFrisch/q7_v34R28_flx.txt")
    datv38  = np.loadtxt("../POST/r1_v38L06_dt3_flx.txt")
    datv58  = np.loadtxt("../POST/r1_v58R19_dt3_flx.txt")

 
 
    #data indexes: dat[:,ii]
   
    ii =  0
    iFi = 1
    ink = 2
    iPi = 3
    iJ =  4


    #-- plot data --

    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6,5.8))
    #inax1 = ax[0,0].inset_axes([0.55, 0.08, 0.40, 0.30])
    #inax2 = ax[0,1].inset_axes([0.10, 0.08, 0.4, 0.30])


    
    MS = 0.5
    LW = 1
    
    def f(j):
        return((-1)**j)

    Fi=datv0[:,iFi];
    Pi0 = 67.65
    
    iax = ax[0,0]
    x1=datv0[:,ii];   y1 = datv0[:,ink]*Fi;
    x2=datv14[:,ii];  y2 = datv14[:,ink]*Fi
    x3=datv38[:,ii];  y3 = datv38[:,ink]*Fi
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x2,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x3,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");
    iax.plot(x1, 40*Fi**(1/3), '--k', lw = LW/2)
    iax.plot(x1, 70*Fi**(1/6), '--k', lw = LW/2)
    iax.plot(x1, 120*Fi**(1/12), '--k', lw = LW/2)

    iax = ax[0,1]
    x1=datv1[:,ii];   y1 = datv1[:,ink]*Fi
    x2=datv34[:,ii];  y2 = datv34[:,ink]*Fi
    x3=datv58[:,ii];  y3 = datv58[:,ink]*Fi
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x2,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x3,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");
    iax.plot(x1, 90*Fi**(-1/3), '--k', lw = LW/2)
    iax.plot(x1, 100*Fi**(-1/6), '--k', lw = LW/2)
    iax.plot(x1, 130*Fi**(-1/12), '--k', lw = LW/2)
 
    iax = ax[1,0]
    x1=datv0[:,ii];   y1 = datv0[:,iJ] * (datv0[:,ink])**(-3/2);
    x2=datv14[:,ii];  y2 = datv14[:,iJ] * (datv14[:,ink])**(-3/2);
    x3=datv38[:,ii];  y3 = datv38[:,iJ] * (datv38[:,ink])**(-3/2);
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x2,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x3,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");
 
    iax = ax[1,1]
    x1=datv1[:,ii];   y1 = datv1[:,iJ] * (datv1[:,ink])**(-3/2);
    x2=datv34[:,ii];  y2 = datv34[:,iJ] * (datv34[:,ink])**(-3/2);
    x3=datv58[:,ii];  y3 = datv58[:,iJ] * (datv58[:,ink])**(-3/2);
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x2,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x3,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");
    
    """
    
    iax = inax1
    x1=datv0[:,ii];   y1 = datv0[:,iPi]/Pi0;
    x2=datv14[:,ii];  y2 = datv14[:,iPi]/Pi0
    x3=datv38[:,ii];  y3 = datv38[:,iPi]/Pi0
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x2,y2, 'bs',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x3,y3, 'gs',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");
 

    iax = inax2
    x1=datv1[:,ii];   y1 = datv1[:,iPi]/Pi0
    x2=datv34[:,ii];  y2 = datv34[:,iPi]/Pi0
    x3=datv58[:,ii];  y3 = datv58[:,iPi]/Pi0
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x2,y2, 'bs',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x3,y3, 'gs',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    """



    
    #-- canvas options --


    iax = ax[0,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(10,1e5)
    iax.set(xlabel='$i$', ylabel='$n_i F_i$')
    iax.set_yscale('log')
    iax.legend(frameon=False)
    iax.yaxis.labelpad = 3  
    iax.xaxis.labelpad = 0


    iax = ax[0,1]
    iax.set_xlim(0,60)
    iax.set_ylim(1e-2,100)
    iax.set(xlabel='$i$', ylabel='$n_i F_i$')
    iax.set_yscale('log')
    iax.legend(frameon=False)
    iax.yaxis.labelpad = -4 
    iax.xaxis.labelpad = 0



    iax = ax[1,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.3,0.05)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    iax.axhline(y=-0.22, color='gray', lw = 0.5)
    iax.axhline(y=-0.11, color='gray', lw = 0.5)
    iax.axhline(y=-0.055, color='gray', lw = 0.5)
    iax.axhline(y=0, color='gray', lw = 0.5)
    ticks = (-0.3, -0.2, -0.1, 0)
    iax.set_yticks(ticks)
    #iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))
    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0



    iax = ax[1,1]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.05,0.3)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    ticks = (0, 0.1, 0.2, 0.3)
    iax.set_yticks(ticks)
    iax.axhline(y=0, color='gray', lw = 0.5)
    iax.axhline(y=0.055, color='gray', lw = 0.5)
    iax.axhline(y=0.11, color='gray', lw = 0.5)
    iax.axhline(y=0.22, color='gray', lw = 0.5)
    iax.yaxis.labelpad = -4
    iax.xaxis.labelpad = 0


    """

    iax = inax1  
    iax.set_xlim(0,60)
    ticks=[0, 20, 40, 60]
    iax.set_xticks(ticks)
    iax.set_xticklabels(["$%g$"%y for y in ticks], fontsize=8, position=(0,0.05))
    iax.set_ylim(-1.25,0.25)
    ticks=[-1, 0]
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks], fontsize=8)
    iax.text(20, 0.35, '$\Pi_i / \Pi_M$')


    
    iax = inax2
    iax.set_xlim(0,60)
    ticks=[0, 20, 40, 60]
    iax.set_xticks(ticks)
    iax.set_xticklabels(["$%g$"%y for y in ticks], fontsize=8, position=(0,0.05))
    iax.set_ylim(-0.25,1.25)
    ticks=[0, 1]
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks], fontsize=8)
    iax.text(20, 1.35, '$\Pi_i / \Pi_M$')

    """

    #-----------------------

    plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.06, hspace=0.25, wspace=0.28)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 
   



#===========================================



#=========================================================
