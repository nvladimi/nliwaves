
#import importlib ; import fibo3858_flux
#importlib.reload(fibo3858_flux) ; fibo3858_flux.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo3858_flux.pdf"


    #-- load data --

    datv0   = np.loadtxt("../POSTFrisch/q7_v0dt1_flx.txt")
    datv1   = np.loadtxt("../POSTFrisch/q7b_g3500dt5_flx.txt")
    datv14  = np.loadtxt("../POSTFrisch/q7_v14L10_flx.txt")
    datv34  = np.loadtxt("../POSTFrisch/q7_v34R28_flx.txt")
    datv38  = np.loadtxt("../POST/r1_v38L06_dt3_flx.txt")
    datv58  = np.loadtxt("../POST/r1_v58R19_dt3_flx.txt")
    
    datv12R  = np.loadtxt("../POST/r1_m10_R180_flx.txt")
    datv12L  = np.loadtxt("../POST/r1_m50_L180_flx.txt")


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
    phi = (1 + np.sqrt(5))/2

    def J(alpha):
        j = np.arange(60) + 1
        return  Pi0/2 * 5**(alpha/2) * phi**(- j * (1 + alpha) + 1 + 2*alpha )
    
    iax = ax[0,0]
    x1=datv0[:,ii];   y1 = datv0[:,iJ] / J(0)
    x2=datv14[:,ii];  y2 = datv14[:,iJ] / J(1/4)
    x3=datv38[:,ii];  y3 = datv38[:,iJ] /J(3/8)
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 0$");
    iax.plot(x2,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/4$");
    iax.plot(x3,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8$");

    iax = ax[0,1]
    x1=datv1[:,ii];   y1 = datv1[:,iJ] /J(1)
    x2=datv34[:,ii];  y2 = datv34[:,iJ] /J(3/4)
    x3=datv58[:,ii];  y3 = datv58[:,iJ] /J(5/8)
    iax.plot(x1,y1, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1$");
    iax.plot(x2,y2, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/4$");
    iax.plot(x3,y3, 'g-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8$");


    iax = ax[1,0]
    x4=datv12L[:,ii];   y4 = datv12L[:,iJ] / J(1/2)
    iax.plot(x4,y4, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/2$");

    iax = ax[1,1]
    x5=datv12R[:,ii];   y5 = datv12R[:,iJ] / J(1/2)
    iax.plot(x5,y5, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 1/2$");





    
    #-- canvas options --


    iax = ax[0,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-1.2,0.2)
    iax.set(xlabel='$i$', ylabel='LHS/RHS')
    iax.legend(frameon=False)
    iax.axhline(y=-1, color='gray', lw = 0.5)
    iax.axhline(y=0, color='gray', lw = 0.5)
    iax.yaxis.labelpad = 0  
    iax.xaxis.labelpad = 0


    iax = ax[0,1]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.2,1.2)
    iax.set(xlabel='$i$', ylabel='LHS/RHS')
    iax.legend(frameon=False)
    iax.axhline(y=1, color='gray', lw = 0.5)
    iax.axhline(y=0, color='gray', lw = 0.5)
    iax.yaxis.labelpad = 0 
    iax.xaxis.labelpad = 0




    iax = ax[1,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-1.2,0.2)
    iax.set(xlabel='$i$', ylabel='LHS/RHS')
    iax.legend(frameon=False)
    iax.axhline(y=-1, color='gray', lw = 0.5)
    iax.axhline(y=0, color='gray', lw = 0.5)
    iax.yaxis.labelpad = 0  
    iax.xaxis.labelpad = 0


    iax = ax[1,1]
    iax.set_xlim(0,60)
    iax.set_ylim(-0.2,1.2)
    iax.set(xlabel='$i$', ylabel='LHS/RHS')
    iax.legend(frameon=False)
    iax.axhline(y=1, color='gray', lw = 0.5)
    iax.axhline(y=0, color='gray', lw = 0.5)
    iax.yaxis.labelpad = 0
    iax.xaxis.labelpad = 0







    

    #-----------------------

    plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.06, hspace=0.25, wspace=0.28)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 
   



#===========================================



#=========================================================
