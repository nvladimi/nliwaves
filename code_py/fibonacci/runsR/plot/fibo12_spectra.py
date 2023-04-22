
#import importlib ; import fibo12_spectra
#importlib.reload(fibo12_spectra) ; fibo12_spectra.all()


import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def all(): 

    outfile = "fibo12_spectra.pdf"


    #-- load data --

    datv38  = np.loadtxt("../POST/r1_v38L06_dt3_flx.txt")
    datv58  = np.loadtxt("../POST/r1_v58R19_dt3_flx.txt")

    dat12Ra  = np.loadtxt("../POST/r1_m10_R180_flx.txt")   #R180  170
    dat12Rb  = np.loadtxt("../POST/r1_m10_R050_flx.txt")   #R050
    dat12Rc  = np.loadtxt("../POST/r1_m10_R400_flx.txt")   #R800

    dat12La  = np.loadtxt("../POST/r1_m50_L180_flx.txt")   #L180  170
    dat12Lb  = np.loadtxt("../POST/r1_m50_L050_flx.txt")   #L050
    dat12Lc  = np.loadtxt("../POST/r1_m50_L400_flx.txt")   #L800


 
    #data indexes: dat[:,ii]
   
    ii =  0
    iFi = 1
    ink = 2
    iPi = 3
    iJ =  4


    #-- plot data --

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(10,6))
    #inax1 = ax[0,0].inset_axes([0.55, 0.08, 0.40, 0.30])
    #inax2 = ax[0,1].inset_axes([0.10, 0.08, 0.4, 0.30])


    
    MS = 0.5
    LW = 1
    
    def f(j):
        return((-1)**j)

    Fi=dat12Ra[:,iFi];
    Pi0 = 67.65
    
    iax = ax[0,0]
    x1a=dat12Ra[:,ii];  y1a = dat12Ra[:,ink]*Fi
    x1b=dat12Rb[:,ii];  y1b = dat12Rb[:,ink]*Fi
    x1c=dat12Rc[:,ii];  y1c = dat12Rc[:,ink]*Fi
    x2a=dat12La[:,ii];  y2a = dat12La[:,ink]*Fi
    x2b=dat12Lb[:,ii];  y2b = dat12Lb[:,ink]*Fi
    x2c=dat12Lc[:,ii];  y2c = dat12Lc[:,ink]*Fi
    x3=datv38[:,ii];  y3 = datv38[:,ink]*Fi
    x4=datv58[:,ii];  y4 = datv58[:,ink]*Fi
    iax.plot(x1a,y1a, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$p=10,  \\gamma_R = 1.8$");
    iax.plot(x1b,y1b, 'r--', mfc='none', ms=MS,  lw=LW,  label="$p=10,  \\gamma_R = 0.5$");
    iax.plot(x1c,y1c, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$p=10,  \\gamma_R = 4.0$");
    iax.plot(x2a,y2a, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$p=50,  \\gamma_L = 1.8$");
    iax.plot(x2b,y2b, 'b--', mfc='none', ms=MS,  lw=LW,  label="$p=50,  \\gamma_L = 0.5$");
    iax.plot(x2c,y2c, 'b:',  mfc='none', ms=MS,  lw=LW,  label="$p=50,  \\gamma_L = 4.0$");
    iax.plot(x4,y4, 'm:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8, \\gamma_R = 2.7$");
    iax.plot(x3,y3, 'g:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8, \\gamma_L = 0.6$");

    
    iax = ax[0,1]
    iax.plot(x1a,y1a, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x1b,y1b, 'r--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x1c,y1c, 'r:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x2a,y2a, 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x2b,y2b, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x2c,y2c, 'b:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x1b,2.8*x1b+40, 'k--',  mfc='none', ms=MS,  lw=0.5*LW);
    iax.plot(x1b,2.8*(60-x1b)+45, 'k--',  mfc='none', ms=MS,  lw=0.5*LW);


    iax = ax[0,2]
    x1a=dat12Ra[:,ii];  y1a = dat12Ra[:,iPi]/Pi0;
    x1b=dat12Rb[:,ii];  y1b = dat12Rb[:,iPi]/Pi0;
    x1c=dat12Rc[:,ii];  y1c = dat12Rc[:,iPi]/Pi0;
    x2a=dat12La[:,ii];  y2a = dat12La[:,iPi]/Pi0
    x2b=dat12Lb[:,ii];  y2b = dat12Lb[:,iPi]/Pi0
    x2c=dat12Lc[:,ii];  y2c = dat12Lc[:,iPi]/Pi0
    x3=datv38[:,ii];  y3 = datv38[:,iPi]/Pi0
    x4=datv58[:,ii];  y4 = datv58[:,iPi]/Pi0
    iax.plot(x1a,y1a, 'r-',  mfc='none', ms=MS,  lw=LW,  label="$p=10,  \\gamma_R = 1.8$");
    iax.plot(x1b,y1b, 'r--', mfc='none', ms=MS,  lw=LW,  label="$p=10,  \\gamma_R = 0.5$");
    iax.plot(x1c,y1c, 'r:',  mfc='none', ms=MS,  lw=LW,  label="$p=10,  \\gamma_R = 4.0$");
    iax.plot(x2a,y2a, 'b-',  mfc='none', ms=MS,  lw=LW,  label="$p=50,  \\gamma_L = 1.8$");
    iax.plot(x2b,y2b, 'b--', mfc='none', ms=MS,  lw=LW,  label="$p=50,  \\gamma_L = 0.5$");
    iax.plot(x2c,y2c, 'b:',  mfc='none', ms=MS,  lw=LW,  label="$p=50,  \\gamma_L = 4.0$");
    iax.plot(x4,y4, 'm:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 5/8, \\gamma_R = 2.7$");
    iax.plot(x3,y3, 'g:',  mfc='none', ms=MS,  lw=LW,  label="$\\alpha = 3/8, \\gamma_L = 0.6$");



    
    
    iax = ax[1,0]
    x1a=dat12Ra[:,ii];  y1a = dat12Ra[:,iJ] * (dat12Ra[:,ink])**(-3/2);
    x1b=dat12Rb[:,ii];  y1b = dat12Rb[:,iJ] * (dat12Rb[:,ink])**(-3/2);
    x1c=dat12Rc[:,ii];  y1c = dat12Rc[:,iJ] * (dat12Rc[:,ink])**(-3/2);
    x2a=dat12La[:,ii];  y2a = dat12La[:,iJ] * (dat12La[:,ink])**(-3/2);
    x2b=dat12Lb[:,ii];  y2b = dat12Lb[:,iJ] * (dat12Lb[:,ink])**(-3/2);
    x2c=dat12Lc[:,ii];  y2c = dat12Lc[:,iJ] * (dat12Lc[:,ink])**(-3/2);
    x3=datv38[:,ii];  y3 = datv38[:,iJ] * (datv38[:,ink])**(-3/2);
    x4=datv58[:,ii];  y4 = datv58[:,iJ] * (datv58[:,ink])**(-3/2);
    iax.plot(x1a,y1a, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x1b,y1b, 'r--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x1c,y1c, 'r:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x2a,y2a, 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x2b,y2b, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x2c,y2c, 'b:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x4,y4, 'm:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x3,y3, 'g:',  mfc='none', ms=MS,  lw=LW);

    iax = ax[1,1]
    iax.plot(60-x1a,y1a, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(60-x1b,y1b, 'r--', mfc='none', ms=MS,  lw=LW);
    iax.plot(60-x1c,y1c, 'r:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x2a-1,-y2a, 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x2b-1,-y2b, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x2c-1,-y2c, 'b:',  mfc='none', ms=MS,  lw=LW);
    
    iax.plot(x1a,0.87/x1a, 'k:',  mfc='none', ms=MS,  lw=LW);
 
    





    
    #-- canvas options --


    iax = ax[0,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(10,1e3)
    iax.set(xlabel='$i$', ylabel='$n_i F_i$')
    iax.set_yscale('log')
    #iax.legend(frameon=False) # location="SouthEast")
    iax.yaxis.labelpad = 2 
    iax.xaxis.labelpad = 0


 
    iax = ax[0,1]  
    iax.set_xlim(0,60)
    iax.set_ylim(0,250)
    iax.set(xlabel='$i$', ylabel='$n_i F_i$')
    #iax.set_yscale('log')
    #iax.legend(frameon=False) # location="SouthEast")
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0
   

    iax = ax[0,2]  
    iax.set_xlim(0,60)
    #ticks=[0, 20, 40, 60]
    #iax.set_xticks(ticks)
    #iax.set_xticklabels(["$%g$"%y for y in ticks]) #, fontsize=8, position=(0,0.05))
    iax.set_ylim(-1.25,1.25)
    ticks=[-1, 0, 1]
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])#, fontsize=8)
    iax.set(xlabel = '$i$', ylabel = '$\Pi_i / \Pi_M$')
    iax.legend(frameon=False)
    iax.axhline(y=1, color='gray', lw = 0.5)
    iax.axhline(y=-1, color='gray', lw = 0.5)
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0





    iax = ax[1,0]  
    iax.set_xlim(0,60)
    iax.set_ylim(-0.15,0.15)
    iax.set(xlabel='$i$', ylabel='$\\xi$')
    #iax.legend(frameon=False)
    iax.axhline(y=0, color='gray', lw = 0.5)
    #ticks = (-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06)
    #iax.set_yticks(ticks)
    #iax.set_yticklabels(["$%g$"%y for y in ticks])
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
    iax.yaxis.labelpad = 2  
    iax.xaxis.labelpad = 0



    iax = ax[1,2]
    iax.text(0.4, 0.5, 'unused')
    






    iax.yaxis.labelpad = 2
    iax.xaxis.labelpad = 0


  
 
    

    #-----------------------

    plt.subplots_adjust(left=0.09, right=0.98, top=0.98, bottom=0.06, hspace=0.25, wspace=0.28)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 
   



#===========================================



#=========================================================
