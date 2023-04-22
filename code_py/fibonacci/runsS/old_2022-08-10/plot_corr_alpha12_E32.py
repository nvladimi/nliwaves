import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=12)


outfile = "plot_corr_alpha12_E.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"

runs = ("s5m90", )

#runs = ("s5m90", "s5m10")
#runs = ("s4v14", "s4v34")
#runs = ("s4v38", "s4v58")


#ymax = ((0.6, 0.06),  (0.15, 0.08),  (0.05, 0.05) );
nmodes=200





leglbl = (("$m = 3,4$",  "$m=5$"),  ("$m=5$", "$m=6$"),  ("$m=6$", "$m=6$") )




#---------------------------------------------------------------------------------


#fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.0,4.5))
fig, ax = plt.subplots(ncols=2, nrows=6, figsize=(10,11))

        
#-- common parameters --

MS = 1
LW = 1

colors = ['#1f77b4',
          '#ff7f0e',
          '#2ca02c',
          '#d62728',
          '#9467bd',
          '#8c564b',
          '#e377c2',
          '#7f7f7f',
          '#bcbd22',
          '#17becf',
          '#1a55FF']



#-- load and plot data --

for run in runs:

    if run in ('s5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
    
    #-- E23 ~ C210-00 --
    
    iax = ax[0][0]
    icolor = 0;

    rname = "Post02/" + run  + "_210-00_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_21-0_mcr.txt"
    dat2 = np.loadtxt(rname)

    x  = dat1[:,0]
    sc = 1/dat1[:,2] * np.abs(x-d) / np.sqrt(6)
    y1 = dat1[:,3] * sc
    y2 = 2*dat2[:,3]*dat2[:,1] * sc
    y = y2-y1

    y1avg = np.average(y1[20:180])
    y2avg = np.average(y2[20:180])
    yavg  = np.average(y[20:180])
    yrms  = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

  
    print("E23 = ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  "$E_{23} = " + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'o',   mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
    #iax.plot(x, y1, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    #iax.plot(x, y2, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
    icolor = icolor + 1


    #-- E33 ~ C2211-00 --
    
    iax = ax[0][1]
    icolor = 0;

    rname = "Post02/" + run  + "_2211-00_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_21-0_mcr.txt"
    dat2 = np.loadtxt(rname)

    x  = dat1[:,0]
    sc = 1/dat1[:,2] * np.abs(x-d) / np.sqrt(8)
    y1 = dat1[:,3] * sc
    y2 = dat2[:,3]**2 * sc
    y = y2-y1

    y1avg = np.average(y1[20:180])
    y2avg = np.average(y2[20:180])
    yavg  = np.average(y[20:180])
    yrms  = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

  
    print("E33 = ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  "$E_{33} = " + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'o',   mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
    iax.plot(x, y1, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, y2, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
    icolor = icolor + 1


#-- E24 ~ C3222-20 --
    
    iax = ax[1][0]
    icolor = 0;

    rname = "Post02/" + run  + "_3222-20_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_322-0_mcr.txt"
    dat2 = np.loadtxt(rname)

    nk2 = dat2[:,1] 
    nk2 = np.hstack(( np.zeros(2), nk2[:-2] ))
    
    x  = dat1[:,0]
    sc = 1/dat1[:,2] * np.abs(x-d)/ np.sqrt(24)
    y1 = dat1[:,3] * sc
    y2 = 6*dat2[:,3]*nk2 * sc
    y = y2-y1

    y1avg = np.average(y1[20:180])
    y2avg = np.average(y2[20:180])
    yavg  = np.average(y[20:180])
    yrms  = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

    print("E24 = ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  "$E_{24} = " + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'o',   mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
    iax.plot(x, y1, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, y2, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
    icolor = icolor + 1


#-- E24a ~ C4310-00 --
    
    iax = ax[1][1]
    icolor = 0;

    rname = "Post02/" + run  + "_4310-00_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_431-0_mcr.txt"
    dat2 = np.loadtxt(rname)
    
    x  = dat1[:,0]
    sc = 1/dat1[:,2] * np.abs(x-d)/ np.sqrt(6)
    y1 = dat1[:,3] * sc
    y2 = 2*dat2[:,3]*dat2[:,1] * sc
    y = y2-y1

    y1avg = np.average(y1[20:180])
    y2avg = np.average(y2[20:180])
    yavg  = np.average(y[20:180])
    yrms  = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

    print("E24a = ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  "$E_{24a} = " + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'o',   mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
    iax.plot(x, y1, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, y2, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
    icolor = icolor + 1




#-- E24b ~ C3220-00 --
    
    iax = ax[2][0]
    icolor = 0;

    rname = "Post02/" + run  + "_3220-00_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_322-0_mcr.txt"
    dat2 = np.loadtxt(rname)

    x  = dat1[:,0]
    sc = 1/dat1[:,2] * np.abs(x-d)/ np.sqrt(12)
    y1 = dat1[:,3] * sc
    y2 = 6*dat2[:,3]*dat2[:,1] * sc
    y = y2-y1

    y1avg = np.average(y1[20:180])
    y2avg = np.average(y2[20:180])
    yavg  = np.average(y[20:180])
    yrms  = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

    print("E24b = ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  "$E_{24b} = " + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'o',   mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
    iax.plot(x, y1, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, y2, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
    icolor = icolor + 1


#-- E33b ~ C110-300 --
    
    iax = ax[2][1]
    icolor = 0;

    rname = "Post02/" + run  + "_110-300_mcr.txt"
    dat1 = np.loadtxt(rname)

    rname = "Post02/" + run  + "_11-30_mcr.txt"
    dat2 = np.loadtxt(rname)
    
    x  = dat1[:,0]
    sc = 1/dat1[:,2] * np.abs(x-d)/ np.sqrt(12)
    y1 = dat1[:,3] * sc
    y2 = 2*dat2[:,3]*dat2[:,1] * sc
    y = y2-y1

    y1avg = np.average(y1[20:180])
    y2avg = np.average(y2[20:180])
    yavg  = np.average(y[20:180])
    yrms  = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

    print("E33b = ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
        y1avg, y2avg, yavg, yrms)) 

    label =  "$E_{33b} = " + "{:>+7.4f}".format(yavg) + "$"

    iax.plot(x, y, 'o',   mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
    iax.plot(x, y1, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, y2, ':',  mfc='none', ms=MS, color=colors[icolor], lw=LW)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
    icolor = icolor + 1






    
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

plt.subplots_adjust(left=0.06, right=0.85, top=0.99, bottom=0.09, hspace=0.15, wspace=0.66)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
