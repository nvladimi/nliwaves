import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=12)


outfile = "plot_corr_alpha12_E33_vs_shift.pdf"
datfile = "plot_corr_alpha12_E33_vs_shift.txt"

#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"



#---------------------------------------------------------------------------------


fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(3.4,3.0))

        
#-- common parameters --

MS = 6
LW = 1


#-- load and plot data --

   
      
iax = ax

dat = np.loadtxt(datfile)


#       INVERSE             DIRECT             INVERSE            DIRECT
#
# j   J(i)J(i-j)          J(i)J(i-j)         J(i)J'(i-j)       J(i)J'(i-j)


j         = dat[:,0]
JJinv     = dat[:,1]
JJinverr  = dat[:,2]
JJdir     = dat[:,3]
JJdirerr  = dat[:,4]
JSinv     = dat[:,5]
JSinverr  = dat[:,6]
JSdir     = dat[:,7]
JSdirerr  = dat[:,8]


JJinvalt  = dat[:,9]
JJdiralt  = dat[:,10]
JSinvalt  = dat[:,11]
JSdiralt  = dat[:,12]



#print(c1 + ":  ({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-  {:>+7.4f} ".format(
#    y1avg, y2avg, yavg, yrms)) 
#label =  c1 + "{:>+7.4f}".format(yavg) + "$"

#iax.errorbar(j, np.abs(JJinv),  yerr=JJinverr, color="r", ls=":", mfc='none', ms=MS, lw=LW, label="$JJ$ inv")
#iax.errorbar(j, np.abs(JJdir),  yerr=JJdirerr, color="b", ls=":", mfc='none', ms=MS, lw=LW, label="$JJ$ dir")
#iax.errorbar(j, np.abs(JSinv),  yerr=JSinverr, color="r", ls="-", mfc='none', ms=MS, lw=LW, label="$JJ^*$ inv")
#iax.errorbar(j, np.abs(JSdir),  yerr=JSdirerr, color="b", ls="-", mfc='none', ms=MS, lw=LW, label="$JJ^*$ dir")

iax.plot(j, np.abs(JJinv),  "r:s", mfc='none', ms=MS, lw=LW, label="$JJ$ inverse")
iax.plot(j, np.abs(JJdir),  "b:s", mfc='none', ms=MS, lw=LW, label="$JJ$ direct")
iax.plot(j, np.abs(JSinv),  "r-o", mfc='none', ms=MS, lw=LW, label="$JJ^*$ inverse")
iax.plot(j, np.abs(JSdir),  "b-o", mfc='none', ms=MS, lw=LW, label="$JJ^*$ direct")

iax.plot(j, np.abs(JJinvalt),  "rs", mfc='none', ms=4, lw=LW)
iax.plot(j, np.abs(JJdiralt),  "bs", mfc='none', ms=4, lw=LW)
iax.plot(j, np.abs(JSinvalt),  "ro", mfc='none', ms=4, lw=LW)
iax.plot(j, np.abs(JSdiralt),  "bo", mfc='none', ms=4, lw=LW)



i = np.arange(-1,9)
iax.plot(i, 0.4*np.exp(-i/2), '-k',    mfc='none', ms=MS, lw=LW/2, label="$0.4 \, e^{-j/2}$" )

iax.set_yscale('log')
iax.set(xlabel='$j$') 
iax.set(ylabel='$\langle J_i J_{i-j} \\rangle ,   \quad  \langle J_i J^*_{i-j} \\rangle $') 
iax.legend(frameon=True)
iax.set_xticks(np.arange(0,8))
iax.set_xlim(-0.3,7.3)
iax.set_ylim(8e-3,4e-1)
    

#-----------------------

plt.subplots_adjust(left=0.20, right=0.98, top=0.98, bottom=0.15, hspace=0.15, wspace=0.15)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
