import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "mirror-returned-graining.pdf"

#-- common parameters --

MS = 2
LW = 1
  
#---------------------------------



fig, ax = plt.subplots(ncols=5, nrows=1, figsize=(15,2.8))


fname = "../2022-06-19_graining/data/lega.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]/2
y1 = dat[:,8]
y2 = dat[:,9]
y3 = dat[:,10]
y4 = dat[:,14]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0].plot(x,y1, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[1].plot(x,y2, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[2].plot(x,y3, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[3].plot(x,y4, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[4].plot(x,y5a,'-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[4].plot(x,y5b,'--ob',                       mfc='none', ms=MS,  lw=LW);



fname = "../2022-06-19_graining/data/zero.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]/2
y1 = dat[:,8]
y2 = dat[:,9]
y3 = dat[:,10]
y4 = dat[:,14]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0].plot(x,y1, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[1].plot(x,y2, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[2].plot(x,y3, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[3].plot(x,y4, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[4].plot(x,y5a,'-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[4].plot(x,y5b,'--or',                   mfc='none', ms=MS,  lw=LW);




        
#-- canvas options --

#ax[0][0].text(1.5, 0.73, "legacy parameters, atmosphere", fontsize=12)
#ax[0][1].text(2.5, 0.52, "legacy parameters, vacuum",  fontsize=12)

for iax in (ax[0], ax[1], ax[2], ax[3], ax[4]):
    iax.set_xlim(0,10)
    iax.set_ylim(0,1.0)
    iax.set(xlabel='$\sigma / \Delta x$')
    iax.legend(frameon=False, labelspacing = 0.3, loc="lower right")

ax[3].set_ylim(0,70)
ax[4].set_ylim(0,100)
ax[2].legend(frameon=False, labelspacing = 0.3, loc="upper right")
ax[3].legend(frameon=False, labelspacing = 0.3, loc="upper left")
ax[4].legend(frameon=False, labelspacing = 0.3, loc="upper right")

    
ax[0].set(ylabel='$|\Delta E|^2$')
ax[1].set(ylabel='$|\Delta I|$')
ax[2].set(ylabel='$P_{\\rm diff}/P$')
#ax[3].set(ylabel='$w$')
ax[3].set(ylabel='$M^2$')
ax[4].set(ylabel='$r(0.5P)$,  $r(0.9P)$')


#ax[1][1].legend(frameon=False, labelspacing = 0.3)


#iax.set_yticks(np.arange(-20, 1, step=5))
#ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)


    
#-----------------------

plt.subplots_adjust(left=0.04, right=0.99, top=0.98, bottom=0.14, wspace=0.26, hspace=0.26)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
