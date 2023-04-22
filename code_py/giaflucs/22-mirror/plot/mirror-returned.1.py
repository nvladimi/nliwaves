import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "mirror-returned.pdf"

#-- common parameters --

MS = 2
LW = 1
  
fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(9,6))

#------------------------------------------------------

fname = "../2022-06-19_grainphi/data/lega.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]/2
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0][0].plot(x,y2, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[0][1].plot(x,y3, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5a,'-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5b,'--ob',                       mfc='none', ms=MS,  lw=LW);



fname = "../2022-06-19_grainphi/data/zero.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]/2
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[0][0].plot(x,y2, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[0][1].plot(x,y3, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5a,'-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[0][2].plot(x,y5b,'--or',                   mfc='none', ms=MS,  lw=LW);


fname = "../2022-06-19_grainphi15/data/lega.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]/2
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[1][0].plot(x,y2, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[1][1].plot(x,y3, '-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5a,'-ob', label = 'atmosphere',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5b,'--ob',                       mfc='none', ms=MS,  lw=LW);



fname = "../2022-06-19_grainphi15/data/zero.txt"
dat = np.loadtxt(fname,  comments='%')
        
x = dat[:,0]/2
y2 = dat[:,9]
y3 = dat[:,10]
y5a = dat[:,15]
y5b = dat[:,16]

ax[1][0].plot(x,y2, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[1][1].plot(x,y3, '-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5a,'-or', label = 'vacuum',  mfc='none', ms=MS,  lw=LW);
ax[1][2].plot(x,y5b,'--or',                   mfc='none', ms=MS,  lw=LW);




        
#-- canvas options --

#ax[0][0].text(1.5, 0.73, "legacy parameters, atmosphere", fontsize=12)
#ax[0][1].text(2.5, 0.52, "legacy parameters, vacuum",  fontsize=12)


for iax in (ax[0][0], ax[0][1], ax[0][2], ax[1][0], ax[1][1], ax[1][2]):
    iax.set_xlim(0,10)
    iax.set_ylim(0,1.0)
    iax.set(xlabel='$\sigma / \Delta x$')
    iax.legend(frameon=False, labelspacing = 0.3, loc="upper right")
    
for i in (0,1):
    ax[i][2].set_ylim(0,100)
    ax[i][0].set(ylabel='$|\Delta I|$')
    ax[i][1].set(ylabel='$P_{\\rm diff}/P$')
    ax[i][2].set(ylabel='$r(0.5P)$,  $r(0.9P)$')

for iax in (ax[0][0], ax[1][1]):
    iax.legend(frameon=False, labelspacing = 0.3, loc="lower right")
    
    

#ax[1][1].legend(frameon=False, labelspacing = 0.3)


#iax.set_yticks(np.arange(-20, 1, step=5))
#ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)
    
#-----------------------

plt.subplots_adjust(left=0.08, right=0.99, top=0.98, bottom=0.10, wspace=0.26, hspace=0.26)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
