import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "mirror-forth-n-back-grain15cm.pdf"

#-- common parameters --

MS = 0.5
LW = 0.5

runs = (
    ("00", "r--", "$\sigma/\Delta x =  0$"),
    ("02", "r--", "$\sigma/\Delta x =  1$"),
    ("04", "r--", "$\sigma/\Delta x =  2$"),
    ("08", "r--", "$\sigma/\Delta x =  4$"),
    ("16", "r--", "$\sigma/\Delta x =  8$"),
    ("24", "r--", "$\sigma/\Delta x = 12$"),
)

        
#---------------------------------



fig, ax = plt.subplots(ncols=2, nrows=4, figsize=(6,12))


for rset in runs:
    run = rset[0]
    fname = "../2022-05-05_grain15cm/data/lega" + run + "/lega" + run + ".000000" 
    print(fname)
    
    dat = np.loadtxt(fname,  comments='%')
        
    x = dat[:,0]
    y1 = dat[:,8]
    y2 = dat[:,9]
    y3 = dat[:,10]
    y4 = dat[:,14]
    ax[0][0].plot(x,y1, '-o',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[1][0].plot(x,y2, '-o',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[2][0].plot(x,y3, '-o',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[3][0].plot(x,y4, '-o',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);


 
for rset in runs:
    run = rset[0]
    fname = "../2022-05-05_grain15cm/data/zero" + run + "/zero" + run + ".000000" 
    print(fname)
    label = rset[2]
    
    dat = np.loadtxt(fname,  comments='%')
        
    x = dat[:,0]
    y1 = dat[:,8]
    y2 = dat[:,9]
    y3 = dat[:,10]
    y4 = dat[:,14]
    ax[0][1].plot(x,y1, '-o',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[1][1].plot(x,y2, '-o',  mfc='none', ms=MS, lw=LW, label=label) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[2][1].plot(x,y3, '-o',  mfc='none', ms=MS, lw=LW, label=label) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[3][1].plot(x,y4, '-o',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);


        
#-- canvas options --

ax[0][0].text(1.5, 1.56, "15 cm beam, atmosphere", fontsize=12)
ax[0][1].text(2.5, 0.021, "15 cm beam, vacuum",  fontsize=12)



for iax in (ax[0][0], ax[0][1], ax[1][0], ax[1][1], ax[2][0], ax[2][1], ax[3][0], ax[3][1]):
    iax.set_xlim(0,20)
    iax.set_ylim(0,1.0)
    iax.set(xlabel='screen')

ax[0][0].set_ylim(0,1.5)
ax[0][1].set_ylim(0,0.02)
ax[1][1].set_ylim(0,0.15)

ax[3][0].set_ylim(0,50)
ax[3][1].set_ylim(0,10)



for iax in (ax[0][0], ax[0][1]):
    iax.set(ylabel='$|\Delta E|^2$')

for iax in (ax[1][0], ax[1][1]):
    iax.set(ylabel='$|\Delta I|$')

for iax in (ax[2][0], ax[2][1]):
    iax.set(ylabel='$P_{\\rm diff}/P$')
    
for iax in (ax[3][0], ax[3][1]):
    iax.set(ylabel='$M^2$')


ax[2][1].legend(frameon=False, labelspacing = 0.3)


#iax.set_yticks(np.arange(-20, 1, step=5))
#ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)


    
#-----------------------

plt.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.08, wspace=0.25, hspace=0.25)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
