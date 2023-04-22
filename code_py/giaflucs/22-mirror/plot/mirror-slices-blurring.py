import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "mirror-slices-blurring.pdf"

#-- common parameters --

MS = 0.5
LW = 0.5

runs = (
    ("00", 0, "$\sigma/\Delta x =  0$"),
    ("02", 1, "$\sigma/\Delta x =  1$"),
    ("04", 2, "$\sigma/\Delta x =  2$"),
    ("08", 3, "$\sigma/\Delta x =  4$"),
)

        
#---------------------------------



fig, ax = plt.subplots(ncols=4, nrows=4, figsize=(12,12))


for rset in runs:
    run = rset[0]
    row = rset[1]
    lbl = rset[2]
    
    fname = "../2022-04-10_blurring/data/lega" + run + "_last.txt" 
    print(fname)
    
    dat = np.loadtxt(fname,  comments='%')
        
    x = dat[:,0]
    y1 = dat[:,1]**2 + dat[:,2]**2
    y2 = dat[:,3]**2 + dat[:,4]**2
    
    ax[row][0].plot(x,y1, '-r',  mfc='none', ms=MS, lw=LW, label=lbl) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[row][0].plot(x,y2, '-b',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[row][1].plot(x,y1, '-r',  mfc='none', ms=MS, lw=LW, label=lbl) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[row][1].plot(x,y2, '-b',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);

    
    fname = "../2022-04-10_blurring/data/zero" + run + "_last.txt" 
    print(fname)
    
    dat = np.loadtxt(fname,  comments='%')
        
    x = dat[:,0]
    y1 = dat[:,1]**2 + dat[:,2]**2
    y2 = dat[:,3]**2 + dat[:,4]**2
    
    ax[row][2].plot(x,y1, '-r',  mfc='none', ms=MS, lw=LW, label=lbl) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[row][2].plot(x,y2, '-b',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[row][3].plot(x,y1, '-r',  mfc='none', ms=MS, lw=LW, label=lbl) #, style, mfc='none', ms=MS,  lw=LW,  label=label);
    ax[row][3].plot(x,y2, '-b',  mfc='none', ms=MS, lw=LW) #, style, mfc='none', ms=MS,  lw=LW,  label=label);

    


        
#-- canvas options --

for i in (0,1,2,3):
    ax[i][1].set_yscale('log')
    ax[i][3].set_yscale('log')
    
    ax[i][1].set_ylim(1e-8,1)
    ax[i][3].set_ylim(1e-8,1)
    ax[i][1].set_xlim(-120,120)   
    ax[i][3].set_xlim(-120,120)
    
    #ax[i][0].set_ylim(1e-8,1)
    #ax[i][2].set_ylim(1e-8,1)
    ax[i][0].set_xlim(-10,10)
    ax[i][2].set_xlim(-10,10)

    for j in (0,1,2,3):
        ax[i][j].legend(frameon=False, handlelength=0)




#iax.set_yticks(np.arange(-20, 1, step=5))
#ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)


    
#-----------------------

plt.subplots_adjust(left=0.04, right=0.98, top=0.98, bottom=0.04, wspace=0.25, hspace=0.25)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
