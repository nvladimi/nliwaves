import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

outfile = "plot_corr_triple_cd.pdf"; fbase =  "../cdataD3_4-10_" ;
ymax=((0.10, 0.05), (0.04, 0.03)); ylabel="$\langle a_1 a_2 a_3^* \\rangle $"


#-- canvas options --

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,8))

runs = (("c1", "c5"),("d1", "d5"))

      
for i in (0,1):
    for j in (0,1):


        run=runs[i][j]
        
        fname = fbase + run +".txt"

        print(fname)
        iax = ax[i][j]
        dat = np.loadtxt(fname, usecols=range(1,12))

        w    = dat[:,0]**2 + dat[:,1]**2
        corr = dat[:,2]
        bg   = corr/dat[:,3]
        phi   = dat[:,4]
        w2    = dat[:,5]**2 + dat[:,6]**2
        corr2 = dat[:,7]
        group = dat[:,10]

        ind = (group == 0 )

        iax.plot(w[ind], corr[ind],  'or',  mfc='none', ms=5,  mew=1, label="resonant")
        iax.plot(w[ind], bg[ind],    'og',  mfc='none', ms=5,  mew=1, label="background")
        iax.plot(w[ind], corr2[ind], 'db',  mfc='none', ms=5,  mew=1, label="non resonant")

        ind = (group == 1 )

        iax.plot(w[ind], corr[ind],  '+r',  mfc='none', ms=4,  mew=0.2)
        iax.plot(w[ind], corr2[ind], '+b',  mfc='none', ms=4,  mew=0.2)
 
        ind = (group == 2 )

        iax.plot(w[ind], corr[ind],  'xr',  mfc='none', ms=4,  mew=0.2)
        iax.plot(w[ind], corr2[ind], 'xb',  mfc='none', ms=4,  mew=0.2)


        ind = (group == 3 )

        iax.plot(w[ind], corr[ind],  '+r',  mfc='none', ms=4,  mew=0.2)
        iax.plot(w[ind], corr2[ind], '+b',  mfc='none', ms=4,  mew=0.2)

        ind = (group > 3 )

        iax.plot(w[ind], corr[ind],  'xr',  mfc='none', ms=4,  mew=0.2)
        iax.plot(w[ind], corr2[ind], 'xb',  mfc='none', ms=4,  mew=0.2)

        iax.set_title(run)
        iax.set(xlabel="$|k_1|^2 + |k_2|^2$", ylabel=ylabel)
        iax.legend() #frameon=False)

        iax.set_xlim(25, 175)
        ax[i][j].set_ylim(0, ymax[i][j])



#-----------------------

plt.subplots_adjust(left=0.15, right=0.98, top=0.9, bottom=0.12, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)


plt.close()

   



#===========================================



#=========================================================
