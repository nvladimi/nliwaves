import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


outfile = "plot_corr_quad_cd.pdf"; fbase =  "../cdataD4_4-10_"  ;
ymax=((0.015, 0.01), (0.03, 0.02)); ylabel="$\langle a_1 a_2 a_3 a_4^* \\rangle $"


#-- canvas options --

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,8))

runs = (("c1", "c5"),("d1", "d5"))
    
for i in (0,1):
    for j in (0,1):


        run=runs[i][j]
        
        fname = fbase + run +".txt"

        print(fname)
        iax = ax[i][j]
        dat = np.loadtxt(fname, usecols=range(1,11))

        w    = dat[:,0]**2 + dat[:,1]**2
        corr = dat[:,2]
        bg   = corr/dat[:,3]
        phi   = dat[:,4]
        w2    = dat[:,5]**2 + dat[:,6]**2
        corr2 = dat[:,7]

        iax.plot(w, corr,  'or',  mfc='none', ms=5,  mew=1, label="resonant")
        iax.plot(w, bg,    'og',  mfc='none', ms=5,  mew=1, label="background")
        iax.plot(w, corr2, 'db',  mfc='none', ms=5,  mew=1, label="non resonant")

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
