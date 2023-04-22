import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

outfile = "plot_corr_rel.pdf"



#-- canvas options --

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,8))

runs = (("c1", "c5"),("d1", "d5"))

      

for i in (0,1):
    for j in (0,1):


        run=runs[i][j]
        
        fname = "../cdataD3_4-10_" + run +".txt"

        print(fname)
        iax = ax[i][j]
        dat = np.loadtxt(fname, usecols=range(1,11))

        w    = dat[:,0]**2 + dat[:,1]**2
        corr = dat[:,2]
        bg   = corr/dat[:,3]
        phi   = dat[:,4]
        w2    = dat[:,5]**2 + dat[:,6]**2
        corr2 = dat[:,7]

        iax.plot(w, corr/bg, 'or',  mfc='none', ms=5,  lw=1, label="resonant")
        #iax.plot(w, bg,   'og',  mfc='none', ms=5,  lw=1, label="background")
        iax.plot(w, corr2/bg, 'db',  mfc='none', ms=5,  lw=1, label="non resonant")

        iax.set_title(run)
        iax.set(xlabel="$|k_1|^2 + |k_2|^2$", ylabel="$\langle a_1 a_2 a_3^* \\rangle  / \langle a_1 a_2 \\bar{a}_3^* \\rangle  $")
        iax.legend(frameon=False)

        iax.set_xlim(25, 175)


ax[0][0].set_ylim(0, 175)
ax[0][1].set_ylim(0, 75)
ax[1][0].set_ylim(0, 75)
ax[1][1].set_ylim(0, 40)




#-----------------------

plt.subplots_adjust(left=0.15, right=0.98, top=0.9, bottom=0.12, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)


plt.close()

   



#===========================================



#=========================================================
