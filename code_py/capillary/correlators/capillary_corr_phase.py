import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

outfile = "capillary_corr_phase.pdf"; fbase =  "../cdataD3_4-10_" ;
ymax=(0.10, 0.10, 0.05); ylabel="$\langle a_1 a_2 a_3^* \\rangle $"

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
matplotlib.rc('axes', labelsize=11)

#-- canvas options --

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(4,2))

runs = ("c1", "c5")
nicknames = ("Low forcing", "High forcing")



for i in (0,1):
 

    run=runs[i]

    fname = fbase + run +".txt"

    print(fname)
    iax = ax[i]
    dat = np.loadtxt(fname, usecols=range(1,12))

    w    = dat[:,0]**2 + dat[:,1]**2
    corr = dat[:,2]
    bg   = corr/dat[:,3]
    phi   = dat[:,4]
    group = dat[:,10]

    
    ind = (group == 0 )
    iax.plot(w[ind], phi[ind],  'or',  mfc='none', ms=5,  mew=1.25)

    iax.set(xlabel="$|k_1|^2 + |k_2|^2$")
    iax.yaxis.set_label_coords(-0.17, 0.5)
    iax.set_xlim(25, 175)
    iax.set_ylim(-1, 1)
    #iax.set_title(nicknames[i] + ":  $C_{res}$" )

ax[0].text(90,0.8, "low forcing", fontsize=12)
ax[1].text(90,0.8, "high forcing", fontsize=12)
ax[0].set(ylabel="$\\theta / \pi$")

    
#-----------------------

plt.subplots_adjust(left=0.11, right=0.99, top=0.96, bottom=0.21, wspace=0.29) # hspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()


#=========================================================
