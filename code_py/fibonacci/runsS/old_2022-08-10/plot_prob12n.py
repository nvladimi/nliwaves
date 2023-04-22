import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
matplotlib.rc('axes', labelsize=14)


outfile = "plot_prob12n.pdf"

#-- common parameters --


# modes:
#
# s3m10:  15, 20, 25, 30, 35, 40
# s4m10:  15, 20, 25, 30, 35, 40,  50, 60, 70, 80
# s5m10:      20,     30,     40,  50, 60, 70, 80,   100, 110, 120, 140, 160, 180
#
# s3m50:                                            20, 25, 30, 35, 40, 45
# s4m90:                        20,  30,  40,  50,  60,     70,     80
# s5m90:  20, 40, 60, 80, 100, 120, 130, 140, 150, 160,    170,    180
         

runs= ("MIx32x1/s5m90", "MIx32x1/s5m10")

modes= ((("20", "$\,20$", "-r"), 
         ("40", "$\,40$", "-b"), 
         ("80", "$\,80$", "--c"), 
         ("120","$120$","--y"), 
         ("160","$160$", ":m"), 
         ("180","$180$", ":g")), 
        (("180", " $20$", "-r"), 
         ("160", " $40$", "-b"), 
         ("120", " $80$", "--c"), 
         ("80", "$120$", "--y"), 
         ("40", "$160$", ":m"), 
         ("20", "$180$", ":g")) 
)
 

#---------------------------------




MS = 0.5
LW = 1.5


alpha = 1/2
phi = (1 + np.sqrt(5))/2



fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(6,2.8))


for i in (0,1):

    #iax = ax[0][i]

    run = runs[i]

    for r in modes[i]:

        fname = run + "_i" + r[0] + "_prob.txt"
        
        label = r[1]
        style = r[2]
        dat = np.loadtxt(fname)
        
        x = dat[:,0]
        y = np.log(dat[:,3])
        ax[i].plot(x,y, style, mfc='none', ms=MS,  lw=LW,  label=label);

        
#-- canvas options --

x=np.arange(0,21)
ax[0].plot(x,-x, '-k', lw=LW/3)
ax[1].plot(x,-x, '-k', lw=LW/3)

ax[0].text(1.5, -18.5, "inverse", fontsize=14)
ax[1].text(1.5, -18.5, "direct",  fontsize=14)

for iax in (ax[0], ax[1]):
    iax.set_xlim(0,20)
    iax.set_ylim(-20,0)
    iax.set(xlabel='$|b_i|^2 / \langle |b_i|^2 \\rangle $')
    iax.set_yticks(np.arange(-20, 1, step=5)) 
    #iax.legend(frameon=False)

ax[1].text(14, -1.5, "$|i-d|$", fontsize=12)
ax[1].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.94), labelspacing = 0.4)

    
ax[0].set(ylabel='$\ln {\cal P}$')

    
#-----------------------

plt.subplots_adjust(left=0.11, right=0.96, top=0.98, bottom=0.18, wspace=0.23)


plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
