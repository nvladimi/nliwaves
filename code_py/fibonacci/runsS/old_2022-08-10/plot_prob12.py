import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes




outfile = "plot_prob12.pdf"

#-- common parameters --

prefix = "MIx32x1/"

fbasedir = ["s3m10", "s4m10", "s5m10"]
fbaseinv = ["s3m50", "s4m90", "s5m90"]


# modes:
#
# s3m10:  15, 20, 25, 30, 35, 40
# s4m10:  15, 20, 25, 30, 35, 40,  50, 60, 70, 80
# s5m10:      20,     30,     40,  50, 60, 70, 80,   100, 110, 120, 140, 160, 180
#
# s3m50:                                            20, 25, 30, 35, 40, 45
# s4m90:                        20,  30,  40,  50,  60,     70,     80
# s5m90:  20, 40, 60, 80, 100, 120, 130, 140, 150, 160,    170,    180
         

curves1= (
#    ("s3m10", "20",  "g--"),
#    ("s3m10", "40",  "g-"),
    ("s4m10", "20",  "r:"),
    ("s4m10", "60",  "r--"),
    ("s4m10", "80",  "r-"),
    ("s5m10", "20",  "b:"),
    ("s5m10", "160",  "b--"),
    ("s5m10", "180",  "b-"),
)

curves2= (
#    ("s3m50", "20",  "g-"),
#    ("s3m50", "40",  "g--"),
    ("s4m90", "20",  "r-"),
    ("s4m90", "40",  "r--"),
    ("s4m90", "80",  "r:"),
    ("s5m90", "20",  "b-"),
    ("s5m90", "40",  "b--"),
    ("s5m90", "180", "b:"),
)



MS = 0.5
LW = 1


alpha = 1/2
phi = (1 + np.sqrt(5))/2



fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(10,6))


 #-- load and plot data --

for c in range(len(curves1)):
    try:
        base=curves1[c][0]; i=curves1[c][1]; style=curves1[c][2]
        
        fname = prefix + base + "_i" + i + "_prob.txt"
        dat = np.loadtxt(fname)
        x = dat[:,0]
        y = np.log(dat[:,3])
        ax[0][0].plot(x,y, style, mfc='none', ms=MS,  lw=LW,  label=i);
        ax[0][1].plot(x,y+x, style, mfc='none', ms=MS,  lw=LW,  label=i);

        fname = prefix + base + "_i" + i + "_probphi.txt"
        dat = np.loadtxt(fname)
        x = dat[:,0]/np.pi
        y = dat[:,1]*np.pi - 0.5
        ax[0][2].plot(x,y, style, mfc='none', ms=MS,  lw=LW,  label=i);
        
    except:
        print(fname)
        pass


for c in range(len(curves2)):
    try:
        base=curves2[c][0]; i=curves2[c][1]; style=curves2[c][2]
        
        fname = prefix + base + "_i" + i + "_prob.txt"
        dat = np.loadtxt(fname)
        x = dat[:,0]
        y = np.log(dat[:,3])
        ax[1][0].plot(x,y, style, mfc='none', ms=MS,  lw=LW,  label=i);
        ax[1][1].plot(x,y+x, style, mfc='none', ms=MS,  lw=LW,  label=i);

        fname = prefix + base + "_i" + i + "_probphi.txt"
        dat = np.loadtxt(fname)
        x = dat[:,0]/np.pi
        y = dat[:,1]*np.pi - 0.5
        ax[1][2].plot(x,y, style, mfc='none', ms=MS,  lw=LW,  label=i);
        
    except:
        print(fname)
        pass



        
#-- canvas options --

x=np.arange(0,21)
ax[0][0].plot(x,-x, '-k', lw=LW/3)
ax[1][0].plot(x,-x, '-k', lw=LW/3)

ax[0][0].text(9,-2, "DIRECT CASCADE")
ax[1][0].text(9,-2, "INVERSE CASCADE")

#ax[0][0].text(2,-14, "green: 60 modes")
ax[0][0].text(2,-16, "red: 100 modes")
ax[0][0].text(2,-18, "blue: 200 modes")
#ax[1][0].text(2,-14, "green: 60 modes")
ax[1][0].text(2,-16, "red: 100 modes")
ax[1][0].text(2,-18, "blue: 200 modes")



for iax in (ax[0][0], ax[1][0]):
    iax.set_xlim(0,20)
    iax.set_ylim(-20,0)
    iax.set(xlabel='$\\rho^2_i / n_i$', ylabel='$ln {\cal P}$')
    iax.set_yticks(np.arange(-20, 1, step=5)) 
    #iax.legend(frameon=False)

for iax in (ax[0][1], ax[1][1]):  
    iax.set_xlim(0,20)
    iax.set_ylim(-1,3)
    iax.set(xlabel='$\\rho^2_i / n_i$', ylabel='$\ln {\cal P} + \\rho^2_i / n_i$')
    iax.set_yticks(np.arange(-1, 4, step=1))
    iax.yaxis.labelpad = 0
    #iax.legend(frameon=False)

for iax in (ax[0][2], ax[1][2]):  
    iax.set_xlim(-1,1)
    iax.set_ylim(-0.03,0.03)
    iax.set(xlabel='$\\theta_i / \pi$', ylabel='${\cal P}  - 1/2$')
    iax.set_yticks(np.arange(-0.03, 0.031, step=0.01))
    iax.yaxis.labelpad = -2
    iax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    iax.grid()

    
#-----------------------

plt.subplots_adjust(left=0.08, right=0.92, top=0.96, bottom=0.1, hspace=0.3, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
