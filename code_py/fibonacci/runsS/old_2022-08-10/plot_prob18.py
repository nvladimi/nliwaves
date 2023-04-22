import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes




outfile = "plot_prob18.pdf"

#-- common parameters --

prefix = "MIx32x1/"

fbasedir = ["s3v58", "s4v58"]
fbaseinv = ["s3v38", "s4v38"]


# modes:
#
# s3v58:  20, 40, 50
# s4v58:  20, 40, 60, 80, 90
#
# s3v38:  10, 20, 40
# s4v38:  10, 20, 40, 60, 80

         

curves1= (
    ("s3v58", "20",  "g-"),
    ("s3v58", "40",  "g--"),
    ("s3v58", "50",  "g:"),
    ("s4v58", "20",  "r-"),
    ("s4v58", "40",  "r--"),
    ("s4v58", "90",  "r:"),
)


curves2= (
    ("s3v38", "10",  "g:"),
    ("s3v38", "20",  "g--"),
    ("s3v38", "40",  "g-"),
    ("s4v38", "10",  "r:"),
    ("s4v38", "60",  "r--"),
    ("s4v38", "80",  "r-"),
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

ax[0][0].text(9,-2, "$\\alpha = 5/8$  (direct)")
ax[1][0].text(9,-2, "$\\alpha = 3/8$  (inverse)")

ax[0][0].text(2,-16, "green: 60 modes")
ax[0][0].text(2,-18, "red: 100 modes")
#ax[0][0].text(2,-18, "blue: 200 modes")
ax[1][0].text(2,-16, "green: 60 modes")
ax[1][0].text(2,-18, "red: 100 modes")
#ax[1][0].text(2,-18, "blue: 200 modes")



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
