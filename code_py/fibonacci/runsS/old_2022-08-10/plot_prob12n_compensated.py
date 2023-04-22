import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes




outfile = "plot_prob12n_compensated.pdf"

#-- common parameters --

prefix = "MIx32x1/"

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
    ("s5m10", "20",  "b:"),
    ("s5m10", "160",  "b--"),
    ("s5m10", "180",  "b-"),
)

curves2= (
    ("s5m90", "20",  "b-"),
    ("s5m90", "40",  "b--"),
    ("s5m90", "180", "b:"),
)


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
LW = 1


alpha = 1/2
phi = (1 + np.sqrt(5))/2



fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6,6))


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
        ax[i][0].plot(x,y, style, mfc='none', ms=MS,  lw=LW,  label=label);
        ax[i][1].plot(x,y+x, style, mfc='none', ms=MS,  lw=LW,  label=label);

        
#-- canvas options --

x=np.arange(0,21)
ax[0][0].plot(x,-x, '-k', lw=LW/3)
ax[1][0].plot(x,-x, '-k', lw=LW/3)

ax[0][0].text(1,-18, "direct")
ax[1][0].text(1,-18, "inverse")


for iax in (ax[0][0], ax[1][0]):
    iax.set_xlim(0,20)
    iax.set_ylim(-20,0)
    iax.set(xlabel='$\\rho^2_i / n_i$', ylabel='$ln {\cal P}$')
    iax.set_yticks(np.arange(-20, 1, step=5)) 
    iax.legend(frameon=False)

for iax in (ax[0][1], ax[1][1]):  
    iax.set_xlim(0,20)
    iax.set_ylim(-1,3)
    iax.set(xlabel='$\\rho^2_i / n_i$', ylabel='$\ln {\cal P} + \\rho^2_i / n_i$')
    iax.set_yticks(np.arange(-1, 4, step=1))
    iax.yaxis.labelpad = 0
    iax.legend(frameon=False)


    
#-----------------------

plt.subplots_adjust(left=0.08, right=0.92, top=0.96, bottom=0.1, hspace=0.3, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
