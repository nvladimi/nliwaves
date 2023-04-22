import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

correlators = "dbl1"

#----------------------
if correlators == "dbl1":

    outfile = "plot_corr_dblnew1.pdf"
    
    corr = ( ("20",  1),
             ("40",  1),
             ("50",  1),
             ("60",  1),
    )


    tiles =( ((-220,220), (-80,80),  (-80,80)),
             ((-60,60),   (-30,30),  (-30,30)),
             ((-30,30),   (-12,12),  (-12,12)),
             ((-30,30),     (-8,8),  (-8,8)))
"""
    
    fits = ( ((1, 2.24),  ( 40.0, -50.0), (20.0, -50.0)),
             ((1, 0.55),  ( 10.0, -12.0), ( 8.0, -15.0)),
             ((1, 0.26),  (  4.0,  -8.0), ( 2.0,  -6.0)),
             ((1, 0.24),  (  4.0,  -7.0), ( 2.0,  -5.0)))

elif correlators == "dbl2":

    outfile = "plot_corr_dbl2.pdf"
    
    corr = ( ("5",  1),
             ("6",  1),
             ("7",  1),
             ("8",  1),
    )
    
    tiles =( ((-20,20),   (-6,6),  (-6,6)),
             ((-15,15),   (-5,5),  (-5,5)),
             ((-8,8),     (-4,4),  (-4,4)),
             ((-6,6),     (-3,3),  (-3,3)))

    fits = ( ((1, 0.165),  ( 3.0, -4.5), ( 2.5, -3.5)),
             ((1, 0.100),  ( 2.2, -2.5), ( 2.0, -1.5)),
             ((1, 0.063),  ( 1.7, -1.5), ( 1.8, -1.6)),
             ((1, 0.037),  ( 1.5, -0.7), ( 1.5, -1.0)))

else:
    quit()

"""    

print(outfile)
#---------------------------------------------------------------------------------

runs = ( ('s4m90', 1/2, 100,  90, '-or'),
         ('s5m90', 1/2, 200, 190, '-oy'),
         ('s4m10', 1/2, 100,  10, '-ob'),
         ('s5m10', 1/2, 200,  10, '-oc'),
         ('s3v38', 3/8,  60,  50, '-or'),
         ('s4v38', 3/8, 100,  90, '-oy'),
         ('s3v58', 5/8,  60,  10, '-ob'),
         ('s4v58', 5/8, 100,  10, '-oc'),
         ('s3v14', 1/4,  60,  50, '-or'),
         ('s4v14', 1/4, 100,  90, '-oy'),
         ('s3v34', 3/4,  60,  10, '-ob'),
         ('s4v34', 3/4, 100,  10, '-oc') )

#--------------------------------




#-- canvas options --

fig, ax = plt.subplots(ncols=3, nrows=4, figsize=(10,12))

ax[0][0].set_title("$\\alpha = 1/2$")
ax[0][1].set_title("$\\alpha - 1/2 = \pm 1/8$")
ax[0][2].set_title("$\\alpha - 1/2 = \pm 1/4$")


for j in range(0,4):
    for i in range(0,3):
        ax[j][i].set(xlabel='$|i-k|$')
        ax[j][i].set_xlim(0,100)
#        ax[j][i].set_yscale('log')
#        ax[j][i].set_xscale('log')
        #ax[j][i].set_ylim(0, 300)
        #ax[j][i].set_ylim(tiles[j][i][0], tiles[j][i][1])
    ax[j][0].set(ylabel="$|qq_{k=" + corr[j][0]  +"}|$")


#-- common parameters --

MS = 2
LW = 0.5
#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, 


def cfactor(alpha,n):    
    phi = (1 + np.sqrt(5))/2
#    return  phi**(1 + alpha)
    return  -phi**((1 + alpha)*2)
#   return  phi**((1 + alpha)*n/3)


#-- load and plot data --

for j in range(0,4):
    
    c = corr[j][0];
    n=5
    for i in range(0,3):
        iax = ax[j][i]
        for r in range(i*4,i*4+4):
            try:
                run = runs[r]
                cf = cfactor(run[1],n)
                rname = "Post3j/" + run[0]  + "_" + c + "_dcn.txt"
                dat = np.loadtxt(rname)
                #x = abs(dat[:,0] - int(c));
                x = dat[:,0];
                y = dat[:,3];
                iax.plot(x, y, run[4],  mfc='none', ms=MS,  lw=LW)
            except:
                pass

#-----------------------

plt.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.05, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
