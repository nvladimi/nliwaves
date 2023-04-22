import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

correlators = "dbl1"

#----------------------
if correlators == "dbl1":

    outfile = "plot_corr_dbl1.pdf"
    
    corr = ( ("1",  1),
             ("2",  1),
             ("3",  1),
             ("4",  1),
    )

    tiles =( ((-0.6,0.6),   (-0.054,0.04),  (-0.2,0.2)),
             ((-0.15,0.15), (-0.015,0.015),  (-0.05,0.05)),
             ((-0.08,0.08), (-0.006,0.006),  (-0.04,0.04)),
             ((-0.08,0.08), (-0.006,0.006),  (-0.02,0.02)))

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


print(outfile)
#---------------------------------------------------------------------------------

runs = ( ('s4m90', 1/2, 100,  90, '-r'),
         ('s5m90', 1/2, 200, 190, '-y'),
         ('s4m10', 1/2, 100,  10, '-b'),
         ('s5m10', 1/2, 200,  10, '-c'),
         ('s3v38', 3/8,  60,  50, '-r'),
         ('s4v38', 3/8, 100,  90, '-y'),
         ('s3v58', 5/8,  60,  10, '-b'),
         ('s4v58', 5/8, 100,  10, '-c'),
         ('s3v14', 1/4,  60,  50, '-r'),
         ('s4v14', 1/4, 100,  90, '-y'),
         ('s3v34', 3/4,  60,  10, '-b'),
         ('s4v34', 3/4, 100,  10, '-c') )

#--------------------------------




#-- canvas options --

fig, ax = plt.subplots(ncols=3, nrows=4, figsize=(10,12))

ax[0][0].set_title("$\\alpha = 1/2$")
ax[0][1].set_title("$\\alpha - 1/2 = \pm 1/8$")
ax[0][2].set_title("$\\alpha - 1/2 = \pm 1/4$")


for j in range(0,4):
    for i in range(0,3):
        ax[j][i].set_ylim(tiles[j][i][0], tiles[j][i][1])
    ax[j][0].set(xlabel='$i$', ylabel="$JJ_{" + corr[j][0]  +"}$")
    ax[j][0].set_xlim(0,200)
    ax[j][1].set_xlim(0,100)
    ax[j][2].set_xlim(0,100)
      
#-- straight lines --

"""
x = np.arange(0,201);
LW = 0.5

for j in range(0,4):
    
    a = fits[j][0][0]; b = fits[j][0][1];
    ax[j][0].plot(x, a + x*b, '--k',  mfc='none', lw=LW)
    ax[j][0].text(101, 0, "{}".format(b))
    
    a = fits[j][1][0]; b = fits[j][1][1];
    ax[j][1].plot(x, a + x*0, '--r',  mfc='none', lw=LW)
    ax[j][1].plot(x, b + x*0, '--b',  mfc='none', lw=LW)
    ax[j][1].text(101, a, "{}".format(a))
    ax[j][1].text(101, b, "{}".format(b))
    
    a = fits[j][2][0]; b = fits[j][2][1];
    ax[j][2].plot(x, a + x*0, '--r',  mfc='none', lw=LW)
    ax[j][2].plot(x, b + x*0, '--b',  mfc='none', lw=LW)
    ax[j][2].text(101, a, "{}".format(a))
    ax[j][2].text(101, b, "{}".format(b))
"""
        
#-- common parameters --

MS = 0.5
LW = 1
#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, 


def cfactor(alpha,n):    
    phi = (1 + np.sqrt(5))/2
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
                #cf = cfactor(run[1],n)
                rname = "Post2j/" + run[0]  + "_" + c + "_dcr.txt"
                dat = np.loadtxt(rname)
                x = dat[:,0];
                y = dat[:,3]/dat[:,2];
                if (run[0]=="s4m90" or run[0]=="s5m90"):
                    y = y * x
                if (run[0]=="s4m10"):
                    y = y * (100-x)
                if (run[0]=="s5m10"):
                    y = y * (200-x)
                iax.plot(x, y, run[4],  mfc='none', ms=MS,  lw=LW)
            except:
                pass

#-----------------------

plt.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.05, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
