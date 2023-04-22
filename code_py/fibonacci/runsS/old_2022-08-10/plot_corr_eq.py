import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

correlators = "eq1"

#------------

if correlators == "eq1":

    outfile = "plot_corr_eq1.pdf"
    
    corr = ( ("_21",   1,  0,   1),
             ("_21",   3,  (2+1)/3, 3 ),
             ("_322",  3,  (3+2+2)/3, 4),
             ("_431",  3,  (4+3+1)/3, 4 ) )
    
    tiles =( ((0,100, 0, 6,    "$C_0$"),     (0,100, 0,10, ""),   (0,100,  0, 5, "")),
             ((0,100,-0.1,0.1, "$C_{21}$"),  (0,100,-1.2,1.2, ""),  (0,100,-1.2,1.2, "")),
             ((0,100,-0.1,0.1, "$C_{322}$"), (0,100,-6,6, ""),  (0,100,-6,6, "")),
             ((0,100,-0.1,0.1, "$C_{431}$"), (0,100,-3,3, ""),  (0,100,-3,3, "")))


    
#--------------------------------

runs = ( ('eq12', 1/2, 0, 0, '-r'),
         ('eq38', 3/8, 0, 0, '-y'),
         ('eq58', 5/8, 0, 0, '-c') )

#--------------------------------




#-- canvas options --

fig, ax = plt.subplots(ncols=3, nrows=4, figsize=(10,12))

ax[0][0].set_title("$\\alpha = 1/2$")
ax[0][1].set_title("$\\alpha - 1/2 = \pm 1/8$")
ax[0][2].set_title("$\\alpha - 1/2 = \pm 1/4$")

for j in range(0,4):
    for i in range(0,3):

        iax = ax[j][i]  
        iax.set_xlim(tiles[j][i][0], tiles[j][i][1])
        iax.set_ylim(tiles[j][i][2], tiles[j][i][3])
        iax.set(xlabel='$i$', ylabel=tiles[j][i][4])
        #iax.legend(frameon=False)

        
#-- common parameters --

MS = 0.5
LW = 1
#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, 


def cfactor(alpha,n):
    
    phi = (1 + np.sqrt(5))/2
    return  phi**((1 + alpha)*n)


#-- load and plot data --

for j in range(0,4):
    
    c = corr[j][0]; ind=corr[j][1];  n=corr[j][2]; 
    for i in range(0,3):
        iax = ax[j][i]
        for r in range(i*4,i*4+4):
            try:
                run = runs[r]
                cf = cfactor(run[1],n)
                rname = "PostEq/" + run[0]  + c + "_mcr.txt"
                #print(c, cf, rname)
                dat = np.loadtxt(rname)
                x = dat[:,0];
                y = dat[:,ind]*cf;
                iax.plot(x, y, run[4],  mfc='none', ms=MS,  lw=LW)
            except:
                pass

#-----------------------

plt.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.05, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
