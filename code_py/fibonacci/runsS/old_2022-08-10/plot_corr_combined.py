import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes


output=1

#outfile = "plot_corr_combined.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"


if output == 1:

    corr = ( ("21",),
             ("322", "431"),
             ("5544443",),
             ("55544444",),
             )

    ymax = ((1, 1, 1),
            (1, 1, 2),
            (100, 100, 400),
            (600, 1500, 8000),
            )

    outfile = "plot_corr_combined1.pdf"

    
elif  output == 2:

    corr = (
        ("4332",  "5422", "5441", "6531"),   
        ("44333","54432","65332","65522"),
        ("65541","76422","76441","76631"),
        ("87531", "63333") 
    )

    ymax = ((1, 2, 2),
            (3, 6, 12),
            (3, 6, 12),
            (3, 6, 12),
    )

    outfile = "plot_corr_combined2.pdf"

    
elif output == 3:
    
    corr = (
        ("544433", "554442", "654333", "655432"),
        ("665551", "764432", "766332", "766522"),
        ("766541", "875332", "875522", "877422"),
        ("875541", "877441", "877631"),
    )

    ymax = ((8, 20, 60),
            (15, 30, 150),
            (12, 20, 150),
            (15, 30, 150),
    )

    outfile = "plot_corr_combined3.pdf"
    
else:
    quit()


#---------------------------------------------------------------------------------


fig, ax = plt.subplots(ncols=3, nrows=4, figsize=(10,12))

        
#-- common parameters --

MS = 1
LW = 1

colors = ['#1f77b4',
          '#ff7f0e',
          '#2ca02c',
          '#d62728',
          '#9467bd',
          '#8c564b',
          '#e377c2',
          '#7f7f7f',
          '#bcbd22',
          '#17becf',
          '#1a55FF']

def cfactor(alpha,n):
    
    phi = (1 + np.sqrt(5))/2
    return  -phi**((1 + alpha)*n/3)



runs = ( ('s5m90', 1/2),
         ('s5m10', 1/2),
         ('s4v38', 3/8),
         ('s4v58', 5/8),
         ('s4v14', 1/4),
         ('s4v34', 3/4) )


#-- load and plot data --


for j in (0,1,2,3):
    color=0
    for  c in corr[j]:
        color=color+1
        mc=colors[color]
        n=0
        for cc in c:
            n += int(cc)

        for i in range(0,3):
            iax = ax[j][i]
            for direction in (0,1):
                
                run = runs[i*2 + direction]
                alpha=run[1]
                rname=run[0]
                
                cf = cfactor(alpha,n)
                rname = "Post1j/" + rname  + "_" + c + "_mcr.txt"
                m=len(c)+1
                r=(m-3)/3
                #print(c, cf, rname)
                dat = np.loadtxt(rname)
                nk = dat[:,1]
                x = dat[:,0]
                y = dat[:,3]*cf
                y = y / nk**(0.5*m-1)

                if direction == 0:
                    label = c
                else:
                    label = ""
                
                iax.plot(x, y, 'o', c=mc,  mfc='none', ms=MS,  lw=LW, label=label)


#-- canvas options --


ax[0][0].set_title("$\\alpha = 1/2$")
ax[0][1].set_title("$\\alpha - 1/2 = \pm 1/8$")
ax[0][2].set_title("$\\alpha - 1/2 = \pm 1/4$")


for j in (0,1,2,3):
    for i in range(0,3):
        ax[j][i].set_ylim(-ymax[j][i],ymax[j][i])
    #ax[j][0].set(xlabel='$i$', ylabel="$C_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
    ax[j][0].set_xlim(0,200)
    ax[j][1].set_xlim(0,100)
    ax[j][2].set_xlim(0,100)
    ax[j][0].legend(frameon=False, handlelength=1, loc="center left", bbox_to_anchor=(-0.5,0.5))

    
#if (correlators == "hi1"):
#    ax[0][0].set_xlim(0,200)
#    ax[0][0].set(xlabel='$i$', ylabel="$C_0$")


#-----------------------

plt.subplots_adjust(left=0.12, right=0.99, top=0.95, bottom=0.05, hspace=0.2, wspace=0.2)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
