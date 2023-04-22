import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('xtick', labelsize=12) 
matplotlib.rc('ytick', labelsize=12) 
matplotlib.rc('axes', labelsize=13)

output=1

outfile = "plot_corr_alpha12_fibo.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"

nmodes=200

corrs = ("21-0", "322-0", "431-0", "11-30")
labels = ("${21\\bar{0}}$", "${322\\bar{0}}$", "${431\\bar{0}}$", "${\\bar{3}11\\bar{0}}$")

ffactor= {
    
    "21-0":   1,
    
    "322-0":  2,
    "431-0":  1,
    "11-30":  2,
    
    "2211-00" : 8,
    "3221-10" : 4,
    "211-320" : 4,
    "4321-20" : 2,
    "221-430" : 2,

}




#---------------------------------------------------------------------------------


#fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.0,4.6))
fig, ax = plt.subplots(ncols=2, figsize=(5.5,2.4))

        
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


runs = ("s5m90", "s5m10")
#runs = ("s4v14", "s4v34")
#runs = ("s4v38", "s4v58")



#-- load and plot correlators vs i --

iax = ax[0]

for run in runs:

    if run in ( 's5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
    
    icolor = 0;
    for c in corrs:

        rname = "Post02/" + run  + "_" + c + "_mcr.txt"
        m=len(c)+1
        r=(m-3)/3

        dat = np.loadtxt(rname)
        nk = dat[:,1]
        x = dat[:,0]
        y = dat[:,3]/dat[:,2] * np.abs(x-d) /  np.sqrt(ffactor[c])

        yavg = np.average(y[20:180])
        yrms = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )


        if run == "s5m90":
            label = labels[icolor]
        else:
            label = ""


        iax.plot(x, y, 'o',  mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
        iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
        icolor = icolor +1

iax.set_ylim(-0.6,0.6)
iax.set(xlabel='$|i - d|$') # ylabel="$D_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
iax.set_xlim(0,nmodes)
iax.legend(frameon=False, handlelength=-0.5, loc="upper left", bbox_to_anchor=(1.01,1.05))
iax.set_yticks(np.arange(-0.6, 0.7, step=0.2))


#-- load and plot correlators vs i --


corrs = ("2211-00", "3221-10", "211-320", "4321-20", "221-430")
labels = ("${2211\\bar{0}\\bar{0} }$",
          "${3221\\bar{1}\\bar{0}}$",
          "${\\bar{3}\\bar{2}211\\bar{0}}$",
          "${432\\bar{2}1\\bar{0}}$",
          "${\\bar{4}\\bar{3}221\\bar{0}}$")

iax = ax[1]

for run in runs:

    if run in ( 's5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
    
    icolor = 4;
    for c in corrs:

        rname = "Post02/" + run  + "_" + c + "_mcr.txt"
        m=len(c)+1
        r=(m-3)/3

        dat = np.loadtxt(rname)
        nk = dat[:,1]
        x = dat[:,0]
        y = dat[:,3]/dat[:,2] * np.abs(x-d) /  np.sqrt(ffactor[c])

        yavg = np.average(y[20:180])
        yrms = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

        if run == "s5m90":
            label = labels[icolor-4]
        else:
            label = ""


        iax.plot(x, y, 'o',  mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
        iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
        icolor = icolor +1


iax.set_ylim(-0.6,0.6)
iax.set(xlabel='$|i - d|$') # ylabel="$D_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
iax.set_xlim(0,nmodes)
iax.legend(frameon=False, handlelength=-0.5, loc="upper left", bbox_to_anchor=(-0.54,0.58))
iax.set_yticks(np.arange(-0.6, 0.7, step=0.2))

        
#------------------------------------------------------------------------------------




plt.subplots_adjust(left=0.07, right=0.98, top=0.97, bottom=0.19, hspace=0.15, wspace=0.55)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
