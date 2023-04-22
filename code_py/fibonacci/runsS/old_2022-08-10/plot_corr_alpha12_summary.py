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

outfile = "plot_corr_alpha12_summary.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"

nmodes=200

corrs = ("21", "322", "431", "4332")



ffactor= {
    
    "21":   1,
    
    "322":  2,
    "431":  1,
    
    "4332": 2,
    "5422": 2,
    "5441": 2,
    "6531": 1,

    "44333": 12,
    "54432": 2,
    "63333": 24,
    "65332": 2,
    "65522": 4,
    "65541": 2,

    "76422": 2,
    "76441": 2,
    "76631": 2,
    "87531": 1,

    "544433": 12,
    "554442": 12,
    "654333": 6,
    "655432": 2,
    "665551": 12,

    "764432": 2,
    "766332": 4,
    "766522": 4,
    "766541": 2,

    "875332": 2,
    "875522": 4,
    "877422": 4,
    "875541": 2,
    "877441": 4,
    "877631": 2,

    "5544443": 48,
    "55544444": 720
}




#---------------------------------------------------------------------------------


#fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.0,4.6))
fig, ax = plt.subplots(ncols=2, figsize=(6,2.5))

        
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

def cfactor(c):
    alpha=1/2
    phi = (1 + np.sqrt(5))/2
    n = 0
    for cc in c:
        n += int(cc)
    #return  -phi**((1 + alpha)*n/3)
    return n

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

        rname = "Post1j/" + run  + "_" + c + "_mcr.txt"
        m=len(c)+1
        r=(m-3)/3

        dat = np.loadtxt(rname)
        nk = dat[:,1]
        x = dat[:,0]
        y = dat[:,3]/dat[:,2] * np.abs(x-d) / ffactor[c]

        yavg = np.average(y[20:180])
        yrms = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

        print(c, m, cfactor(c), ffactor[c], yavg, yrms)


        #label =  "$D_{"  + "{}".format(c) + "}" + "=" + "{:>+7.4f}".format(yavg) + "$"

        if run == "s5m90":
            label =  "$D_{"  + "{}".format(c) + "}" + "$"
        else:
            label = ""


        iax.plot(x, y, 'o',  mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
        iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
        icolor = icolor +1

#-- load and plot correlarors vs m --

iax = ax[1]

dat = np.loadtxt("corr_alpha12_all.txt")

allm = np.arange(3,8)
sumDinv = np.zeros(5)
sumDdir = np.zeros(5)

for m in allm:
    ind = (dat[:,1] == m)
    sumDinv[m-3] = np.sqrt(np.sum( (dat[ind,4])**2  ))
    sumDdir[m-3] = np.sqrt(np.sum( (dat[ind,6])**2  ))
    

#iax.errorbar(dat[:,1], -dat[:,4], 2*dat[:,5],
#         marker='^',  mfc='none', ms=6, color=colors[3], lw=LW)

iax.plot(dat[:,1], -dat[:,4],
         'v',  mfc='none', ms=3, color=colors[2], lw=LW, label="inverse")
iax.plot(dat[:,1],  dat[:,6],
         '^',  mfc='none',  ms=3, color=colors[3], lw=LW, label="direct")

iax.plot(allm, sumDinv,
         'o',  mfc='none', ms=6, color=colors[2], lw=LW, label="$\sqrt{\sum D_j^2}$")
iax.plot(allm, sumDdir,
         's',  mfc='none', ms=6, color=colors[3], lw=LW, label="$\sqrt{\sum D_j^2}$")

i=np.arange(2.5,8)
#iax.plot(i, 10*np.exp(-i), color="k", lw=LW/2)
#iax.plot(i, 10*np.exp(-i), color="k", lw=LW/2)
iax.plot(i, 45*i**(-4), color="k", lw=LW/2, label="$\propto m^{-4}$")
#iax.plot(i, 10*i**(-3), color="k", lw=LW/2, label="$\propto m^{-3}$")
iax.set_yscale('log')
iax.set_xscale('log')


        
#-- canvas options --

iax = ax[0]
iax.set_ylim(-0.6,0.6)
iax.set(xlabel='$|i - d|$') # ylabel="$D_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
iax.set_xlim(0,nmodes)
iax.legend(frameon=False, handlelength=0.1, loc="upper left", bbox_to_anchor=(-0.47,1.05))
iax.set_yticks(np.arange(-0.6, 0.7, step=0.2))


iax = ax[1]
iax.set_xlim(2.5,7.5)
iax.set_ylim(1e-4,1)
iax.set(xlabel='$m$') # ylabel="$D_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
iax.legend(frameon=False,  loc="lower left",  handlelength=1,  bbox_to_anchor=(0.05,0.03))
iax.set_xticks(np.arange(3, 8, step=1))
#iax.ticklabel_format(axis='x', style="plain")
iax.set_xticklabels(["$3$", "$4$", "$5$", "$6$", "$7$"])

#-----------------------

plt.subplots_adjust(left=0.17, right=0.99, top=0.97, bottom=0.18, hspace=0.15, wspace=0.30)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
