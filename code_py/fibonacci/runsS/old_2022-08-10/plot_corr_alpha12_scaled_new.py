import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

matplotlib.rc('xtick', labelsize=11) 
matplotlib.rc('ytick', labelsize=11) 
matplotlib.rc('axes', labelsize=12)

output=1

outfile = "plot_corr_alpha12_new_inverse.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"

runs = ("s5m10", )

#runs = ("s5m90", "s5m10")
#runs = ("s4v14", "s4v34")
#runs = ("s4v38", "s4v58")







ymax = ((0.6, 0.06),  (0.15, 0.08),  (0.05, 0.05) ); nmodes=200
#ymax = ((3, 0.6),  (0.6, 0.6) ); nmodes = 100



corrs = ( ( ("21-0", "322-0", "431-0", "11-30"),
            ("4332-0",  "5422-0", "5441-0", "6531-0", "411-20", "331-50", '222-40', '11-540'),),        
          ( ("3332-50", "4422-60",  "5421-30", "4441-60", "5531-70", "6511-20"),
            ("6331-40", "5222-30", "4111-00", "8721-60", "7621-50", "6521-40", "2211-00"),),
          ( ("611-440", "311-520", "331-760", "222-650", "621-870", "521-760", "421-650"),
            ( "321-540", "221-430","11-7640",  "11-6550",
              "210-00", "3222-20", "4310-00", '3220-00','110-300') )
)
#          ( ("44333", "54432", "65332", "65522", "76422", "76631"),        
#            ("544433", "554442",  "655432", "665551", "764432", "766522") ),
#          ( ("63333", "65541",  "76441",  "87531"), 
#            ("654333", "766332", "766541", "875332", "875522",
#             "877422",  "875541",  "877441",    "877631") )


ffactor= {
    
    "21-0":   1,
    
    "322-0":  2,
    "431-0":  1,
    "11-30":  2,
    
    "4332-0": 2,
    "5422-0": 2,
    "5441-0": 2,
    "6531-0": 1,

    "11-540": 2,
    "411-20": 2,
    "331-50": 2,
    "222-40": 6,
    
    "210-00": 6,
    
    "3332-50": 6,
    "4422-60": 4,
    "5421-30": 1,
    "4441-60": 6,
    "5531-70": 2,
    
    "6511-20": 2,
    "6331-40": 2,
    "5531-70": 2,
    "5222-30": 6,
    "5531-70": 2,

    "6511-20": 2,
    "6331-40": 2,
    "5222-30": 8,
    "4111-00": 8,
    "8721-60": 1,

    "7621-50": 1,
    "6521-40": 1,
    "2211-00": 8,
    "611-440": 4, 
    "311-520": 2,
    
    "331-760": 2,
    "222-650": 6,
    "621-870": 1,
    "521-760": 1,
    "421-650": 1,

    "321-540": 1,
    "221-430": 2,
    "11-7640": 2,
    "11-6550": 4,

    "3222-20": 24,
    "3220-00": 12,
    "4310-00": 6,
    "110-300": 12,
    
    "44333-0": 12,
    "54432-0": 2,
    "63333-0": 24,
    "65332-0": 2,
    "65522-0": 4,
    "65541-0": 2,

    "76422-0": 2,
    "76441-0": 2,
    "76631-0": 2,
    "87531-0": 1,

    "544433-0": 12,
    "554442-0": 12,
    "654333-0": 6,
    "655432-0": 2,
    "665551-0": 12,

    "764432-0": 2,
    "766332-0": 4,
    "766522-0": 4,
    "766541-0": 2,

    "875332-0": 2,
    "875522-0": 4,
    "877422-0": 4,
    "875541-0": 2,
    "877441-0": 4,
    "877631-0": 2,

    "5544443-0": 48,
    "55544444-0": 720
}



leglbl = (("$m = 3,4$",  "$m=5$"),  ("$m=5$", "$m=6$"),  ("$m=6$", "$m=6$") )




#---------------------------------------------------------------------------------


#fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.0,4.5))
fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(10,11))

        
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



#-- load and plot data --

for run in runs:

    if run in ( 's5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
    
    
    for j in range(0,3):
        for i in range(0,2):
            iax = ax[j][i]
            icolor = 0;
            for c in corrs[j][i]:

                try:
                    rname = "Post02/" + run  + "_" + c + "_mcr.txt"
                    m=len(c)-1
                    r=(m-3)/3
                    #print(rname)
                
                
                    dat = np.loadtxt(rname)
                    nk = dat[:,1]
                    x = dat[:,0]
                    y = dat[:,3]/dat[:,2] * np.abs(x-d) / np.sqrt(ffactor[c])

                    yavg = np.average(y[20:180])
                    yrms = np.sqrt(np.average( (y[20:180] - yavg)**2 ) )

                    #print(c, m, yavg, yrms)  #cfactor(c), ffactor[c],
                    print("{}  m={}:  {:>+7.4f} +/-  {:>+7.4f} ".format(c, m, yavg, yrms)) 
                          
                    label =  "$D_{"  + "{}".format(c) + "}" + "=" + "{:>+7.4f}".format(yavg) + "$"
                    #label =  "$D_{"  + "{}".format(c) + "}" + "$"

                    iax.plot(x, y, 'o',  mfc='none', ms=MS, color=colors[icolor], lw=LW, label=label)
                    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=colors[icolor])
                    icolor = icolor + 1
                except:
                    print(rname)

#-- canvas options --


for j in range(0,3):
    for i in range(0,2):
        ax[j][i].set_ylim(-ymax[j][i],ymax[j][i])
        ax[1][i].set(xlabel='$i$') # ylabel="$D_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
        ax[j][i].set_xlim(0,nmodes)
        ax[j][i].legend(frameon=False, handlelength=0.1, loc="upper left", bbox_to_anchor=(1.00,0.95))
        ax[j][i].text(1.05*nmodes, 0.9*ymax[j][i], leglbl[j][i], fontsize=12) 

#-----------------------

plt.subplots_adjust(left=0.06, right=0.85, top=0.99, bottom=0.09, hspace=0.15, wspace=0.66)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
