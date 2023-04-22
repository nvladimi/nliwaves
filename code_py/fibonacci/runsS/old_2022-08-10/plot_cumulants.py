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

outfile = "plot_cumulants.pdf"
#outfmt = "{:>8.4f}{:>8.4f}"
#lblfmt = "{:>6.3f} i {:>+5.2f}"


#---------------------------------------------------------------------------------

fig, ax = plt.subplots(ncols=2, figsize=(5.5,2.5))
#fig, ax = plt.subplots(ncols=2, figsize=(7,3))


#-- common parameters --

MS = 1
LW = 1

iax = ax[0]

#---------------------------------------------------------------------------------



nmodes=200

corrs = ("21-0", "322-0", "431-0", "11-30", "11-540", "221-20", "221-430", "3221-320")

#corrs = ("21-0", "221-20", "221-430", "3221-320")

param= {
    
    "21-0":     (1, "${21\\bar{0}}$", '#d00000', 50, -0.09),
    
    "322-0":    (2, "${322\\bar{0}}$", '#e07000', 50, 0.03), 
    "431-0":    (1, "${431\\bar{0}}$", '#e07000', 70, -0.08),
    "11-30":    (2, "${\\bar{3}11\\bar{0}}$", '#e07000', 100, 0.03),

    "11-540":   (2, "${\\bar{5}\\bar{4}11\\bar{0}}$", '#40a000', 140, -0.13),
    "221-20":   (6, "${22\\bar{2}1\\bar{0}}$", '#00a0f0', 120, 0.05),
    "221-430":  (2, "${\\bar{4}\\bar{3}221\\bar{0}}$", '#0000f0', 50, -0.09),
    "3221-320": (12,"${3\\bar{3}22\\bar{2}1\\bar{0}}$", '#c000f0', 20, 0.04),

}


#run = "s5m10"; d = nmodes
run = "s5m90"; d = 1

rname = "Links/" + run  + "_" + "21-0" + "_mcr.txt"
dat = np.loadtxt(rname)
y210 = dat[:,3]


for c in corrs:

    rname = "Links/" + run  + "_" + c + "_mcr.txt"
    ffactor = param[c][0] 
    label   = param[c][1] 
    color   = param[c][2] 
    xoffset = param[c][3] 
    yoffset = param[c][4] 

    dat = np.loadtxt(rname)
    nk = dat[:,1]
    x = dat[:,0]
    y = dat[:,3]

    if c == "221-20":
        #print("#221-20	 (21-0)+(2-2)x2")
        n2 = np.hstack(( nk[2:], (0,0) ))
        y = y - y210*n2*2
        
    if c == "221-430":
        #print("#221-430  (21-0)+(2-43)x2")
        yy = np.hstack(( y210[2:], (0,0) ))
        y = y - y210*yy*2
        
    if c == "3221-320":
        #print("#3221-320   (3-3)(2-2)(21-0)x2; (322-0)(1-32) small")
        n2 = np.hstack(( nk[2:], (0,0) ))
        n3 = np.hstack(( nk[3:], (0,0,0) ))
        y = y - y210*n2*n3*2

    y = y/dat[:,2] * np.abs(x-d) /  np.sqrt(ffactor)

    yavg = np.average(y[20:180])
    print(c, yavg)
    iax.text(xoffset, yavg + yoffset, label,color=color)
    #t.set_bbox(dict(facecolor='white', alpha=0.75,  edgecolor='white'))
    
    iax.plot(x, y, 'o',  mfc='none', ms=MS, color=color, lw=LW, label=label)
    iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=color)

iax.axhline(y=0, xmin=0, xmax=0.58, color='k', lw=LW/2, ls='--')
iax.axhline(y=0, xmin=0.78, xmax=1, color='k', lw=LW/2, ls='--')

    
#iax.axhline(y=0, color='k', lw=LW/2, ls='--')
    
#iax.set_ylim(-0.,0.6)
iax.set_ylim(-0.5,0.5)
iax.set(xlabel='$|i - d|$')
iax.set(ylabel='$|i - d| D_m$')
iax.yaxis.set_label_coords(-0.2,0.5)

iax.set_xlim(0,nmodes)
#iax.legend(frameon=False, handlelength=-0.5, loc="upper left", bbox_to_anchor=(1.01,1.05))
iax.set_yticks(np.arange(-0.4, 0.5, step=0.2))

        
#------------------------------------------------------------------------------------
iax = ax[1]
iax.set_ylim(-0.3,0.3)
iax.set_xlim(0,nmodes)
iax.set(xlabel='$|i - d|$') 
iax.text(25,0, "DOUBLING MODEL DATA")


#------------------------------------------------------------------------------------




plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.19, wspace=0.25)

plt.savefig(outfile, pad_inches=0)

plt.close()

   
#------------------------------------------------------------------------------------

