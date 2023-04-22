import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import scipy.io
import math
import fiboPost


#output=1

outfile = "plot_cumulants_retest.pdf"


#---------------------------------------------------------------------------------

fig, ax = plt.subplots(ncols=2, figsize=(5.5,2.5))


#-- common parameters --

MS = 1
LW = 1
clr = ('#f00000', '#40a000', '#0000f0', '#00a0f0', '#c000f0')

iax = ax[0]

#---------------------------------------------------------------------------------


nmodes=200

corrs = ("21-0", "322-0", "431-0", "11-30")

param= {
    
    "21-0":     (1, "${21\\bar{0}}$",                   clr[0],  50, -0.09),
    
    "322-0":    (2, "${322\\bar{0}}$",                  clr[1],  50,  0.03),
    "431-0":    (1, "${431\\bar{0}}$",                  clr[2],  80, -0.09),
    "11-30":    (2, "${\\bar{3}11\\bar{0}}$",           clr[3], 100,  0.03),

    "11-540":   (2, "${\\bar{5}\\bar{4}11\\bar{0}}$",   clr[2], 140, -0.13),
}


#run = "s5m10"; d = nmodes
run = "s5m90"; d = 1



#-------------------------

x=np.arange(1,nmodes+1)
xd = np.abs(x-d)
prefactor=1.13
xd13 = (prefactor*xd)**(1/3)


for c in corrs:

    Cm = fiboPost.Cmodel(c, x, d)
    
    rname = "../runsS/Post05/" + run  + "_" + c + "_mcr.txt"
    ffactor = param[c][0] 
    label   = param[c][1] 
    color   = param[c][2] 
    xoffset = param[c][3] 
    yoffset = param[c][4]
    m=len(c)-1

    combfactor = 1
    denom = np.ones(200)
    for cc in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"):   # "A"):
        nu = c.count(cc)
        ishift = int(cc)
        if ishift > 0:
            denom = denom * np.hstack( (np.zeros(ishift), xd13[0:-ishift]**nu) )
        else:
            denom = denom * xd13**nu
        combfactor = combfactor * math.factorial(nu)
        #print(cc, nu)

    combfactor = np.sqrt(combfactor)
    denom = denom*combfactor

    dat = np.loadtxt(rname)
    y = dat[:,3]
      
    y = y/denom * xd
    ym = Cm/denom *xd

    ax[0].plot(x, y, 'o',  mfc='none', ms=MS, color=color, lw=LW, label=label)
    ax[0].plot(x, ym, color='k', lw=LW/2, zorder=10)
    
    yavg = np.average(y[20:180])
    print(c, yavg)
    iax.text(xoffset, yavg + yoffset, label,color=color)
    iax.plot(x, y, 'o',  mfc='none', ms=MS, color=color, lw=LW, label=label)

     
iax.axhline(y=0, color='k', lw=LW/2, ls='--')

iax.set_ylim(-0.5,0.5)
iax.set(xlabel='$|j - d|$')
iax.set(ylabel='$|j - d| D_m$')
iax.yaxis.set_label_coords(-0.2,0.5)
iax.set_xlim(0,nmodes)
iax.set_yticks(np.arange(-0.4, 0.5, step=0.2))
   
#------------------------------------------------------------------------------------


plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.19, wspace=0.25)

plt.savefig(outfile, pad_inches=0)

plt.close()

   
#------------------------------------------------------------------------------------

