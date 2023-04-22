import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import scipy.io
import math



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


#-- common parameters --

MS = 1
LW = 1
clr = ('#f00000', '#40a000', '#0000f0', '#00a0f0', '#c000f0')

iax = ax[0]



#---------------------------------------------------------------------------------



nmodes=200

corrs = ("21-0", "11-30", "11-540", "221-430", "3221-320")
#corrs = ("21-0", "322-0", "431-0", "11-30", "11-540", "221-20", "221-430", "3221-320")
#corrs = ("21-0", "221-20", "221-430", "3221-320")

param= {
    
    "21-0":     (1, "${21\\bar{0}}$",                   clr[0],  50, -0.09),
    
    "322-0":    (2, "${322\\bar{0}}$",                  clr[1],  50,  0.03),
    "431-0":    (1, "${431\\bar{0}}$",                  clr[1],  80, -0.09),
    "11-30":    (2, "${\\bar{3}11\\bar{0}}$",           clr[1], 100,  0.03),

    "11-540":   (2, "${\\bar{5}\\bar{4}11\\bar{0}}$",   clr[2], 140, -0.13),
    "221-20":   (6, "${22\\bar{2}1\\bar{0}}$",          clr[2], 120,  0.05),

    "221-430":  (2, "${\\bar{4}\\bar{3}221\\bar{0}}$",  clr[3],  50, -0.09),
    "3221-320": (12,"${3\\bar{3}22\\bar{2}1\\bar{0}}$", clr[4],  20,  0.02),

}


#run = "s5m10"; d = nmodes
run = "s5m90"; d = 1
x=np.arange(1,201)
xd = np.abs(x-d)
prefactor=1.13
xd13 = (prefactor*xd)**(1/3)

rname = "Post05/" + run  + "_" + "21-0" + "_mcr.txt"
#rname = "Links/" + run  + "_" + "21-0" + "_mcr.txt"
dat = np.loadtxt(rname)
y210 = dat[:,3]

        
for c in corrs:

    rname = "Post05/" + run  + "_" + c + "_mcr.txt"
    ffactor = param[c][0] 
    label   = param[c][1] 
    color   = param[c][2] 
    xoffset = param[c][3] 
    yoffset = param[c][4]
    m=len(c)-1

#------------

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
        print(cc, nu)

    combfactor = np.sqrt(combfactor)

    dat = np.loadtxt(rname)
    nk = dat[:,1]
    y = dat[:,3]

    #denom = dat[:,2]*combfactor;  iax.set_title("old")  #-- old
    #denom = dat[:,5]; iax.set_title("new")              #-- new
    denom = denom*combfactor; #iax.set_title("gauss")     #-- gauss
   
 
#-------------
 
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

    y = y/denom * xd

    
    yavg = np.average(y[20:180])
    print(c, yavg)
    iax.text(xoffset, yavg + yoffset, label,color=color)
   
    iax.plot(x, y, 'o',  mfc='none', ms=MS, color=color, lw=LW, label=label)
    #iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=color)

#iax.axhline(y=0, xmin=0, xmax=0.58, color='k', lw=LW/2, ls='--')
#iax.axhline(y=0, xmin=0.78, xmax=1, color='k', lw=LW/2, ls='--')
iax.axhline(y=0, color='k', lw=LW/2, ls='--')



iax.set_ylim(-0.5,0.5)
iax.set(xlabel='$|j - d|$')
iax.set(ylabel='$|j - d| D_m$')
iax.yaxis.set_label_coords(-0.2,0.5)
iax.set_xlim(0,nmodes)
iax.set_yticks(np.arange(-0.4, 0.5, step=0.2))
   
#------------------------------------------------------------------------------------
iax = ax[1]
iax.set_ylim(-0.4,0.4)
iax.set_xlim(0,200)
iax.set(xlabel='$|j - d|$') 
#iax.set_yticks(np.arange(-0, 0.16, step=0.05))
iax.axhline(y=0, xmin=0, xmax=0.58, color='k', lw=LW/2, ls='--')
iax.axhline(y=0, xmin=0.78, xmax=1, color='k', lw=LW/2, ls='--')


#------------------------------------------------------------------------------------

direct=False

#{'JJ','        C11',     'JJc1',         'C221',   'C221J',       'n1J'};
#{'$11110`0`$', '$110`$', '$2`2`1110`$', '$2210`$', '$221110`0`$', '$1111`0`$'};


if (direct == True):
    DCorr = scipy.io.loadmat('Irreducible_cumulants_doubling/NormD.mat')
    x = np.arange(200,0,-1)-1
else:
    DCorr = scipy.io.loadmat('Irreducible_cumulants_doubling/NormI.mat')
    x = np.arange(0,200)


Dparam= {
    "C11":    ("${11\\bar{0}}$",                  clr[0], 120,  0.010),
    "C221":   ("${221\\bar{0}}$",                 clr[1], 120,  0.020),
    "n1J":    ("${111\\bar{1}\\bar{0}}$",         clr[2], 120, -0.070),
    "JJ":     ("${1111\\bar{0}\\bar{0}}$",        clr[3], 120,  0.020),
    "JJc1":   ("${\\bar{2}\\bar{2}111\\bar{0}}$", clr[3],   0,  0),
    "C221J":  ("${22111\\bar{0}\\bar{0}}$",       clr[4], 120,  0.020),
}




corr = ('C11', 'C221', 'n1J', 'JJ', 'C221J', 'n1J')

for c in corr:
    y = DCorr['y'+c].reshape(200)
    #print(y)
    
    label   = Dparam[c][0] 
    color   = Dparam[c][1] 
    xoffset = Dparam[c][2] 
    yoffset = Dparam[c][3] 
    y = y*x

    yavg = np.average(y[20:180])
    print(c, yavg)
    iax.text(xoffset, yavg + yoffset, label,color=color)
    
    iax.plot(x, y, 'o',  mfc='none', ms=MS, color=color, lw=LW, label=label)
    #iax.plot(x, x*0 + yavg, '-',  lw=LW/2, color=color)


plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.19, wspace=0.25)

plt.savefig(outfile, pad_inches=0)

plt.close()

   
#------------------------------------------------------------------------------------

