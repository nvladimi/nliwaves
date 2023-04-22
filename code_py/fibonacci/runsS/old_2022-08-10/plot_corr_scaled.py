import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import pickle
  
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

correlators = "hi9"

outfmt = "{:>8.4f}{:>8.4f}"


#----------------------
if correlators == "hi1":

    outfile = "plot_corr_hi1sc.pdf"
    
    corr = ( ("0",    1),
             ("21",   1),
             ("322",  1),
             ("431",  1),
    )

    lblfmt = "{:>6.3f} i {:>+5.2f}"
    
    tiles =( (( 0,40),     ( 0,10),   (  0, 5)),
             ((-1.2,1.2),  (-1.2,1.2),  (-1.2,1.2)),
             ((0.5,2.5),    (-6,6),  (-6,6)),
             ((-0.5,1.5),      (-3,3),  (-3,3)))

    ymax=1
    
#------------------------    
elif correlators == "hi2":

    outfile = "plot_corr_hi2sc.pdf"
    
    corr = ( ("4332",    1),
             ("5422",    1),
             ("5441",    1),
             ("6531",    1),
    )

    lblfmt = "{:>6.3f} i {:>+5.2f}"

    
    tiles =( (( 1,5),  ( -20,20),   ( -20,20)),
             (( 1,5),  ( -40,20),   ( -40,20)),
             (( 1,4),  ( -15,15),   ( -15,15)),
             (( 1,4),  ( -20,10),   ( -20,10)))
    ymax=1
 

#------------------------ 
elif correlators == "hi3":

    outfile = "plot_corr_hi3sc.pdf"
     
    corr = (
        ("44333",   1 ),
        ("54432",   1 ),
        ("65332",   1 ),
        ("65522",   1 ),
    )

    lblfmt = "{:>5.2f} i {:>+6.1f}"
 
    
    tiles =(
        ((2,8),  ( -200,200), ( -200, 200)),
        ((2,8),  ( -250,250), ( -250, 250)),
        ((1,7),  (  -60, 60), (  -60, 60)),
        ((2,8),  ( -400,200), ( -400,200))
    )

    ymax=6

#------------------------ 
elif correlators == "hi4":

    outfile = "plot_corr_hi4sc.pdf"
    
    corr = ( 
        ("65541",   1 ),
        ("76422",   1 ),
        ("76441",   1 ),
        ("76631",   1 ),
    )

    lblfmt = "{:>5.2f} i {:>+6.1f}"
 
    
    tiles =(
        ((2,8),  ( -200,100), ( -200,100)),
        ((2,8),  ( -200,100), ( -200,100)),
        ((0,6),  ( -100,200), ( -100,200)),
        ((2,8),  ( -300,200), ( -300,200))      
    )
 
    ymax=6


#------------------------ 
elif correlators == "hi5":

    outfile = "plot_corr_hi5sc.pdf"
    
    corr = ( 
        ("87531",   1),
        ("63333",   1),
        ("5544443",   1),
        ("55544444",  1)
    )

    lblfmt = "{:>5.2f} i {:>+6.2n}"


    
    tiles =(
        ((1, 7),      ( -60,60), ( -60, 60)),
        ((2, 8),      ( -150,150), ( -150,150)), 
        ((6,14),  ( -30000,20000),   ( -30000, 20000)),
        ((8,18),( -4.e5,4.e5), ( -4e5, 4e5)) 
    )
    
    ymax = (6, 6, 100, 600)
    
#------------------------   
elif correlators == "hi6":

    outfile = "plot_corr_hi6sc.pdf"

    lblfmt = "{:>5.2f} i {:>+6.1f}"
    
    corr = (
        ("544433",    1),
        ("554442",    1),
        ("654333",    1),
        ("655432",    1),
    )

    tiles =(
        ((4,12),  ( -3000,2000),   ( -3000, 2000)),
        ((4,12),  ( -3000,2000),   ( -3000, 2000)),
        ((2,10),  (  -400,1200),   (  -400, 1200)),
        ((4,10),  ( -1500,1000),   ( -1500, 1000)),        
     )

    ymax = (15, 15, 15, 15)


#------------------------   
elif correlators == "hi7":

    outfile = "plot_corr_hi7sc.pdf"


    lblfmt = "{:>5.2f} i {:>+6.1f}"

    
    corr = (
        ("665551",    1),
        ("764432",    1),
        ("766332",    1),
        ("766522",    1),
    )

    tiles =(
        ((4,12),  ( -5000,1000),   ( -5000, 1000)),
        ((4,10),  ( -2000,1000),   ( -2000, 1000)),
        ((2,10),  ( -500,2000),    (  -500, 2000)),
        ((4,12),  ( -4000,2000),   ( -4000, 2000)),        
     )

    ymax = (15, 15, 15, 15)


#------------------------   
elif correlators == "hi8":

    outfile = "plot_corr_hi8sc.pdf"


    lblfmt = "{:>5.2f} i {:>+6.1f}"


    
    corr = (
        ("766541",    1),
        ("875332",    1),
        ("875522",    1),
        ("877422",    1),
    )

    tiles =(
        ((4,10),  ( -1000,500),   ( -1000, 500)),
        ((2,10),  (  -500,500),   (   -500, 500)),
        ((3,11),  ( -2000,1000),   ( -2000, 1000)),
        ((2,12),  ( -1000,2000),   ( -1000, 3000)),        
     )

    ymax = (15, 15, 15, 15)

 

#------------------------   
elif correlators == "hi9":

    outfile = "plot_corr_hi9sc.pdf"

    lblfmt = "{:>5.2f} i {:>+6.1f}"

    
    corr = (
        ("875541",    1),
        ("877441",    1),
        ("877631",    1),
        ("000000",    1),
    )

    tiles =(
        ((2,10),  ( -2000,1000),   ( -2000, 1000)),
        ((4,12),  ( -2000,6000),   ( -2000, 6000)),
        ((4,12),  ( -6000,2000),   ( -6000, 2000)),
        ((4,12),  ( -3000,3000),   ( -3000, 3000)),        
     )

    ymax = (15, 15, 15, 15)

             
#------------------------   
else:
    quit()
    
#---------------------------------------------------------------------------------





runs = ( ('s4m90', 1/2, 100,  90, 'or'),
         ('s5m90', 1/2, 200, 190, 'oy'),
         ('s4m10', 1/2, 100,  10, 'ob'),
         ('s5m10', 1/2, 200,  10, 'oc'),
         ('s3v38', 3/8,  60,  50, 'or'),
         ('s4v38', 3/8, 100,  90, 'oy'),
         ('s3v58', 5/8,  60,  10, 'ob'),
         ('s4v58', 5/8, 100,  10, 'oc'),
         ('s3v14', 1/4,  60,  50, 'or'),
         ('s4v14', 1/4, 100,  90, 'oy'),
         ('s3v34', 3/4,  60,  10, 'ob'),
         ('s4v34', 3/4, 100,  10, 'oc') )

#--------------------------------


fig, ax = plt.subplots(ncols=3, nrows=4, figsize=(10,12))

        
#-- common parameters --

MS = 1
LW = 1
#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, 


def cfactor(alpha,n):
    
    phi = (1 + np.sqrt(5))/2
    return  -phi**((1 + alpha)*n/3)


def leastsq(x0,y0):
    #Fit a line, y = ax + b, through some noisy data-points:

    #i1 = int( 0.2 * len(x0) )
    #i2 = int( 0.8 * len(x0) )
    i1 = 20
    i2 = len(x0)-20

    x = x0[i1:i2]
    y = y0[i1:i2]
    
    A = np.vstack([x, np.ones(len(x))]).T

    a, b = np.linalg.lstsq(A, y, rcond=None)[0]

    return a, b



#-- load and plot data --

for j in range(0,4):
    
    c = corr[j][0];
    n=0
    for cc in c:
        n += int(cc)
    for i in range(0,3):
        iax = ax[j][i]
        for r in range(i*4,i*4+4):
            try:
                run = runs[r]
                cf = cfactor(run[1],n)
                rname = "Post1j/" + run[0]  + "_" + c + "_mcr.txt"
                m=len(c)+1
                r=(m-3)/3
                #print(c, cf, rname)
                dat = np.loadtxt(rname)
                nk = dat[:,1]
                x = dat[:,0]
                y = dat[:,3] *cf
                y = y / nk**(0.5*m-1)
                
                iax.plot(x, y, run[4],  mfc='none', ms=MS,  lw=LW)
                """
                try:
                    a,b = leastsq(x,y)
                except:
                    a=0
                    b=0
                iax.plot(x, a*x+b, ':k', lw=LW/2, label=lblfmt.format(a,b))
                if (i==0):
                    q = y - r*x
                    d = np.average(q[20:-20])
                    iax.plot(x, r*x+d, '--k', lw=LW/2, label=lblfmt.format(r,d))
                else:
                    d=0    
                print(c, rname, outfmt.format(a,b),  outfmt.format(r,d) )
                #print(c, rname, outfmt.format(a,b))
                """
            except:
                pass


if correlators == "hi1":

    for i in range(0,3):
        iax = ax[0][i]
        for r in range(i*4,i*4+4):
            try:
                run = runs[r]
                rname = "Post1j/" + run[0]  + "_21_mcr.txt"
                dat = np.loadtxt(rname)
                x = dat[:,0];
                y = dat[:,1];               
                iax.plot(x, y, run[4],  mfc='none', ms=MS,  lw=LW)
            except:
                pass

#-- canvas options --


ax[0][0].set_title("$\\alpha = 1/2$")
ax[0][1].set_title("$\\alpha - 1/2 = \pm 1/8$")
ax[0][2].set_title("$\\alpha - 1/2 = \pm 1/4$")


for j in range(0,4):
    for i in range(0,3):
        ax[j][i].set_ylim(-ymax[j],ymax[j])
    ax[j][0].set(xlabel='$i$', ylabel="$C_{" + corr[j][0]  +"} \,  n^{1 - m/2}$")
    ax[j][0].set_xlim(0,200)
    ax[j][1].set_xlim(0,100)
    ax[j][2].set_xlim(0,100)


if (correlators == "hi1"):
    ax[0][0].set_xlim(0,200)
    ax[0][0].set(xlabel='$i$', ylabel="$C_0$")


#-----------------------

plt.subplots_adjust(left=0.09, right=0.97, top=0.95, bottom=0.05, hspace=0.25, wspace=0.28)

plt.savefig(outfile, pad_inches=0)

plt.close()

   



#===========================================



#=========================================================
