import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import math

#======================================================================================

def all(filename, rdir="Post04", run="s5m10", method="gauss"):

    
    with open(filename) as file:
        lines = file.readlines()


    for line in lines:

        cc=line.rstrip()
        
        if cc == "":
            print(line)
            continue
        if cc[0] == "#":
            print(line)
            continue
        
        cc = cc.split()

        c0 = cc[0]

        if len(cc)==1:
            mult = 0
        else:
            mult = int(cc[-1][1])

        try:
            cFit(c0, cc[1:-1], mult=mult, rdir=rdir, run=run, method=method)
        except:
            print("Error at:  " + line)



#======================================================================================

def cFit(c0, csub, mult=1, rdir="Post05", run="s5m90", method='gauss', toff=0):

# c0       main correlator 
# csub     list of correlators in subtracted term
# mult     combinatorial factor for subtracted term
# toff     seconds the plot is display; prints to file for toff > 60

    nmodes=200; i1=20; i2=180;
    
    #-- main correlator --
    
    rname = rdir + '/' + run  + "_" + c0 + "_mcr.txt"
    dat0 = np.loadtxt(rname)
    
    if run in ('s5m10', 's4v34',  's4v58'):
        d = nmodes
    else:
        d = 1
        
    x  = np.arange(1,nmodes+1)
    xd = np.abs(x-d)

    denom, combfactor = cDenom(c0, xd, nmodes)

    if method == 'old':
        denom = dat0[:,2]*combfactor
    elif method == 'new':
        denom = dat0[:,5]
    else:
        denom = denom*combfactor

    #-- subtractions --
    
    x  = x[i1:i2] 
    xd = xd[i1:i2] 
    denom = denom[i1:i2] 
    corr   = dat0[i1:i2, 3]

    y1 = xd/denom * corr
    y2 = xd/denom * mult

    
    for ck in csub:
    
        c2,ishift2 = cShift(ck);    # shift for subtracted correlator

        rname = rdir + '/' + run  + "_" + c2 + "_mcr.txt"
        dat2 = np.loadtxt(rname)
        csub2  = dat2[i1-ishift2 : i2-ishift2, 3]
        y2 = y2*csub2

        
    #-- average and plot --
        
    y  = y1-y2

    y1avg = np.average(y1)
    y2avg = np.average(y2)
    yavg  = np.average(y)
    yrms  = np.sqrt(np.average( (y - yavg)**2 ) )

    c = '  ('
    for i in csub:
        c = c + i + ")+("
    c= c[:-2] + ' x' + str(mult) + "  " 
        
    print(c0 + c + "({:>+7.4f}) - ({:>+7.4f}) = {:>+7.4f} +/-{:>7.4f} ".format(y1avg, y2avg, yavg, yrms)) 

    cPlot(c0, x, y, y1, y2, yavg, toff)


#======================================================================================

def cDenom(c0, xd, nmodes):

    xd13 = (1.13*xd)**(1/3)

    combfactor = 1
    denom = np.ones(nmodes)
    for cc in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"):   # "A"):
        nu = c0.count(cc)
        ishift = int(cc)
        if ishift > 0:
            denom = denom * np.hstack( (np.zeros(ishift), xd13[0:-ishift]**nu) )
        else:
            denom = denom * xd13**nu
        combfactor = combfactor * math.factorial(nu)
        #print(cc, nu)
    combfactor = np.sqrt(combfactor)

    return(denom, combfactor)

#======================================================================================

def cPlot(c0, x, y, y1, y2, yavg, toff):
    
    outfile = "fit_corr_tmp.pdf"

    
    if toff > 0:

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5,4))

        MS = 1
        LW = 1

        #label =  c0 + "{:>+7.4f}".format(yavg) + "$"

        ax.plot(x, y, 'or',   mfc='none', ms=MS, lw=LW, label="Diff")
        ax.plot(x, y1, ':g',  mfc='none', ms=MS, lw=LW, label="Corr")
        ax.plot(x, y2, ':b',  mfc='none', ms=MS, lw=LW, label="Csub")
        ax.plot(x, x*0 + yavg, 'r-',  lw=LW/2)
        ax.set_title(c0)
        ax.set_xlabel("|i-d|")
        ax.axhline(y=0, color='k', lw=LW/2, ls='--')
        ax.legend(frameon=False)

        #plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.09, hspace=0.15, wspace=0.15)

        if toff > 60:
            plt.savefig(outfile, pad_inches=0)
        else:            
            plt.ion()
            plt.show()
            plt.pause(toff)
            plt.ioff()

        plt.close()

#======================================================================================

def cShift(corr):
    
    cc=corr.split("-") 

    x = int( min( min(cc[0]), min(cc[1]) ))

    c1=''; c2=''
    for c in cc[0]:
        c1 += str(int(c)-x)
    for c in cc[1]:
        c2 += str(int(c)-x)
        
    c1=''.join(sorted(c1))[::-1]
    c2=''.join(sorted(c2))[::-1]
    
    if c2[-1]=="0":
        c = c1 + '-' +c2
    else:
        c = c2 + '-' +c1

    return(c, x)

#======================================================================================

