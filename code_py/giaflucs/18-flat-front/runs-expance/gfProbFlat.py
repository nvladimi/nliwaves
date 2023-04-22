#import importlib ; import gfProbFlat
#importlib.reload(gfProbFlat) ; gfProbFlat.Test(iz=1)

def Test(iz=1, cfold=10):
   
    fbase = "test/flat"

    fbaseout = fbase+'_iz' + str(iz).zfill(4)

    Iranges = ( (0.125, 0.5), (0.25, 1.0), (0.5, 4), (1,10), (100,520) )
    
    ProbFromCount(fbase, iz, range(20), dI=0.125, Imax=512.0,
                  dlnI=0.01, lnImin=-10, lnImax=6,
                  ranges=Iranges, fbaseout=fbaseout)

    #ProbFromCount(fbase, iz, nseeds=20, dI=0.125, Imax=512.0,
    #              dlnI=0.01, lnImin=-10, lnImax=6,
    #              cfold=cfold, fbaseout=fbaseout)

    
    #PlotProb(fbaseout)

#=================================================================================================

def ProbFromCount(fbase, iz, seeds,
                  dI, Imax, dlnI, lnImin, lnImax,
                  cfold=0, ranges=(), fbaseout=''):
    
    import numpy as np
 
    I = (np.arange(0, np.fix(Imax/dI)) + 0.5)*dI 
    cntI = np.zeros(np.shape(I), dtype=int)

    lnI = (np.arange(np.fix(lnImin/dlnI), np.fix(lnImax/dlnI)) + 0.5)*dlnI 
    cntlnI = np.zeros(np.shape(lnI)) 

    #-- collect total count -- 
    
    nseeds = len(seeds)
    for iseed in seeds:

        fname = fbase + '.' + str(iseed).zfill(6) + '_I_' +  str(iz).zfill(4)
        
        c = np.fromfile(fname, 'uint32')
        
        cntI += c

        fname = fbase + '.' + str(iseed).zfill(6) + '_lnI_' +  str(iz).zfill(4)
        
        c = np.fromfile(fname, 'uint32')

        cntlnI += c

        
    #-- use distribution to compute sigmas and averages --

    ntotI    = np.sum(cntI)
    pI       = cntI/ntotI/dI    

    Iavg     = np.sum(I * pI)*dI
    Isigma   = np.sqrt( np.sum((I - Iavg)**2 * pI) *dI )


    ntotlnI  = np.sum(cntlnI)
    plnI     = cntlnI/ntotlnI/dlnI    

    lnIavg   = np.sum(lnI * plnI)*dlnI
    lnIsigma = np.sqrt(np.sum((lnI - lnIavg)**2 * plnI)*dlnI)

    if cfold > 0:
        I, pI, dI = foldTail(I, cntI, cfold)
    elif len(ranges) > 0 : 
        I, pI, dI = reBin(I, cntI, ranges)
    else:
        I, pI, dI = trimTail(I, cntI)
        
            
    #--  save to files --

    if (len(fbaseout)>0):

        xf = (lnI-lnIavg)/lnIsigma
        yf = plnI*lnIsigma
        dataout = np.vstack((xf, yf)).transpose()

        fnameout = fbaseout + "_problogI.dat"
        
        h = "Run \"" + fname + "\",  screen " + str(iz) + ",  " + str(nseeds) + " seeds:  "
        h = h +  str(ntotlnI) + " elements\n"
        h = h + "# with <ln I> = " + str(lnIavg) + " +/-  " +  str(lnIsigma) + "\n\n"
        h = h + "1.(lnI-lnIavg)/lnIsigma   2.PDF  \n"
        
        np.savetxt(fnameout, dataout, header = h)


        dataout = np.vstack((I, pI)).transpose()

    
        fnameout = fbaseout + "_probI.dat"
    
        h = "Run \"" + fbase + "\":  " + str(ntotI) + " elements"
        h = h + " with <I> = " + str(Iavg) +  " +/-  " +  str(Isigma) + "\n\n"
        h = h + "1.I   2.PDF \n"
        
        np.savetxt(fnameout, dataout, header = h)

        
#=================================================================================================

def foldTail(I, cI, cmin):
    
    import numpy as np

    dI  = I[1] - I[0]
    Is  = I
    dIs = np.ones(I.shape)*dI
    cIs = cI
    
    Iavg = np.sum(I * cI) /  np.sum(cI) 
    
    #print(sum(dIs), sum(cIs))

    iavg = np.nonzero(Is>Iavg)[0][0]
    
    while (np.min(cIs[iavg:-1]) < cmin):
 
        i0 = np.nonzero((cIs<cmin) & (Is>Iavg))[0][0]
        n0 = Is.shape[0]
        i1 = np.arange(i0,n0,2)
        i2 = i1+1
        
        if (i2[-1] == n0-1):
            Inew  = np.hstack((Is[:i0], (Is[i1] + Is[i2])*0.5))
            cInew = np.hstack((cIs[:i0], cIs[i1] + cIs[i2]))
            dInew = np.hstack((dIs[:i0], dIs[i1] + dIs[i2]))
        else:
            Inew  = np.hstack((Is[:i0], (Is[i1[:-1]] + Is[i2[:-1]])*0.5,  Is[-1]))
            cInew = np.hstack((cIs[:i0], cIs[i1[:-1]] + cIs[i2[:-1]], cIs[-1]))
            dInew = np.hstack((dIs[:i0], dIs[i1[:-1]] + dIs[i2[:-1]], dIs[-1]))
    
        Is = Inew
        cIs = cInew
        dIs = dInew

        #print(sum(dIs), sum(cIs), i0,  n0, Is[i0], cI[i0-1:i0+3],  cI[-3:])
        
        #print(' I:',  Is[i0-1:i0+3],  Is[-3:])
        #print('cI:', cIs[i0-1:i0+3], cIs[-3:])
        #print('dI:', dIs[i0-1:i0+3], dIs[-3:])
        
    pIs = cIs/dIs /sum(cIs) 
        
    return Is, pIs, dIs

        
#=================================================================================================

def trimTail(I, cI):
    
    import numpy as np
    
    dIo  = I[1] - I[0]

    q = np.nonzero(cI>0)[0][-1] + 1
    cIs = cI[:q]
    Is  = I[:q]
    dIs = dIo * np.ones(np.shape(Is))
    pIs = cIs / np.sum(cIs) / dIo
  
    return Is, pIs, dIs


#=================================================================================================

def reBin(I, cI, ranges, debug=False):
    
    import numpy as np

    dIo  = I[1] - I[0]

    cIs = list({})
    dIs = list({})
    Is  = list({})
    
    i1 = 0
    
    for r in ranges:
        
        i2 = r[1]/dIo
        n = max(0, int(r[0]/dIo))
        
        while (i1 < i2):
            
            cIr = cI[i1:i1+n]
            c = np.sum(cIr)

            nr = min(n, len(cIr))
            Ir = (np.arange(nr)+0.5)*dIo
            s1 = np.sum(cIr)*dIo
            s2 = np.sum(cIr*Ir)*dIo
            width = dIo * nr
            q = 2*s2/s1/width
           
            #center = i1*dIo + 0.5 * width                   # center of box
            #center = i1*dIo + 0.5 * width * q               # center of mass
            #center = i1*dIo + 0.5 * width * (2-q)/(4-3*q)   # moments
            center = i1*dIo + 0.5 * width * 0.5*(q+1)        # empirical
                
    
            cIs.append(c)    
            Is.append(center)
            dIs.append(width)

            if debug:
                print("first bin = {}, bins = {}, interval = [{}, {}], q = {}, center = {}, count {}".format(
                    i1, n, i1*dIo, (i1+n)*dIo, q, center, c))
            i1 += n
                  
    cIs = np.asarray(cIs)
    Is  = np.asarray(Is)
    dIs = np.asarray(dIs)

    q = np.nonzero(cIs>0)[0][-1] + 1
    cIs = cIs[:q]
    Is  = Is[:q]
    dIs = dIs[:q]
    
    pIs = cIs / np.sum(cIs) / dIs
    pI  = cI  / np.sum(cI)  / dIo


    if debug:
    
        Iavg  = np.sum(I  * pI * dIo) 
        Isavg = np.sum(Is * pIs * dIs) 
    
        print("Totals: ", np.sum(cI), np.sum(cIs)),  
 
        print("Averages: ", Iavg,  Isavg) 

        #-- show plot --
        
        import matplotlib.pyplot as plt
        plt.ion()

        fig, ax1 = plt.subplots(1, 1, figsize=(6,4))

        ax1.grid()
 
        ax1.set_yscale('log')

        q = np.nonzero(cI>0)[0][-1] + 1

        ax1.set_xlim(0,q*dIo)
        ax1.set(xlabel="$I$", ylabel="PDF")
    
        ax1.plot(I, cI/np.sum(cI)/dIo, 'ob',  markersize=2, linewidth=1,  mfc='none')
        ax1.plot(Is, pIs, 'or',  markersize=4, linewidth=1,  mfc='none')
        
        plt.tight_layout()
        
        plt.show()
    
    return Is, pIs, dIs

    
    
        
#=================================================================================================

def ProbFromField(fbase, N, iz, nseeds=1,
                  binsizeI = 0.5, numbinsI = 100,
                  binsize = 0.2, numbins = 100, fbaseout=''):
    
    import numpy as np
   
    binedges = np.arange(-numbins/2, numbins/2+1)*binsize
    binedgesI = np.arange(numbinsI+1)*binsizeI

    fsum = 0
    Isum = 0
    
    #-- compute average -- 
    
    for iseed in range(0,nseeds):

        fname = fbase + '.' + str(iseed).zfill(6) + '_' +  str(iz).zfill(4)
        
        a = np.fromfile(fname, 'float64')
        I = a**2
        f = np.log(I)

        fsum += np.average(f)
        Isum += np.average(I)
        
    favg = fsum/nseeds
    Iavg = Isum/nseeds

    #-- compute distributions and sigmas --

    p  = np.zeros(numbins)
    pI = np.zeros(numbinsI)

    fsum = 0
    Isum = 0
      
    for iseed in range(1,nseeds+1):

        fname = fbase + '.' + str(iseed).zfill(6) + '_' +  str(iz).zfill(4)
        
        a = np.fromfile(fname, 'float64')
        I = a**2
        f = 2*np.log(a)

        fsum += sum((f - favg)**2)
        Isum += sum((I - Iavg)**2)
        
        (p0, x) = np.histogram(f, binedges, density=False)
        p += p0

        (p0, x) = np.histogram(I, binedgesI, density=False)
        pI += p0
                 
    sigma  = np.sqrt(fsum / (N*N*nseeds))
    sigmaI = np.sqrt(Isum / (N*N*nseeds))
    
    x  = 0.5 * (binedges[:-1] + binedges[1:])
    x = (x - favg)/sigma
    p = p / (N*N*nseeds) / binsize * sigma
    
    xI  = 0.5 * (binedgesI[:-1] + binedgesI[1:])
    pI = pI / (N*N*nseeds) / binsizeI


    #--  save to files --

    if (len(fbaseout)>0):
    
        dataout = np.vstack((x, p)).transpose()

        fnameout = fbaseout + "_problogI.dat"
        
        h = "Run \"" + fbase + "\":  " + str(N*N*nseeds/1e6) + "M elements"
        h = h + " with <ln I> = " + str(favg) + " +/-  " +  str(sigmaI) + "\n\n"
        h = h + "1.(lnI-lnIavg)/sigma   2.PDF  \n"
        
        np.savetxt(fnameout, dataout, header = h)


        dataout = np.vstack((xI, pI)).transpose()

    
        fnameout = fbaseout + "_probI.dat"
    
        h = "Run \"" + fbase + "\":  " + str(N*N*nseeds/1e6) + "M elements"
        h = h + " with <I> = " + str(Iavg) +  " +/-  " +  str(sigmaI) + "\n\n"
        h = h + "1.I   2.PDF \n"
        
        np.savetxt(fnameout, dataout, header = h)

#=================================================================================================


def PlotProb(fbase):
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    
    dat = np.loadtxt(fbase+"_probI.dat")
    xI = dat[:,0]
    pI = dat[:,1]

    dat = np.loadtxt(fbase+"_problogI.dat")
    xL = dat[:,0]
    pL = dat[:,1]

    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    ax1.grid()
 
    ax1.set_yscale('log')
    ax1.set_xlim(0,10)
    ax1.set_ylim(1e-8,10)
    ax1.set(xlabel="$I$", ylabel="PDF")
    
    ax1.plot(xI, pI, 'or',  markersize=4, linewidth=1,  mfc='none', label='$z = 5.31$')

    a = 3.1; q = 1; lbl= '$\exp(-3.1 \, I)$'
    #a = 1.8; q = 0.62; lbl = '$\exp(-1.8 \, I^{0.62})$'
    #a = 1.8; q = 0.67; lbl = '$\exp(-1.8 \, I^{0.67})$'
    
    dx = 0.1
    x  = np.arange(0, 40, dx)
    y  = np.exp( -a*x**q)
    y  = y / sum(y)/dx
    
    
    ax1.plot(x, y, '-r', linewidth=1, label=lbl)
    
    ax1.legend(frameon=False)
        
    ax2.set_yscale('log')
    ax2.set_xlim(-5,5)
    ax2.set_ylim(1e-6,1)
    ax2.set_yscale('log')
    ax2.plot(xL, pL, 'o--r',  markersize=4, linewidth=1,  mfc='none')


    ax2.set(xlabel="$(\ln I  - \langle \ln I \\rangle ) /\sigma$", ylabel="PDF")

    
    plt.tight_layout()
    plt.savefig(fbase+'.pdf', pad_inches=0)

    # plt.close()
        
    plt.show()


#=================================================================================================


def probBinEqual(f,  numbins=2000,  bincoef=0.01):

    import numpy as np

    favg = np.average(f)

    binsize = bincoef*favg    
    binedges = np.arange(numbins)*binsize
    
    (p, x) = np.histogram(f, binedges, density=True)
    
    x  = 0.5 * (binedges[:-1] + binedges[1:])
    dx = np.vstack((x-binedges[:-1],  binedges[1:]-x))

    return(p, x, dx)


def probBinStretch(f,  numbins=50,  coef0=1e-2, coeflast=1e2):

    import numpy as np

    favg = np.average(f)

    binsize0 = coef0*favg
    lastedge = coeflast*favg
    c = np.log(1+ lastedge/binsize0)/numbins
    binedges = (np.exp(c * np.arange(numbins+1)) - 1)*binsize0
    
    (p, x) = np.histogram(f, binedges, density=True)
    
    x = np.sqrt( binedges[:-1] * binedges[1:] )
    dx = np.vstack((x-binedges[:-1],  binedges[1:]-x))

    ind = (p>0)
    x = x[ind]
    p = p[ind]
    
    return (p, x, dx)


#=============================================

