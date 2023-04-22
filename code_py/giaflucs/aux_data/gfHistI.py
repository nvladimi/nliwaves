#import importlib ; import gfHistI
#importlib.reload(gfHistI) ; gfHistI.Test()

def Test(inFileName = "disk.npy",  r1=0, r2=30, dx=0.1):
    
    import numpy as np

    Iring =  SelectIntesityInRing(inFileName, r1/dx, r2/dx, debug=False)

    # average intensity - replace by theoretical
    
    Iavg = np.average(Iring)
    
    ringHist(Iring, Iavg, density=False)

#=================================================================


def ringHist(Iring, avgItheory, density=False):
    
    import numpy as np
    import matplotlib.pyplot as plt

    #-- histogram for equal bins --
    
    numbins = 2000                                              # input parameter
    binsize = 0.1*avgItheory                                    # input parameter
    
    binedges = np.arange(numbins+1)*binsize
    (P1, x) = np.histogram(Iring, binedges, density=density)
    x1 = 0.5 * (binedges[:-1] + binedges[1:])                   # bin centers for plotting
    dx1 = binedges[1:] - binedges[:-1]                          # bin sizes for normalization  

    #-- histogram for stretched bins --

    binedges, bincntrs, binsizes = binStretch(avgItheory)
    (P2, x) = np.histogram(Iring, binedges, density=density)
    x2 = bincntrs                                              # bin centers for plotting
    dx2 = binedges[1:] - binedges[:-1]                         # bin sizes for normalization 
    

    #-- normalize histogram --

    ntot = Iring.size
    
    if density == False:
        n1 = np.sum(P1)
        n2 = np.sum(P2)
        P1 = P1/ntot/dx1
        P2 = P2/ntot/dx2
        txt0 = 'Totals: {}, {}, {}'.format(Iring.size, n1, n2)
        
    txt1 = 'Integrals: {}, {}'.format(np.sum(P1*dx1), np.sum(P2*dx2))
    txt2 = 'Averages:  {:8.2e},  {:8.2e}, {:8.2e}, {:8.2e}'.format(
        avgItheory, np.average(Iring), np.sum(P1*x1*dx1), np.sum(P2*x2*dx2)) 
    
    #-- plot probability --

    plt.figure(figsize=(6, 6))
    plt.yscale('log')
    plt.xlabel('intensity')
    plt.ylabel('probability')
    plt.xlim(0,0.0005)
    plt.ylim(0.5,1e5)
    
    plt.plot(x1, P1, 'or', mfc='none',  markersize=2)
    plt.plot(x2, P2, 'ob', mfc='none',  markersize=6)

    a=36000
    plt.plot(x2, a*np.exp(-x2*a), '-g',  linewidth = 1)
    plt.plot(binedges, np.ones(binedges.shape), '|b')    
    
    plt.grid()

    if density == False:
        plt.text(0.00012, 4e4, txt0)

    plt.text(0.00012, 2e4, txt1)
    plt.text(0.00012, 1e4, txt2)
    
    plt.show()


#=================================================================

def binStretch(favg,  numbins=50,  coef0=1e-2, coeflast=1e2):

    # inputs:
    # favg     - average of the distribution
    # coef0    - size of the first bin, in terms of average 
    # coeflast - location of the last bin, in terms of average
    
    import numpy as np

    binsize0 = coef0*favg
    lastedge = coeflast*favg
    c = np.log(1+ lastedge/binsize0)/numbins
    
    binedges = (np.exp(c * np.arange(numbins+1)) - 1)*binsize0
    bincntrs = np.sqrt( binedges[:-1] * binedges[1:] )
    binsizes = np.vstack((bincntrs-binedges[:-1],  binedges[1:]-bincntrs))
  
    return (binedges, bincntrs, binsizes)

#=================================================================


def SelectIntesityInRing(inFileName, r1, r2, debug=True):

    import numpy as np

    f = np.load(inFileName)

    (ixmax, iymax) = f.shape
    ix = np.arange(-np.floor(ixmax/2), np.ceil(ixmax/2))
    iy = np.arange(-np.floor(iymax/2), np.ceil(iymax/2))

    #-- arrays with coordinates
    (x, y) = np.meshgrid(ix,iy) 

    #-- array of distances to the center
    r = np.sqrt(x*x + y*y)

    #-- all intensity in the ring --
    
    rring = ( (r>r1)  & (r<r2) )
    Iring = f[rring];

    if debug:

        from copy import copy
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        #plt.ion()

        graymasked = copy(plt.cm.gray);
        graymasked.set_over('y', 1.0);
        graymasked.set_under('b', 1.0);
        graymasked.set_bad('r', 1.0)

        jetmasked = copy(plt.cm.jet);
        jetmasked.set_under('k', 1.0);
        jetmasked.set_over('m', 1.0);
        jetmasked.set_bad('m', 1.0)

        fr = np.log(f) - 2*(rring == False)

        plt.imshow(fr, cmap='jet', vmin=-12, vmax=-7, origin='lower');
        plt.show();  #debug
    
        #plt.savefig('disk_zoom_ringmax.png', pad_inches=0, dpi=figDPI)
     
    return(Iring)

#=================================================================
