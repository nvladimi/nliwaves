#import importlib ; import giaflucs
#importlib.reload(giaflucs) ; giaflucs.ringPDF()

def ringImax():

    from copy import copy
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    import os
    #import pickle

    inFileName = "disk.npy"                              # intensity data for given z 
    inRingBorders = [200, 400, 600, 800]                 # input for given z, in pixels
    inRingThresholds = [1.0e-4, 1.5e-4, 2.0e-4, 2.5e-4]  # input from analytical curve for given z
    inMinDist = 10                                       # minimal distance b/w specles, in pixels

    debug = False

    x1=300; x2=500; y1=300; y2=500;                      # zoom area for debugging images
    figDPI=216                                           # resolution for debugging images
    
#--------------------------------------------

#-- predefine some colors for visuaizarion --

    if debug:

        #plt.ion()

        graymasked = copy(plt.cm.gray);
        graymasked.set_over('y', 1.0);
        graymasked.set_under('b', 1.0);
        graymasked.set_bad('r', 1.0)

        jetmasked = copy(plt.cm.jet);
        jetmasked.set_under('k', 1.0);
        jetmasked.set_over('m', 1.0);
        jetmasked.set_bad('m', 1.0)

    
#-- read and show field data --

    #-- array of field intensity, keep

    f = np.load(inFileName)

    if debug:
        plt.imshow(np.log(f), cmap='jet', vmin=-12, vmax=-7,  origin='lower'); # plt.show()
        plt.savefig('disk-jet.png', pad_inches=0, dpi=figDPI)

        #plt.imshow(f, cmap='gray', vmin=0, vmax=4e-4, origin='lower');  #plt.show()
        #plt.savefig('disk-gray.png')

    
#-- find all local maxima in whole array --

    f0 = f[1:-1, 1:-1]
    fE = f[1:-1, 2:]
    fW = f[1:-1, 0:-2]
    fN = f[0:-2, 1:-1]
    fS = f[2:,   1:-1]

    q0 = ( (f0>fN) & (f0>fS) & (f0>fE) & (f0>fW) )
    
    del f0, fE, fW, fN, fS

    #-- array of booleans with all local maxima, keep

    q = np.full(f.shape, False)

    q[1:-1,  1:-1] = q0

    del q0

    if debug:

        fm = np.ma.masked_where( (q == True), f)
                
        plt.imshow(np.log(fm[x1:x2,y1:y2]), extent=(x1,x2,y1,y2),
                   cmap=jetmasked, vmin=-12, vmax=-7, origin='lower');
        plt.show();  #debug
        plt.savefig('disk_zoom_allmax.png', pad_inches=0, dpi=figDPI)
 
        
#-- create array of distances to the center in pixels ---
#-- IMPORTANT: MODIFY CODE TO REFLECT CORRECT LOCATION OF ORIGIN ---

    (ixmax, iymax) = f.shape
    ix = np.arange(-np.floor(ixmax/2), np.ceil(ixmax/2))
    iy = np.arange(-np.floor(iymax/2), np.ceil(iymax/2))

    #-- arrays with coordinates, keep
    (x, y) = np.meshgrid(ix,iy) 

    #-- array of distances to the center, keep
    r = np.sqrt(x*x + y*y)

    # plt.imshow(r);  plt.show()  #debug

    fmax = f[q]
    xmax = x[q]
    ymax = y[q]
 
    print(fmax.shape,  q.shape)  # debug
    return([fmax, xmax, ymax])
   
    
    
#-- this should be a loop over rings, testing for one ring only --

    ring=2

    #-- array of booleans with all local maxima in the ring
    
    qring = ( q & (f>inRingThresholds[ring])  & (r>inRingBorders[ring-1])  & (r<inRingBorders[ring]) )

    if debug:
        rring = ( (r>inRingBorders[ring-1])  & (r<inRingBorders[ring]) )

        fm = np.ma.masked_where( (qring == True), f)
        fr = np.log(fm) - 2*(rring == False)
        
        plt.imshow(fr[x1:x2,y1:y2], extent=(x1,x2,y1,y2),
                   cmap=jetmasked, vmin=-12, vmax=-7, origin='lower');
        #plt.show();  #debug
    
        plt.savefig('disk_zoom_ringmax.png', pad_inches=0, dpi=figDPI)
     
    indring = np.nonzero(qring)

    fmax = f[indring]
    xmax = x[indring]
    ymax = y[indring]
    
    ismax = np.full(fmax.shape, True)  # array indicating which maxima are isolated

    #-- find false maxima that are too close to others, set corrsponding ismax to false

    fmaxsort = np.flip(np.argsort(fmax)) 
     
    for i in fmaxsort:
        if (ismax[i]):
             x0 = xmax[i]
             y0 = ymax[i]
             # might be faster if we compute distance only to maxima in the square, something like this  
             # ind = ( (abs(xmax - x0) < inMinDist ) & (abs(ymax - y0) < inMinDist ) & (ismax > 0) )
             dx = xmax - x0
             dy = ymax - y0
             drr = dx*dx + dy*dy
             ind =  ( drr < inMinDist*inMinDist )
             ind[i] = False
             ismax[ind] = False
     
    #-- return array of selected maxima
    
    fmaxout = fmax[ismax]
    
    #print(fmax.shape,  fmaxout.shape)  # debug
    
    #-- images for verification and debugging

    if debug:

        isnot =  np.logical_not(ismax)

        indringtrue  = (indring[0][ismax], indring[1][ismax]) 
        indringfalse = (indring[0][isnot], indring[1][isnot]) 

        #fm = np.ma.masked_where( (indringtrue == True), f)
        fm = f;

        fm[indringfalse] = 1e-20
        fm[indringtrue] = 1

        fr = np.log(fm) - 2*(rring == False)
        
        plt.imshow(fr[x1:x2,y1:y2], extent=(x1,x2,y1,y2),
                   cmap=jetmasked, vmin=-12, vmax=-7, origin='lower');
        #plt.show();  #debug
    
        plt.savefig('disk_zoom_select.png', pad_inches=0, dpi=figDPI)
   
    
    return(fmaxout)

#=================================================================================================
def ringPDF():

    from copy import copy
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    import os
    #import pickle

    inFileName = "disk.npy"                              # intensity data for given z 
    inRingBorders = [0, 300, 600]                           # input for given z, in pixels

    x1=0; x2=400; y1=0; y2=400;                      # zoom area for debugging images
    figDPI=216                                           # resolution for debugging images

    debug = True

    
#--------------------------------------------
#-- predefine some colors for visuaizarion --

    if debug:

        #plt.ion()

        graymasked = copy(plt.cm.gray);
        graymasked.set_over('y', 1.0);
        graymasked.set_under('b', 1.0);
        graymasked.set_bad('r', 1.0)

        jetmasked = copy(plt.cm.jet);
        jetmasked.set_under('k', 1.0);
        jetmasked.set_over('m', 1.0);
        jetmasked.set_bad('m', 1.0)

  
    
#-- read and show field data --

    #-- array of field intensity, keep

    f = np.load(inFileName)

    if False:
        plt.imshow(np.log(f), cmap='jet', vmin=-12, vmax=-7,  origin='lower');  #plt.show()
        #plt.savefig('disk-jet.png', pad_inches=0, dpi=figDPI)

        #plt.imshow(f, cmap='gray', vmin=0, vmax=4e-4, origin='lower');  #plt.show()
        #plt.savefig('disk-gray.png')




        
#-- create array of distances to the center in pixels ---
#-- IMPORTANT: MODIFY CODE TO REFLECT CORRECT LOCATION OF ORIGIN ---

    (ixmax, iymax) = f.shape
    ix = np.arange(-np.floor(ixmax/2), np.ceil(ixmax/2))
    iy = np.arange(-np.floor(iymax/2), np.ceil(iymax/2))

    #-- arrays with coordinates, keep
    (x, y) = np.meshgrid(ix,iy) 

    #-- array of distances to the center, keep
    r = np.sqrt(x*x + y*y)

    # plt.imshow(r);  plt.show()  #debug
    
    #-- this should be a loop over rings, testing for one ring only --

    ring = 0

    #-- all intensity in the ring --
    
    rring = ( (r>inRingBorders[ring])  & (r<inRingBorders[ring+1]) )
    Iring = f[rring];
    
    Iavg = np.average(Iring)

    #-- histogram for equal bins --
    
    numbins = 2000                                              # input parameter
    binsize = 0.1*Iavg                                          # input parameter
    
    binedges = np.arange(numbins)*binsize
    (P1, x) = np.histogram(Iring, binedges, density=True)
    x1 = 0.5 * (binedges[:-1] + binedges[1:])                   # bin centers for plotting

    #-- histogram for stretched bins --

    numbins  = 50                                               # input parameter
    binsize0 = 1e-2*Iavg                                        # input parameter
    lastedge = 1e+3*Iavg                                        # input parameter
    
    c = np.log(1+ lastedge/binsize0)/numbins
    binedges = (np.exp(c * np.arange(numbins+1)) - 1)*binsize0
    (P2, x) = np.histogram(Iring, binedges, density=True)
    x2 = np.sqrt( binedges[:-1] * binedges[1:] )                # bin centers for plotting

    x2err = np.vstack((x2-binedges[:-1],  binedges[1:]-x2))     # errorbars to bin edges
    
    
    #-- plot both histrograms --
    
    fig, ax = plt.subplots();
    plt.yscale('log')   
    ax.set_xlim(0,0.0005)
    
    ax.plot(x1,P1, 'ro', markersize=2, linewidth=1);
    ax.errorbar(x2, P2, xerr=x2err, markersize=4,  fmt='o', mfc='none', linewidth=1)
    
    ax.set(xlabel='intensity', ylabel='probability')
    plt.show()


    

    if False:

        fr = np.log(f) - 2*(rring == False)

        plt.imshow(fr, cmap='jet', vmin=-12, vmax=-7, origin='lower');
        plt.show();  #debug
    
        #plt.savefig('disk_zoom_ringmax.png', pad_inches=0, dpi=figDPI)
     
