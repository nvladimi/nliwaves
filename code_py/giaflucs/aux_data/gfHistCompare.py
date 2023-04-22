#import importlib ; import gfHistCompare
#importlib.reload(gfHistCompare) ; gfHistCompare.Test(35, 7.0)

def Test(nseeds=35, d=28.0):
# fit data only are available only for 3.5 km and 7.0 km
    
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    import pandas as pd

    print("\n")
    
    dx = 0.1     # grid resolution, cm
    Pmin = 1e-5  # min probability, for plotting

    
    #-- data for Io from analytical curves --
    
    Io_theory = {
        3.50:  0.016476207362774452,
        7.00:  0.002582854410799337,
        14.00: 0.0003473574209617051,
        21.00: 0.00010215215320075901,
        28.00: 0.00004219865278260401
    }


    half_width_theory = {
        3.50:    6.0,
        7.00:   15.0,
        14.00:  40.0,
        21.00:  75.0,
        28.00: 120.0
    }

  
    #-- Load data --

    fname = "data_" + str(nseeds) + "_seeds.pkl"
    
    with open(fname, "rb") as fp:
        data_curr = pickle.load(fp)

    #for k in list(data_curr.keys()):
    #    print(k, data_curr[k].keys())       
    
    db_x = data_curr["x"]
    db = data_curr["intensity"]

    xd = db_x[d][0]
    yd = db[d][0]

 
    #-- borders of circles --
    
    rcirc = 0.2*half_width_theory[d]
    ncirc = int(np.floor(np.pi * (rcirc/dx)**2))
    ntot = ncirc*nseeds
    print("Center circle:", str(rcirc), "cm,",  str(ntot), "points for ",  str(nseeds), ' seeds' )


    #-- analyze and histogram --
    
    binsize = xd[1]-xd[0]
    bincenters = (xd[1:]+xd[:-1])/2
    ncount = sum(yd)
    Iavg = sum(yd*bincenters) / ncount
    
    #print(xd.shape, "xd size")
    #print(yd.shape, "yd size")
    #print(binsize, "bin size")
    #print(sum(yd), "total samples")
    #print(Iavg, "I_average")

    print("Samples mismatch:  ", (sum(yd) - ntot)/ntot)
    print("Intensity mismatch:   ", (Iavg - Io_theory[d])/Io_theory[d])


    #-- plot probability --

    plt.figure(figsize=(6, 6))
    plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel('intensity')
    plt.ylabel('probability')
        
    #plt.bar(bincenters, yd, width=xd[1]-xd[0])
    #plt.plot(bincenters, yd, 'or', markersize=1)
    plt.plot(bincenters, yd/ncount/binsize, 'or', markersize=1)


    sb_edges, sb_cntrs, sb_sizes = binStretch(Io_theory[d])

    plt.plot(sb_edges, np.ones(sb_edges.shape), '|b')    
    
    plt.title("PDF of intensity at z = "+str(d)+" km for circle r = " + str(rcirc) + " cm")

    plt.grid()

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

