#import importlib ; import gfHistCompare
#importlib.reload(gfHistCompare) ; gfHistCompare.Test()

def Test(d=3.5, fnameout='test.pdf', debug=True):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    import pandas as pd
          
    prefix = "thick/data_collection"

    nseeds = {
        3.50:  8000,
        7.00:  8000,
        14.00:  300,
        21.00:  300,
        28.00:  300, #120
        56.00:  120,
        84.00:  120,
        112.00: 120
    }

    ngrid = {
        3.50:   2048,
        7.00:   2048,
        14.00:  8192,
        21.00:  8192,
        28.00:  8192, #36864
        56.00:  36864,
        84.00:  36864,
        112.00: 36864
    }

    
    #-- Load data --

    fname = prefix + "_" + str(nseeds[d]) + "_seeds_" + str(ngrid[d]) + ".pkl"   
    label = str(nseeds) + "@" + str(ngrid)
    
    with open(fname, "rb") as fp:
        data_curr = pickle.load(fp)        
    
    xd     = data_curr["x"][d][0]
    yd     = data_curr["y"][d][0]
    numAvg = data_curr["numer_mean"][d][0]  
    thrAvg = data_curr["theor_mean"][d][0]
    totCnt = data_curr["total_count"][d][0]

         
    #-- normalize and analyze and histogram --

    binsizes = xd[1:] - xd[:-1]
    bincenters = np.sqrt( xd[:-1] * xd[1:] )
       
    #pI = yd/np.sum(yd)/binsizes

    #-- plot probability --

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8,4))
    
    iax = ax
 
    iax.plot(bincenters, yd, ':or')
    
    iax.set_yscale('log')
    iax.set_xlim(0,0.5)
    iax.grid()
    
    if debug:

        A = 1
        iax.plot(bincenters/A, bincenters*0+1, 'xg', mfc='none',  markersize=4)
        iax.plot(xd/A, np.ones(xd.shape), '|g',  markersize=8)
       
        sb_edges, sb_cntrs, sb_sizes = binStretch(thrAvg, coeflast=1e3)
        iax.plot(sb_edges/A, np.ones(sb_edges.shape), '|b',  markersize=4)
        print(xd.size, xd[0:1], xd[-1])
        print(sb_edges.size, sb_edges[0:1], sb_edges[-1])

   
    #-----------
    
    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close() 

   
    
#=================================================================

def FitGuess(d=112, fnameout='test.pdf', debug=False):
    import numpy as np
    import matplotlib.pyplot as plt

    
    prefix = "thick/data_collection"

    nseeds = {
        3.50:  8000,
        7.00:  8000,
        14.00:  300,
        21.00:  300,
        28.00:  300, #120
        56.00:  120,
        84.00:  120,
        112.00: 120
    }

    ngrids = {
        3.50:   2048,
        7.00:   2048,
        14.00:  8192,
        21.00:  8192,
        28.00:  8192, #36864
        56.00:  36864,
        84.00:  36864,
        112.00: 36864
    }

    coef = {
        3.50:   (1, 1),
        7.00:   (1, 1),
        14.00:  (1, 1),
        21.00:  (1, 1),
        28.00:  (1, 1),
        56.00:  (1, 1),
        84.00:  (1.154, 1.36, 0.048, 1.45, 1.6),
        112.00: (1.065, 1.15, 0.022, 1.80, 2.0)
    }


    #-- read in data --
        
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,6.5))
    
    marker_style = dict(linestyle=':', color='blue', markersize=6, mfc="none", mec="blue")
    x,y,dx = oneplot(ax[0,0], prefix, ngrids[d], nseeds[d], d, "o", marker_style, debug=debug)
    ax[0,0].cla()
   
    x = x[y>0]
    y = y[y>0]

    y1 =  coef[d][0]*np.exp(-x*coef[d][1])
    y2 =  coef[d][2] * np.exp(-coef[d][3] * (- np.log(x/coef[d][4]))**2 )

    lbl = str(d) + " km, " + str(nseeds[d]) + " seeds @ " + str(ngrids[d]) 
    lbl1 = '{} exp(-{}x)'.format( coef[d][0], coef[d][1] ) 
    lbl2 = '{} exp(-{} ln^2(x/{}))'.format( coef[d][2], coef[d][3], coef[d][4]) 
 
 #-- plot fit in different corrdinates --

    iax=ax[0,0]
    iax.set_xscale('log')
    iax.plot(x, y, ':ob', markersize=4, linewidth=1,  mfc='none', label = lbl)
    iax.plot(x, y1, '-k', markersize=4, linewidth=1,  mfc='none', label = lbl1)
    iax.plot(x, y1+y2, 'r', markersize=4, linewidth=1,  mfc='none')
    iax.set_xlabel('intensity/average')
    iax.set_ylabel('probability and fit')
    iax.grid()
    #iax.legend()


    iax=ax[0,1]
    iax.set_yscale('log')
    iax.plot(x, y, ':ob', markersize=4, linewidth=1,  mfc='none', label = lbl)
    iax.plot(x, y1, '-k', markersize=4, linewidth=1,  mfc='none', label = lbl1)
    iax.plot(x, y1+y2, 'r', markersize=4, linewidth=1,  mfc='none')
    iax.set_xlabel('intensity/average')
    iax.set_ylabel('probability and fit')
    iax.set_ylim(1e-10,10)
    iax.grid()
    iax.legend()


    #-----------------
    
    iax=ax[1,0]
    iax.set_xscale('log')
    iax.plot(x, y-y1, ':ob', markersize=4, linewidth=1,  mfc='none')
    iax.plot(x, y2, '-r', markersize=4, linewidth=1,  mfc='none', label=lbl2)
    iax.set_xlabel('intensity/average')
    iax.set_ylabel('probability - fit')
    iax.grid()
    #iax.legend()
  
    iax=ax[1,1]
    iax.set_xscale('log')
    iax.set_yscale('log')
    iax.plot(x, y-y1, ':ob', markersize=4, linewidth=1,  mfc='none')
    iax.plot(x, y2, '-r', markersize=4, linewidth=1,  mfc='none')
    iax.set_xlabel('intensity/average')
    iax.set_ylabel('probability - fit')
    iax.set_ylim(1e-10,0.1)
    iax.grid()



    #-----------
    
    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close() 


    
#=================================================================


def FitPowerExp(d, fnameout='pdfIdiskFit.pdf', debug=False):
    import numpy as np
    import matplotlib.pyplot as plt

    prefix = "thick/data_collection"

    nseeds = {
        3.50:  8000,
        7.00:  8000,
        14.00:  300,
        21.00:  300,
        28.00:  300, #120
        56.00:  120,
        84.00:  120,
        112.00: 120
    }

    ngrids = {
        3.50:   2048,
        7.00:   2048,
        14.00:  8192,
        21.00:  8192,
        28.00:  8192, #36864
        56.00:  36864,
        84.00:  36864,
        112.00: 36864
    }

    
    xbounds = {
        3.50:   ((0.002,  16), (0.002,  16), (0.6, 16)),
        7.00:   ((0.002,  80), (0.002,  80), (0.2, 80)),
        14.00:  ((0.002, 100), (0.002, 100), (1, 100)),
        21.00:  ((0.002,  80), (0.002,  80), (1, 80)),
        28.00:  ((0.002,  60), (0.002,  60), (1, 60)),
        56.00:  ((0.002,  40), (0.002,  40), (1, 40)),
        84.00:  ((0.002,  40), (0.002,  40), (5, 40)),
        112.00: ((0.002,  40), (0.002,  40), (5, 40))
    }

    #-- read in data --
        
    fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(8,9.5))
    
    marker_style = dict(linestyle=':', color='blue', markersize=6, mfc="none", mec="blue")
    x,y,dx = oneplot(ax[0,0], prefix, ngrids[d], nseeds[d], d, "o", marker_style, debug=debug)
    ax[0,0].cla()
   
    x = x[y>0]
    y = y[y>0]

    lbl = str(d) + " km, " + str(nseeds[d]) + " seeds @ " + str(ngrids[d]) 
    
    #-- plot error as function of left edge of fitting interval --

    x1, e1, q1 = FitError(x,y, xbounds[d][0])

    iax = ax[0,0]
    iax.plot(x1, e1, ':o',  markersize=3, linewidth=1)
    iax.set_xlabel('left edge of fit interval (intensity/average)')
    iax.set_ylabel('fit error')
    iax.set_xlim(0,20)
    iax.grid()

    iax = ax[0,1]
    iax.plot(x1, q1, ':o',  markersize=3, linewidth=1, label = lbl)
    iax.set_xlabel('left edge of fit interval (intensity/average)')
    iax.set_ylabel('q')
    iax.set_xlim(0,20)
    iax.set_ylim(0,1)
    iax.legend()
    iax.grid()

        
    #-- plot fit in different corrdinates --
    
    x1, y1, q1 = FitData(x,y, xbounds[d][1])
    x2, y2, q2 = FitData(x,y, xbounds[d][2])
    lbl1 = 'q = {:4.2f} @'.format(q1) + str(xbounds[d][1])
    lbl2 = 'q = {:4.2f} @'.format(q2) + str(xbounds[d][2])
                 
    iax=ax[1,0]
    iax.plot(x, y, 'ob', markersize=4, linewidth=1,  mfc='none')
    iax.plot(x1, y1, '-y', label=lbl1, linewidth=1)
    iax.plot(x2, y2, '-r', label=lbl2, linewidth=2)
    iax.plot(x1, np.exp(-x1), '--k', linewidth=1)
    iax.set_xlabel('intensity/average)')
    iax.set_ylabel('probability')
    iax.set_xlim(0,2.5)
    iax.legend()
    iax.grid()

        
    iax=ax[1,1]
    iax.set_xscale('log')
    iax.plot(x, y, 'ob', markersize=4, linewidth=1,  mfc='none')
    iax.plot(x1, y1, '-y', label=lbl1, linewidth=1)
    iax.plot(x2, y2, '-r', label=lbl2, linewidth=2)
    iax.plot(x1, np.exp(-x1), '--k', linewidth=1)
    iax.set_xlabel('intensity/average)')
    iax.set_ylabel('probability')
    iax.legend()
    iax.grid()

     
    iax=ax[2,0]
    iax.set_yscale('log')
    iax.plot(x, y, 'ob', markersize=4, linewidth=1,  mfc='none')
    iax.plot(x1, y1, '-y', label=lbl1, linewidth=1)
    iax.plot(x2, y2, '-r', label=lbl2, linewidth=1)
    iax.plot(x1, np.exp(-x1), '--k', linewidth=1)
    iax.set_xlabel('intensity/average)')
    iax.set_ylabel('probability')
    iax.set_ylim(1e-8,10)
    iax.legend()
    iax.grid()

    iax=ax[2,1]
    iax.set_xscale('log')
    iax.set_yscale('log')
    iax.plot(x, y, 'ob', markersize=4, linewidth=1,  mfc='none')
    iax.plot(x1, y1, '-y', label=lbl1, linewidth=1)
    iax.plot(x2, y2, '-r', label=lbl2, linewidth=1)
    iax.plot(x1, np.exp(-x1), '--k', linewidth=1)
    iax.set_xlabel('intensity/average)')
    iax.set_ylabel('probability')
    iax.set_ylim(1e-8,10)
    iax.legend()
    iax.grid()

    #-----------
    
    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close() 

#=================================================================

def FitError (x,y, xbounds):
    from scipy.optimize import least_squares
    import numpy as np
    
    u  = np.log(x)
    v  = np.log(y)
    
    u1 = np.log(xbounds[0])
    u2 = np.log(xbounds[1])
        
    u1range = np.arange(u1,u2,0.1)
    x1range = np.exp(u1range)
    e1range = u1range*0
    q1range = u1range*0
    
    for i in np.arange(0,len(u1range)):
        
        u1 = u1range[i]
        ind = ((u>u1) & (u<u2))
        vb = v[ind]
        ub = u[ind]
    
        c0   = np.array([1, -5, 0.5])
        
        res = least_squares(fun, c0, jac=jac,  args=(ub, vb), verbose=0)

        dv =  model(res.x, ub) - vb
        err = np.sqrt(np.average(dv*dv))
        q1range[i] = res.x[2]
        e1range[i] = err


    print(x1range[-4:], e1range[-4:], e1range.shape[-4:] ) 

    return(x1range, e1range, q1range)

#=================================================================

def FitData (x,y, xbounds):
    from scipy.optimize import least_squares
    import numpy as np
    
    u  = np.log(x)
    v  = np.log(y)
    
    u1 = np.log(xbounds[0])
    u2 = np.log(xbounds[1])
        
    ind = ((u>u1) & (u<u2))
    vb = v[ind]
    ub = u[ind]
        
    c0   = np.array([1, -5, 0.5])
        
    res = least_squares(fun, c0, jac=jac,  args=(ub, vb), verbose=0)

    dv =  model(res.x, ub) - vb
    err = np.sqrt(np.average(dv*dv))

    u_fit = np.arange(u1,u2,0.1)
    v_fit = model(res.x, u_fit)
 
    x_fit = np.exp(u_fit)
    y_fit = np.exp(v_fit)

    q = res.x[2]
    
    return(x_fit, y_fit, q)

                  
#=================================================================

def All(fnameout='pdfIdisk.pdf', iset=0, debug=False, scaled=True, xlog=False):
    import matplotlib.pyplot as plt

    dist = ((3.5, 7), (14, 21), (28, 56), (84, 112))
   
    fig, ax = plt.subplots(ncols=2, nrows=4, figsize=(8,12))

    setcanvas(ax, dist, scaled, xlog)

    #fname = "thick/data_collection_100_seeds_4096.pkl";   nseeds = 100
    #fname = "thick/data_collection_120_seeds_36864.pkl";  nseeds = 120  #OK
    #fname = "thick/data_collection_300_seeds_8192.pkl";   nseeds = 300  #OK
    #fname = "thick/data_collection_8000_seeds_2048.pkl";  nseeds = 8000  #OK

    #fname = "thin/data_collection_300_seeds_8192.pkl";  nseeds = 300  # OK
    #fname = "thin/data_collection_8000_seeds_2048.pkl"; nseeds = 8000 # OK

    if iset == 0:
        showbox = False
    else:
        showbox = True

    if iset == 1 or iset == 0:       
        marker_style = dict(linestyle=':', color='red', markersize=6, mfc="none", mec="red")
        plot_data_set(ax, "thick/data_collection", 8192, 300, dist, "o", marker_style,
                      showbox=showbox, debug=debug, scaled=scaled)
        
    if iset == 2 or iset == 0:
        marker_style = dict(linestyle=':', color='green', markersize=6, mfc="none", mec="green")
        plot_data_set(ax, "thick/data_collection", 2048, 8000, dist, "s", marker_style,
                      showbox=showbox, debug=debug, scaled=scaled)

    if iset == 3 or iset == 0:
        marker_style = dict(linestyle=':', color='blue', markersize=6, mfc="none", mec="blue")
        plot_data_set(ax, "thick/data_collection", 36864, 120, dist, "^", marker_style,
                      showbox=showbox, debug=debug, scaled=scaled)

    if iset == 4 or iset == 0:
        marker_style = dict(linestyle=':', color='y', markersize=6, mfc="none", mec="y")
        plot_data_set(ax, "thick/data_collection", 4096, 100, dist, "v", marker_style,
                      showbox=showbox, debug=debug, scaled=scaled)
        
    if iset == 5 or iset == 0:
        marker_style = dict(linestyle='', color='red', markersize=8, mfc="none", mec="red")
        plot_data_set(ax, "thin/data_collection", 8192, 300, dist, "+", marker_style,
                      showbox=showbox, debug=debug, scaled=scaled)

    if iset == 6 or iset == 0:
        marker_style = dict(linestyle='', color='green', markersize=8, mfc="none", mec="green")
        plot_data_set(ax, "thin/data_collection", 2048, 8000, dist, "+", marker_style,
                      showbox=showbox, debug=debug, scaled=scaled)
    
    
    add_Io500_data(ax[0,0], 'nx1024nz40_z3p50km_probIo.dat', 3.50, scaled)
    add_Io500_data(ax[0,1], 'nx1024nz40_z7p00km_probIo.dat', 7.00, scaled)

    add_Io500_fit(ax[0,0], 3.5, scaled)
    add_Io500_fit(ax[0,1], 7.0, scaled)

    if showbox or xlog:
        loc="lower left"
    else:
        loc="upper right"
        
    add_legends(ax, loc)

    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()

#=================================================================

def plot_data_set(ax, fprefix, ngrid, nseeds, dist, marker, marker_style,
                  showbox=False, debug=False, scaled=True):
    
    import numpy as np

    
    Io_theory = TheoryData("Io")


    for i in (0,1,2,3):
        for j in (0,1):
            
            d = dist[i][j]
    
            oneplot(ax[i,j], fprefix, ngrid, nseeds, d,
                    marker=marker, marker_style=marker_style,
                    showbox=showbox, debug=debug, scaled=scaled)
 
            
#=================================================================

def add_Io500_data(iax, fname, d, scaled):
    import numpy as np


    dat = np.loadtxt(fname)

    x = dat[:,0]
    y = dat[:,1]

    Io_theory = TheoryData("Io")
    
    if scaled:
        A = Io_theory[d]
    else:
        A = 1
    
    iax.plot(x/A, y*A, 'ok', mfc='none',  markersize=2)


#=================================================================

def add_Io500_fit(iax, d, scaled):
    import numpy as np

    Io_theory = TheoryData("Io")
    
    if scaled:
        A = Io_theory[d]
    else:
        A = 1

    ymin, ymax = iax.get_ylim()
           
    P_model, I_model = fit500(d, ymin)
        
    iax.plot(I_model/A, P_model*A, '-k', linewidth=0.5)


#=================================================================

def setcanvas(ax, dist, scaled, xlog):
    import matplotlib.pyplot as plt
    
    Imax = {
        3.50:   30,
        7.00:  100,
        14.00: 120,
        21.00: 120,
        28.00:  80,
        56.00:  50,
        84.00:  50,
        112.00: 50
    }

    Pmin = {
        3.50:   1e-8,
        7.00:   1e-8,
        14.00:  1e-8,
        21.00:  1e-8,
        28.00:  1e-8,
        56.00:  1e-10,
        84.00:  1e-10,
        112.00: 1e-10
    }

    Io_theory = TheoryData("Io")
    
    for i in (0,1,2,3):
        for j in (0,1):
            
            d = dist[i][j]
            A = Io_theory[d]
            
            ax[i,j].set_title("z = {:5.2f} km".format(d))
            ax[i,j].set_yscale('log')
            ax[i,j].set_ylabel('probability')
            #ax[i,j].grid()

            if xlog:
                Imin=1e-2
                ax[i,j].set_xscale('log')
            else:
                Imin=0
            
            if scaled:
                ax[i,j].set_xlabel('intensity / average')
                ax[i,j].set_xlim(Imin,Imax[d])
                ax[i,j].set_ylim(Pmin[d],10)
               
            else:
                ax[i,j].set_xlabel('intensity')
                ax[i,j].set_xlim(0,Imax[d]*A)
                ax[i,j].set_ylim(Pmin[d]/A,10/A)
                ax[i,j].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 

#=================================================================

def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

   
#=================================================================


def oneplot(iax, fprefix, ngrid, nseeds, d, marker, marker_style,
            showbox=False, debug=False, scaled=True):

# fit data only are available only for 3.5 km and 7.0 km
    
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    import pandas as pd
      
    #-- Load data --

    fname = fprefix + "_" + str(nseeds) + "_seeds_" + str(ngrid) + ".pkl"   
    label = str(nseeds) + "@" + str(ngrid)
    
    with open(fname, "rb") as fp:
        data_curr = pickle.load(fp)        
        
    if debug:
        for k in list(data_curr.keys()):
            print(k, data_curr[k].keys())

    if not d in data_curr["x"]:
        return(0)

    xd     = data_curr["x"][d][0]
    yd     = data_curr["y"][d][0]
    numAvg = data_curr["numer_mean"][d][0]  
    thrAvg = data_curr["theor_mean"][d][0]
    totCnt = data_curr["total_count"][d][0]

         
    #-- normalize and analyze and histogram --

    binsizes = xd[1:] - xd[:-1]
    bincenters = np.sqrt( xd[:-1] * xd[1:] )
       
    pI = yd/np.sum(yd)/binsizes

    #-- plot probability --
    
    Io_theory = TheoryData("Io")
    
    if scaled:
        A = Io_theory[d]
    else:
        A = 1

    x = bincenters/A
    y = pI*A
    dx = binsizes/A

    
    iax.plot(x, y, marker, **marker_style, label=label)

    if showbox:
        
        textbox = '\n'.join((
            'Integral:         {:10.1f}'.format(np.sum(pI*binsizes)),
            'Count (data):     {:10.0f}'.format(totCnt),
            'Count (pdf):      {:10.0f}'.format(int(np.sum(yd))),
            'Average (theory): {:10.2e}'.format(thrAvg),
            'Average (direct): {:10.2e}'.format(numAvg),
            'Average (PDF):    {:10.2e}'.format(np.sum(pI*bincenters*binsizes)),
        ))
 
        xmin, xmax = iax.get_xlim()
        ymin, ymax = iax.get_ylim()
        xtext = 0.95*xmax 
        ytext = 0.50*ymax 
    
        props = dict(boxstyle='round', facecolor='white', alpha=1.0)    
        iax.text(xtext, ytext, textbox, fontsize=8, fontfamily='monospace',
                 verticalalignment='top', horizontalalignment='right', bbox=props)
   
    
    if debug:

        iax.plot(bincenters/A, pI*A*0+1, 'xg', mfc='none',  markersize=4)
        iax.plot(xd/A, np.ones(xd.shape), '|g',  markersize=8)
       
        sb_edges, sb_cntrs, sb_sizes = binStretch(thrAvg, coeflast=1e3)
        iax.plot(sb_edges/A, np.ones(sb_edges.shape), '|b',  markersize=4)
        print(xd.size, xd[0:1], xd[-1])
        print(sb_edges.size, sb_edges[0:1], sb_edges[-1])

    return(x, y, dx)

#=================================================================

def model(c, u):
    
    import numpy as np
    return c[0]  +  c[1] * np.exp(c[2]*u)

def fun(c, u, v):
    return model(c, u) - v

def jac(c, u, v):
    import numpy as np

    J = np.empty((u.size, c.size))
    J[:, 0] = 1
    J[:, 1] = np.exp(c[2]*u)
    J[:, 2] = c[1] * u *  np.exp(c[2]*u)
    return J




#=================================================================

def fit500(d, Pmin):
    import numpy as np

   #-- fit coeffeicients for distribution of Io in old 500k set --
    
    coef = {
        2.80: [ 4.68485559,  -30.91654415,   0.8073269],
        3.50: [ 6.125217,    -27.2301761,    0.55827967],
        4.20: [ 5.19359986,  -29.94183211,   0.60971486],
        4.90: [ 5.84976489,  -29.86700398,   0.52638347],
        5.60: [ 6.46256312,  -30.36966598,   0.46885353],
        6.30: [ 7.12712106,  -30.45285702,   0.41602589],
        7.00: [ 7.58833131,  -31.49145846,   0.38884645]
    }


    c = coef[d]
    
    Imax = ((np.log(Pmin) - c[0])/c[1])**(1/c[2])
    
    dI = Imax/1000;

    I = np.arange(dI, Imax, dI)
    
    p = np.exp(  model(c, np.log(I)) )

    return(p, I)

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

    binedges[-1] = 1e8
    return (binedges, bincntrs, binsizes)

#=================================================================

def TheoryData(name):
    
    #-- data for Io from analytical curves --
    
    Io_theory = {
        3.50:  0.016476207362774452,
        7.00:  0.002582854410799337,
        14.00: 0.0003473574209617051,
        21.00: 0.00010215215320075901,
        28.00: 0.00004219865278260401,
        56.00: 4.840305005520772e-6,
        84.00: 1.3441721107398225e-6,
        112.0: 5.395028392456496e-7
    }
   

    half_width_theory = {
        3.50:    6.0,
        7.00:   15.0,
        14.00:  40.0,
        21.00:  75.0,
        28.00: 120.0
    }

    if name == "Io":
        return(Io_theory)
    else:
        return(0)


#=================================================================

