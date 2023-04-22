#import importlib ; import minimi
#importlib.reload(minimi) ; minimi.All()


#=================================================================

def plotAll(ampname="a1e3", fnameout='test.pdf', every = 1):
    import matplotlib.pyplot as plt
    import numpy as np

    prefix = "DAT/"
    modenames = (("m05", "m07"), ("m09", "m13"), ("m19", "m21"))
    seeds = np.arange(0,4)

    
    fig, ax = plt.subplots(ncols=2, nrows=4, figsize=(9,12))

    setcanvas(ax)

    marker_etamin = dict(linestyle='-', linewidth=0.5, color='b', markersize=6, mfc="none", mec="y")
    marker_etamax = dict(linestyle='-', linewidth=0.5, color='r', markersize=6, mfc="none", mec="y")
    marker_etarms = dict(linestyle='-', linewidth=0.5, color='g', markersize=6, mfc="none", mec="y")

    for i in (0,1,2):
        for j in (0,1):
            runname = prefix + modenames[i][j] + ampname + ".dat"
            plot_data(ax[i,j], runname, 'eta_min', '-', marker_etamin, every=every)
            plot_data(ax[i,j], runname, 'eta_max', '-', marker_etamax, every=every)
            plot_data(ax[i,j], runname, 'eta_rms', '-', marker_etarms, every=every)

            plot_data(ax[3,0], runname, 'H',  '-', marker_etamin, every=10)
            plot_data(ax[3,1], runname, 'NL', '-', marker_etamin, every=10)

            
    #add_legends(ax, loc)

    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()


#=================================================================

def setcanvas(ax):
    import matplotlib.pyplot as plt

    modes = (("5x5", "7x7"), ("9x9", "13x13"),  ("19x19", "21x21"))
    
    for i in (0,1,2):
        for j in (0,1):

            ax[i,j].set_title(modes[i][j])
            ax[i,j].set_xlabel('$t$')
            ax[i,j].set_ylabel('$\eta$')
            ax[i,j].grid()
            ax[i,j].set_ylim(-1,1)
            ax[i,j].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 

    ax[3,0].set_title("all runs")
    ax[3,0].set_xlabel('$t$')
    ax[3,0].set_ylabel('$H$')
    ax[3,0].grid()
    ax[3,0].set_ylim(0,1.5)
    ax[3,0].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 
    
    ax[3,1].set_title("all runs")
    ax[3,1].set_xlabel('$t$')
    ax[3,1].set_ylabel('$H_{NL} / H$')
    ax[3,1].grid()
    ax[3,1].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 

            
#=================================================================
   
def plot_data(iax, fname, var, marker, marker_style, every=1):
    import numpy as np
    
    v = data_index()
        
    dat = np.loadtxt(fname, comments="#")   
    ind = np.arange(0,len(dat), every)
    dat = dat[ind,:]
 
    t = dat[:,v["time"]]
        
    if (var == "eta_min" or var == "eta_max" or var == "eta_rms"):
        f = dat[:,v[var]]
    if (var == "H"):
        f = dat[:,v["E_pot"]] + dat[:,v["E_kin"]] + dat[:,v["E_nl"]]
    if (var == "NL"):
        h = dat[:,v["E_pot"]] + dat[:,v["E_kin"]] + dat[:,v["E_nl"]]
        f = dat[:,v["E_nl"]]/h

    #print(dat.shape, t.shape, f.shape)
    iax.plot(t, f, marker, **marker_style)
   
#=================================================================

def avg_seeds_all():
    import numpy as np

    prefix = "DAT/"
    modenames = ("m05", "m07", "m09", "m13", "m19", "m21")
 
    tmax_test = (100, 100,  100, 100,  100, 100)
    tmax1e3   = (0.9e+5, 0.9e+5, 1.3e+5, 0.9e+5, 1.1e+5, 0.8e+5)
    tmax2e3   = (0.9e+5, 1.3e+5, 1.2e+5, 1.1e+5, 0.7e+5, 1.3e+5)

    tmax2e3   = (1.3e+5, 1.2e+5, 1.1e+5, 0.7e+5, 1.3e+5)
    modenames = ("m07", "m09", "m13", "m19", "m21")

    
    #ampname = "a1e3";  tmax = tmax1e3
    ampname = "a2e3";  tmax = tmax2e3

    seeds = np.arange(0,10)
    every = 100
    
    for i in range(len(modenames)):
        runname = prefix + modenames[i] + ampname
        avg_seeds(runname, seeds, tmax[i], every=every)

#=================================================================

    
def avg_seeds(runname, seeds, tmax, every=1):
    import numpy as np
  
    v = data_index()
 
    for seed in seeds:
        
        fname = runname + "_s{:02d}.dat".format(seed)

        print(fname)
        print("tmax = ", tmax)
        
        dat = np.loadtxt(fname, comments="%")

        tt = dat[:,v["time"]]      
        dat = dat[tt<tmax, :]
        tt = tt[tt<tmax]
 
        ind = clean_time(tt)       
        dat = dat[ind,:]
        
        ind = np.arange(0,len(dat), every)
        dat = dat[ind,:]

        if 't' in locals():
            eta_min = np.minimum(eta_min, dat[:,v["eta_min"]])
            eta_max = np.maximum(eta_max, dat[:,v["eta_max"]])
            eta_sq  += (dat[:,v["eta_rms"]])**2
            E_pot   += dat[:,v["E_pot"]]
            E_kin   += dat[:,v["E_kin"]]
            E_nl    += dat[:,v["E_nl"]]
        else:
            t = dat[:,v["time"]]
            eta_min = dat[:,v["eta_min"]]
            eta_max = dat[:,v["eta_max"]]
            eta_sq  = (dat[:,v["eta_rms"]])**2
            E_pot   = dat[:,v["E_pot"]]
            E_kin   = dat[:,v["E_kin"]]
            E_nl    = dat[:,v["E_nl"]]

    nseeds = len(seeds)
    
    eta_rms = np.sqrt(eta_sq/nseeds)
    E_pot = E_pot/nseeds
    E_kin = E_kin/nseeds
    E_nl  = E_nl/nseeds
    
    dat = np.zeros((len(t), len(v)))
    print("data out:", dat.shape, "\n")
    
    
    dat[:,v["time"]] = t
    dat[:,v["eta_min"]] = eta_min
    dat[:,v["eta_max"]] = eta_max
    dat[:,v["eta_rms"]] = eta_rms
    dat[:,v["E_pot"]]   = E_pot
    dat[:,v["E_kin"]]   = E_kin
    dat[:,v["E_nl"]]    = E_nl

    h=""
    for k in list(v.keys()):
        h = h + "  " + str(v[k]) + "." + k
        
    fnameout = runname + ".dat"
    np.savetxt(fnameout, dat, header = h)

#=================================================================

def clean_time(t, debug = False):
# input:   t - 1D array of time from *.dat file
# output:  1D array of indexes for consecutive time

    import numpy as np

    #-- find left edges of intervals  --
    
    i1 = np.asarray(np.nonzero(t[:-1] > t[1:])).squeeze() + 1

    if debug:
        print("i1 = ", i1, ",   shape = ", i1.shape)

    if i1.shape == (0,):
        return(np.array(len(t)))

        
    #-- find rigth edges of intervals --
        
    i2 = []

    if i1.shape == ():
        try:
            q = np.max(np.nonzero(t<t[i1]))
        except:
            q = 0
        i2.append(q)

    else:
        for i in i1:
            try:
                q = np.max(np.nonzero(t<t[i]))
            except:
                q = 0
            i2.append(q)

    if debug:
        print("i2 = ", i2)

    #-- find intervals between jumps --
        
    if i2[0] == 0:
        i2 = i2[1:]
        i2.append(len(t)-1)
    else:
        i2.append(len(t)-1)
        i2 = np.asarray(i2) 
        i1 = np.hstack((np.asarray(0),i1))
    
    ind = np.asarray([],dtype=int)
    for i in range(len(i1)):
        ind = np.hstack((ind, np.arange(i1[i], i2[i]+1)))

    if debug:
        print("i1 = ", i1)
        print("i2 = ", i2)
        print(t[ind])

    return(ind)
    
#=================================================================

def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

   
#=================================================================

def data_index():
    
    v = {
        "time"    :   0,
        "E_pot"   :   1,  
        "E_kin"   :   2,
        "E_nl"    :   3,
        "eta_rms" :   4,
        "eta_min" :   5,
        "eta_max" :   6,
        "sumA_sq" :   7,
        "sumAk_sq":   8,
        "Px"      :   9,
        "Py"      :  10
    }

    return(v)

#=================================================================

#=================================================================

