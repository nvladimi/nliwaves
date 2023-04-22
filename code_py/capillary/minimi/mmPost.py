#import importlib ; import mmPost
#importlib.reload(mmPost) ; mmPost.Test()


#=================================================================

def avg_klo_all():
    import numpy as np
    import nlwTools

    prefix = "DAT/klo_"
    
    runs= ("m05a1e3", "m07a1e3", "m09a1e3", "m13a1e3", "m19a1e3", "m21a1e3",
           "m05a2e3", "m07a2e3", "m09a2e3", "m13a2e3")

    seeds   = np.arange(0,10)
    seeds19 = (0, 1, 2, 3,    5, 6,    8, 9)
    seeds21 = (0,    2, 3, 4, 5, 6,    8, 9)

    No = 32
    n  = 10000
    
    for r in runs:
        
        fbase = r + '/' + r
        fout = prefix + r + ".npy"
        
        if r == "m19a1e3":
            seeds = seeds19
        elif r == "m21a1e3":
            seeds = seeds21

        spc = nlwTools.AvgSpectrumKlo(fbase, n, No, seeds, istart = 1,  iend = -1)
        np.save(fout, spc)
        
#=================================================================

def avg_seeds_all():
    import numpy as np

    prefix = "dat/sub_"
    
    runs= ("m05a1e3", "m07a1e3", "m09a1e3", "m13a1e3", "m19a1e3", "m21a1e3",
           "m05a2e3", "m07a2e3", "m09a2e3", "m13a2e3")

    seeds   = np.arange(0,10)
    seeds19 = (0, 1, 2, 3,    5, 6,    8, 9)
    seeds21 = (0,    2, 3, 4, 5, 6,    8, 9)
    
    every = 1

    for r in runs:
        
        runname = prefix + r
        outname = prefix + r + ".dat"
        
        if r == "m19a1e3":
            seeds = seeds19
        elif r == "m21a1e3":
            seeds = seeds21
            
        avg_seeds(runname, outname, seeds,   every=every)


#=================================================================

def subset_all():
    import numpy as np

    prefix = "DAT/sub_"
    
    runs= ("m05a1e3", "m07a1e3", "m09a1e3", "m13a1e3", "m19a1e3", "m21a1e3",
           "m05a2e3", "m07a2e3", "m09a2e3", "m13a2e3")

    seeds   = np.arange(0,10)
    seeds19 = (0, 1, 2, 3,    5, 6,    8, 9)
    seeds21 = (0,    2, 3, 4, 5, 6,    8, 9)
    
    every = 100
    
    for r in runs:
        
        runname = r + '/' + r
        outname = prefix + r
        
        if r == "m19a1e3":
            seeds = seeds19
        elif r == "m21a1e3":
            seeds = seeds21
            
        subset_time(runname, outname, seeds,   every=every)

        
#=================================================================

def subset_time(runname, outname, seeds, tmax=1e99, every=100):
    import numpy as np
    import nlwTools
  
    v = nlwTools.data_index()
    
    h=""
    for k in list(v.keys()):
        h = h + "  " + str(v[k]) + "." + k
  
    for seed in seeds:
        
        fname =     runname + "_s{:02d}.dat".format(seed)
        fnameout =  outname + "_s{:02d}.dat".format(seed)

        #print(fname)
        
        dat = np.loadtxt(fname, comments="%")

        tt = dat[:,v["time"]]      
        dat = dat[tt<tmax, :]
        tt = tt[tt<tmax]
 
        ind = nlwTools.CleanTime(tt)       
        dat = dat[ind,:]
        
        ind = np.arange(0,len(dat), every)
        dat = dat[ind,:]

        np.savetxt(fnameout, dat, header = h)

#=================================================================

def avg_seeds(runname, fnameout, seeds, tmax=1e99, every=1):
    import numpy as np
    import nlwTools
  
    v =  nlwTools.data_index()
   
    for seed in seeds:
        
        fname = runname + "_s{:02d}.dat".format(seed)

        print(fname)
        
        dat = np.loadtxt(fname)

        tt = dat[:,v["time"]]      
        dat = dat[tt<tmax, :]
        tt = tt[tt<tmax]
         
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
        
    np.savetxt(fnameout, dat, header = h)

#=================================================================


#=================================================================
