#import importlib ; import nlwTools; import cap512wexac
#importlib.reload(nlwTools) ;  importlib.reload(cap512wexac) ; nlwTools.Test()

#==============================================================================

def Test():

    #test_avg_klo("test/a1", 80, 256,  fnum=1)
    #test_avg_klo("test/a1s55", 160, 512,  fnum=124)
    #test_avg_klo("test/a1s55", 160, 512,  fnum=125)    
    
    #test_klo_vs_psi("test/a1",  80,  fnum=1)
    #test_klo_vs_psi("test/tmp_m05", 6, fnum=2)

    #test_interp_klo("test/a1s55", 160, 512,  fnum=124)

    pass


#==============================================================================



#==============================================================================

def avgpsi_wexac(fbaseout, istart = 1000, iend = 3000):
#
# fbaseout  file name base for averaged interpolated spectrum 
# istart    first time to process
# iend      last file reading to attempt (even non-existing)
    
    import os
    import numpy as np
    import nlwTools as nlw
    
    seeds = ('s50', 's51', 's52', 's53')
 
    No     = 512         # size of original simulation
     
    nfiles = 0

    nk = np.zeros((No,No))
    
    for s in seeds:

        fbase_seed =  s + "/" + "a1" + s
        
        for ifile in np.arange(istart, iend):

            fname = fbase_seed + '.psi.' + str(ifile).zfill(4)
            #print(fname)

            if os.path.isfile(fname):

                nk += nlw.ReadPsi(fname, spectrum = True)
                
                nfiles += 1             
                msg = "For seed {} last file read: {}".format(s, fname)
            
        print("{}. Total {} file read".format(msg, nfiles))

    nk = nk/nfiles
    
    np.save(fbaseout+"_nk.npy", nk)

    print("Number of files read: {}.".format(nfiles) )


#==============================================================================

def average_wexac(fbaseout, istart = 1, iend = 10, dNk=1, dNphi=120):
#
# fbaseout  file name base for averaged interpolated spectrum 
# istart    first time to process
# iend      last file reading to attempt (even non-existing)
# dNk       spacing in |k| for polar grid
# dNphi     number of angles in polar grid
    
    import os
    import numpy as np
    import nlwTools as nlw
    
    seeds = ('s50', 's51', 's52', 's53')
 
    No     = 512         # size of original simulation
    Nklo   = 160         # size of reduced data
    
    nfiles = 0
    nslices  = 0

    nlw.KloAverageClear()

    for s in seeds:

        fbase_seed =  s + "/" + "a1" + s
        
        for ifile in np.arange(istart, iend):

            fname = fbase_seed + '.klo.' + str(ifile).zfill(4)
            #print(fname)

            if os.path.isfile(fname):

                klo = nlw.ReadKlo(fname, Nklo, No, spectrum = False)

                nlw.KloAverageAdd(klo, dNk=dNk, Nphi=32)
                
                nslices += klo.shape[0]
                nfiles += 1             
                msg = "For seed {} last file read: {}".format(s, fname)
            
        print("{}. Total {} slices in {} files".format(msg, nslices, nfiles))
           


    a, b, count = nlw.KloAverageGet()
    
    np.save(fbaseout+"_kxky.npy", a)
    np.save(fbaseout+"_kphi.npy", b)

    print("Number of slices in {} files read: {} ".format(nfiles,  nslices) )

#==============================================================================

def evolution_wexac(foutbase, iend=10, every=1):
#
#  fbaseout  text output
#  iend      last file reading to attempt (even non-existing) 
#  every     process one in every so many slices
    
    import numpy as np
    import os
    import nlwTools as nlw

    seeds = ('s50', 's51', 's52', 's53',  's54',  's55')
 
    No     = 512         # size of original simulation
    Nklo   = 160         # size of reduced data
    dt     = 0.1         # time between snapshots
    
    for s in seeds:

        fbase_seed =  s + "/" + "a1" + s
        fname_out = foutbase + "_" + s + ".txt"
        nslices = 0
        nfiles  = 0
        
        for ifile in np.arange(1, iend):

            fname = fbase_seed + '.klo.' + str(ifile).zfill(4)
            #print(fname)

            if os.path.isfile(fname) :

                klo = nlw.ReadKlo(fname, Nklo, No, spectrum = False)
                t0 = nslices * dt
                
                nlw.KloEvolution(klo, t0, dt, every, fname_out)

                nslices += klo.shape[0]
                nfiles += 1
                msg = "For seed {} last file read: {}".format(s, fname)
                

        print("{}. Number of slices in {} files read: {}".format(msg, nfiles, nslices))


#==============================================================================
