
#import importlib ; import fiboPX
#importlib.reload(fiboPX) ; fiboPX.Test()

def Test():

    fbase = "q7_v0dt1"
    fbaseout= 'out/q7_v0dt1'    
    #probX(fbase, istart=1, iend=-1, binsize=1, numbins=20, fbaseout=fbaseout)


    fbase = "q7_v0dt1_s1"
    fbaseout= 'out/q7_v0dt1_s1'
    imode = 25
    #corrSigma(fbase, imode, jmax=7, tmaxM=1/16., ntaumax = 1024, istart=1, iend=-1,  fbaseout=fbaseout)

    corrZero(fbase, imode, jmax=11, istart=1, iend=-1,  fbaseout='')

    
#==============================================================================
    
def datasets(fbase):
    
    import numpy as np
    import os
    
    sets = list()
    
    for file in os.listdir("."):
        if file.startswith(fbase):
            sets.append(file)
            #print(os.path.join("/mydir", file))

    [m, nsave] = ReadParam(sets[0])

    return(m, nsave, sets)
    

#==============================================================================


def probX(fbase, istart=1, iend=-1, binsize=0.1, numbins=200, fbaseout=''):
#
#   
    import numpy as np
    import pickle
    import os

    nmodes, nsave, sets = datasets(fbase)
    
    edg = (np.arange(-numbins/2, numbins/2+1))*binsize;
    cnt = (np.arange(-numbins/2, numbins/2)+0.5)*binsize;
    edg[0]  = -1.e16
    edg[-1] =  1.e16
 
    
    P    = np.zeros((numbins, nmodes), 'uint32')
    nk   = np.zeros(nmodes)

    ntot = 0

    ire = np.arange(nmodes)*2 + 1
    iim = ire + 1;
   

    #-- first round: compute average --
        
    for s in sets:

        ifile = istart;
        
        while True:

            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break
  
            aak = dat[:,ire]**2 + dat[:,iim]**2

            nk += np.sum(aak,0)
            ntot += aak.shape[0]
   
    Xk0 = 0.5 * np.log(nk/ntot)

    
   #-- second round: compute probability --
        
    for s in sets:

        ifile = istart;
        
        while True:

            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break
            
            Xk = 0.5 * np.log(dat[:,ire]**2 + dat[:,iim]**2) #  Xk[:, nmodes]
            
            
            for m in range(nmodes):
                h, x = np.histogram(Xk[:,m] - Xk0[m], edg)
                P[:,m] += h.astype('uint32')

    ncount = np.sum(P,0)


    #-- output --

    dictout = {
              "Xhist":     P,
              "centers":   cnt,
              "edges":     edg,
              "ntot":      ntot
              }
    
    with open(fbaseout + "_probX.pkl", 'wb') as handle:
        pickle.dump(dictout, handle, protocol=pickle.HIGHEST_PROTOCOL)

        
    #-- images for gebugging --
  
    if False:
        import matplotlib.pyplot as plt
             
        plt.figure(figsize=(6,6))

        plt.yscale('log')
        
        plt.plot(cnt, P[:,10], 'bo-');
        plt.plot(cnt, P[:,20], 'rx-');
        plt.plot(cnt, P[:,30], 'gs-');
        plt.grid(True)

        plt.show()

    

#==============================================================================



def corrZero(fbase, imode, jmax=11, istart=1, iend=-1,  fbaseout=''):
#
#   
    import numpy as np
    import pickle
    import os

    nmodes, nsave, sets = datasets(fbase)
 
    Savg   = np.zeros(jmax)
    Corr   = np.zeros(jmax)

    ntot = 0

    ire = np.arange(nmodes)*2 + 1
    iim = ire + 1;
   

    
    #-- first round: compute average --
        
    for s in sets:

        ifile = istart;
        
        while True:

            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break


            rho = np.sqrt(dat[:,ire]**2 + dat[:,iim]**2)

            S = np.log(  rho[:, imode:imode+jmax] / rho[:, imode-1:imode+jmax-1]) 
 
            Savg += np.sum(S,0)
            ntot += S.shape[0]
   
    Savg = Savg/ntot

    
   #-- second round: compute correlation --
        
    for s in sets:

        ifile = istart;
        
        while True:

            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break

            rho = np.sqrt(dat[:,ire]**2 + dat[:,iim]**2)

            S = np.log(  rho[:, imode:imode+jmax] / rho[:, imode-1:imode+jmax-1]) 

            for j in np.arange(jmax):
                Corr[j] += sum((S[:,j] - Savg[j])*(S[:,0] - Savg[0]), 0)

   
    Corr = Corr/ntot

  
    #-- output --

    dictout = {
              "imode":     imode,
              "jmax":   jmax,
              "Sigma":     Savg,
              "Corr":      Corr,
              "ntot":      ntot,
              }


    with open(fbaseout + "_sigmacorrzero.pkl", 'wb') as handle:
        pickle.dump(dictout, handle, protocol=pickle.HIGHEST_PROTOCOL)

 



#==============================================================================


def corrSigma(fbase, imode, jmax=7, tmaxM=1, ntaumax = 1024, istart=1, iend=-1,  fbaseout=''):
    
    import numpy as np
    import os
    import pickle
    
    [nmodes, nsave] = ReadParam(fbase)
    
    ire = np.arange(nmodes)*2 + 1
    iim = ire + 1;

    ntmax = round(tmaxM*1024*1024)

    S = np.zeros((ntmax, jmax))
    
    #-- read in data --

    n0 = 0
    
    ifile = istart;

    while n0 < ntmax:

        fname = fbase + '/'+ fbase + '.' + str(ifile).zfill(4) + '.ak'

        if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
            dat = np.fromfile(fname, 'float64')
            dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
            ifile += 1
        else:
            #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
            break

        rho = np.sqrt(dat[:,ire]**2 + dat[:,iim]**2)

        n1 = rho.shape[0]
        if n0+n1 > ntmax:
            n1 = ntmax - n0

        S[n0:n0+n1, :]  = np.log(  rho[:n1, imode:imode+jmax] / rho[:n1, imode-1:imode+jmax-1]) 
        n0 += n1

    #-- subtract average --

    for j in range(jmax):
        S[:,j] += -np.average(S[:,j])
    
            
    #-- correlations --

    CC  = np.zeros((2*ntaumax, jmax))
    ckj = np.zeros(jmax)
    cjj = np.zeros(jmax)

    n1 = round(ntmax/2) - ntaumax
    n2 = round(ntmax/2) + ntaumax
    
    for j in range(jmax):
        C = np.correlate(S[:,0], S[:,j], 'same')
        CC[:,j] = C[n1:n2]
        ckj[j] = np.sum(S[:,0] * S[:,j])
        cjj[j] = np.sum(S[:,j] * S[:,j])


    #-- output --

    dictout = {"imode":    imode,
               "ntmax":    ntmax,
               "corrfun":  CC,
               "ckj":      ckj,
               "cjj":      cjj
    }
    
    with open(fbaseout + "_sigmacorr.pkl", 'wb') as handle:
        pickle.dump(dictout, handle, protocol=pickle.HIGHEST_PROTOCOL)
 

    #-- images for gebugging --
        
    if False:
        import matplotlib.pyplot as plt
             
        plt.figure(figsize=(10,4))
        t = np.arange(0,ntmax,100)
        print (t.shape, S.shape)
        
        plt.plot(t, S[t,0], 'b-', lw=0.5);
        plt.plot(t, S[t,1], 'r-', lw=0.5);
        plt.plot(t, S[t,2], 'g-', lw=0.5);
        plt.grid(True)

        plt.show()

    if False:
        import matplotlib.pyplot as plt
             
        plt.figure(figsize=(10,4))
        t = np.arange(0,2*ntaumax) - ntaumax
        print (t.shape, S.shape)
        
        plt.plot(t, CC[:,0], 'b-', lw=0.5);
        plt.plot(t, CC[:,1], 'r-', lw=0.5);
        plt.plot(t, CC[:,2], 'g-', lw=0.5);
        plt.grid(True)

        plt.show()


        


                       
   
#================================================================
    
#================================================================
    
def ReadParam(fbase):
  
    fname = fbase + '/'+ fbase + '.param'

    
    f = open(fname, 'r')
    s = f.read()
    f.close()
    
    ss = s.split("\n\n")

    for s in ss:
        
        s1 = s.replace("\n#", ";")
        s2 = s1.split("\n")        
        
        if ( s2[0].count('name: nsave') > 0):
            nsave = int(s2[1])
        if ( s2[0].count('name: m') > 0):
            m = int(s2[1])

                      
    return([m, nsave])

    
#================================================================


def ReadData(fbase, m, istart = 1,  iend = -1):

    import os
    import numpy as np
 
    i = istart;
    dat = np.empty( shape=(0, 2*m+1) )
    
    while True:

        fname = fbase + '/'+ fbase + '.' + str(i).zfill(4) + '.ak'
        
        if (os.path.isfile(fname) &  ((iend < 0) | (i < iend)) ) :
            dat0 = np.fromfile(fname, 'float64')
            dat0 = dat0.reshape(2*m+1, round(len(dat0)/(2*m+1))).transpose()
            dat = np.vstack((dat, dat0))
            i = i+1
        else:
            #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
            break

    t = dat[:,0]
    
    ire = np.arange(m)*2 + 1
    iim = ire + 1;
    
    phi = np.angle(dat[:,ire] + 1j*dat[:,iim])

    rhosq = dat[:,ire]**2 + dat[:,iim]**2
       
    return(rhosq, phi, t)    


#================================================================

