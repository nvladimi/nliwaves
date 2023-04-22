
#import importlib ; import fiboKmult
#importlib.reload(fiboKmult) ; fiboKmult.Test('q5_v0m30', iend=4, fbaseout='tmp')

def Test(fbase, istart=1, iend=-1, fbaseout=''):

    i = 25

    probability2D(i, fbase, istart, iend, binsize=0.2, numbins=11, fbaseout=fbaseout)

    #MI(fbaseout)
    plot2Dprob(fbaseout)

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


def probability2D(imode, fbase, istart=1, iend=-1, binsize=0.1, numbins=200, fbaseout=''):
#
# 4D probability for modes [i-2,i-1,i]
    
    import numpy as np
    import os
    
    pi = np.pi

    nmodes, nsave, sets = datasets(fbase)

    numbins1=numbins
    numbins2=numbins
    
    edg1 = (np.arange(-numbins1/2, numbins1/2+1))*binsize;
    edg2 = (np.arange(-numbins2/2, numbins2/2+1))*binsize;
    cnt1 = (np.arange(-numbins1/2, numbins1/2)+0.5)*binsize;
    cnt2 = (np.arange(-numbins2/2, numbins2/2)+0.5)*binsize;
    edg1[0]  = -1.e16
    edg2[0]  = -1.e16
    edg1[-1] =  1.e16
    edg2[-1] =  1.e16

    
    P   = np.zeros((numbins1, numbins2), 'int32')

    smmin =  1.e+16
    smmax = -1.e+16
    smavg =  0.0
    spmin =  1.e+16
    spmax = -1.e+16
    spavg =  0.0
    ntot  = 0
     
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

            ire = np.arange(nmodes)*2 + 1
            iim = ire + 1;
    
            rho = np.sqrt(dat[:,ire]**2 + dat[:,iim]**2)
            rho = rho[:, imode-2:imode+1]
            
            sm = np.log(rho[:,1]/rho[:,0])
            sp = np.log(rho[:,2]/rho[:,1])

            smmin = min(smmin, sm.min())
            smmax = max(smmax, sm.max())
            smavg += np.average(sm)
            
            spmin = min(spmin, sp.min())
            spmax = max(spmax, sp.max())
            spavg += np.average(sp)

            ntot += rho.shape[0]
 
            u = np.array([sm, sp]).T      
            (h, edges) = np.histogramdd(u, (edg1, edg2))
            
            P += h.astype('int32')

    ncount = np.sum(P)
    
    #print("ntot = ", ntot,  ",   ncount = ",  ncount )


    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_Kmult.param"
       
        f = open(fnameout, 'w')
        
        f.write("\n")

        f.write("fbase       " + fbase + "\n")
        f.write("istart      " + str(istart) + "\n")
        f.write("iend        " + str(iend) + "\n")
       
        f.write("\n")

        f.write("nmodes      " + str(nmodes)     + "\n")
        f.write("imode       " + str(imode)      + "\n")
        f.write("binsize     " + str(binsize)    + "\n")
        f.write("numbins     " + str(numbins)    + "\n")
 
        f.write("\n")

        f.write("ntot        " + str(ntot)   + "\n")
        f.write("ncount      " + str(ncount) + "\n")
        f.write("sigma_min   " + str((smmin, spmin))  + "\n")
        f.write("sigma_max   " + str((smmax, spmax))  + "\n")
        f.write("sigma_avg   " + str((smavg, spavg))  + "\n")

        f.write("\n")

        f.close()

        P.tofile(fbaseout + "_Pkmult.dat")
 
  

#================================================================

def MI(fbase):
    
    import numpy as np

    pi = np.pi
    
    imode, numbins, binsize, ncount =  ReadParamMI(fbase)
   
    P12  = np.fromfile(fbase + "_Pkmult.dat",  'int32')

    #-- normalize and reshape probabilities --
    
    ntot = P12.sum()

    #print (ntot, ncount)
    
    P12  = P12 / (ncount * binsize* binsize)
    P12  = P12.reshape(numbins, numbins)
    P1   = P12.sum(1)*binsize
    P2   = P12.sum(0)*binsize


    #-- Entropies --
    
    ind = (P12>0);   S12 = -np.sum(P12[ind] *np.log2(P12[ind]) ) *binsize*binsize
    ind = (P1>0);    S1  = -np.sum(P1[ind] * np.log2(P1[ind])) * binsize
    ind = (P2>0);    S2  = -np.sum(P2[ind] * np.log2(P2[ind])) * binsize
 
    #print(P12.shape, P1.shape, P2.shape)

    
    #-- Mutual information -- 

    I12 = S1 + S2 - S12


    #-- Output to text file --

    cnt = (np.arange(-numbins/2, numbins/2)+0.5)*binsize;
    edg = (np.arange(-numbins/2, numbins/2+1))*binsize;

    dataout = np.vstack((cnt, P1,  P2)).transpose()

    h = "Header:  1.i  2.I12   3.S12  4.S1  5.S2\n"
    h = h + '{:2d}  {:e}  {:e}  {:e}  {:e}\n\n'.format(imode, I12, S12, S1, S2)
    h = h + 'Data: 1.sigma  2.P(sigma_n)  3.P(sigma_n+1)\n\n'
    

    fnameout = fbase + "_Kmult.txt"
    np.savetxt(fnameout, dataout, header = h)


    #-- plot 1D distributions of probabilities ---

    
    if True:
        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6,4))

        plt.subplot()
        plt.yscale('log')
        plt.xticks(edg)
      
        plt.plot(cnt,  P1, 'bo-');
        plt.plot(cnt,  P2, 'ro-');
        plt.grid(True)
      
        plt.show()

    return
    
#================================================================

def plot2Dprob(fbase, alpha=-1):
    
    import numpy as np

    pi = np.pi

    equil = {0: -0.1604039416865561,    # -1/3 ln (phi)
             1: -0.3208078833731123,    # -2/3 ln (phi)
             -1: 0}

    #lev0 = {0: -0.55,   1: -0.50,   -1: 0}

    
    imode, numbins, binsize, ncount =  ReadParamMI(fbase)
    numbins1=numbins
    numbins2=numbins
 

    
    P12  = np.fromfile(fbase + "_Pkmult.dat",  'int32')

    #-- normalize and reshape probabilities --
    
    ntot = P12.sum()

    #print (ntot, ncount)
    
    P12  = P12 / (ncount * binsize* binsize)
    P12  = P12.reshape(numbins1, numbins2)
    P1   = P12.sum(1)*binsize
    P2   = P12.sum(0)*binsize

    #-- plot 2D distributions of probabilities ---

    levels = np.hstack((np.arange(-6,0,1), -0.55))
    
    if True:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6,6))
        plt.xlim(-6,6)
        plt.ylim(-6,6)
        plt.grid(True)


        
        logP12 =np.log10(P12)
        print(np.min(logP12), np.max(logP12))
        
        plt.imshow(logP12, cmap='jet', origin='lower');

    
        x=np.arange(-10,11,1)
        plt.plot(x,-x, '-k', linewidth=0.5)
        plt.contour(logP12, levels, cmap='jet',
                    extent=(-10 - equil[alpha], 10 - equil[alpha],
                            -10 - equil[alpha], 10 - equil[alpha]))
        plt.title(fbase.replace("_", " "))

        plt.savefig(fbase + ".pdf", pad_inches=0)
        
        plt.show()
        
        plt.close()
        
    return
    


  
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
    
def ReadParamMI(fbase):
  
    fname = fbase + '_Kmult.param'

    f = open(fname, 'r')
    s = f.read()
    f.close()
    
    ss = s.split("\n")

    for s in ss:
 
        s2 = s.split()        

        if (len(s2)>0):

            if (s2[0] == 'ncount'):
                ncount= int(s2[1])
            if (s2[0] == 'imode'):
                imode = int(s2[1])
            if (s2[0] == 'numbins'):
                numbins = int(s2[1])
            if (s2[0] == 'binsize'):
                binsize = float(s2[1])


    return(imode, numbins, binsize, ncount)
    
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

