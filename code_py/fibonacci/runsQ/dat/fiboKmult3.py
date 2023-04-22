
#import importlib ; import fiboKmult3
#importlib.reload(fiboKmult3) ; fiboKmult3.Test('q5_v0m30', iend=4)

def Test(fbase, istart=1, iend=-1, fbaseout=''):

    i = 25

    probability2D(i, fbase, istart, iend, binsize=0.2, numbins=11, fbaseout=fbaseout)

    MI(fbaseout)

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

    edg = (np.arange(-numbins/2, numbins/2+1))*binsize;
    cnt = (np.arange(-numbins/2, numbins/2)+0.5)*binsize;
    edg[0]  = -1.e16
    edg[-1] =  1.e16
    P   = np.zeros((numbins, numbins), 'int32')

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
            rho = rho[:, imode-2:imode+3]
            
            sm = np.log(rho[:,1]/rho[:,0])
            sp = np.log(rho[:,4]/rho[:,3])

            smmin = min(smmin, sm.min())
            smmax = max(smmax, sm.max())
            smavg += np.average(sm)
            
            spmin = min(spmin, sp.min())
            spmax = max(spmax, sp.max())
            spavg += np.average(sp)

            ntot += rho.shape[0]
 
            u = np.array([sm, sp]).T      
            (h, edges) = np.histogramdd(u, (edg, edg))
            
            P += h.astype('int32')

    ncount = np.sum(P)
    
    #print("ntot = ", ntot,  ",   ncount = ",  ncount )


    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_Kmult3.param"
       
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

        P.tofile(fbaseout + "_Pkmult3.dat")
 
  

#================================================================

def MI(fbase):
    
    import numpy as np

    pi = np.pi
    
    imode, numbins, binsize, ncount =  ReadParamMI(fbase)
   
    P12  = np.fromfile(fbase + "_Pkmult3.dat",  'int32')

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
    

    fnameout = fbase + "_Kmult3.txt"
    np.savetxt(fnameout, dataout, header = h)


    #-- plot 1D distributions of probabilities ---

    
    if False:
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
  
    fname = fbase + '_Kmult3.param'

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

