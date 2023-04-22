
#import importlib ; import fiboMI
#importlib.reload(fiboMI) ; fiboMI.Test('q5_v0m30', fbaseout='mitest/q5_i10')

def Test(fbase, istart=1, iend=-1, fbaseout=''):

    i = 10

    #probability4D(i, fbase, istart, iend, binsize=1, numbins=32, numbinsPhi=16, fbaseout=fbaseout)

    #MI(fbaseout)
    probability1D(fbaseout)

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

def avgOccupation(nmodes, sets, istart, iend):

    import numpy as np
    import os
    
    nsum = np.zeros(nmodes)
    ntot = 0
    for s in sets:

        ifile = istart
        
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
    
            phi = np.angle(dat[:,ire] + 1j*dat[:,iim])
            rhosq = dat[:,ire]**2 + dat[:,iim]**2
        
            nsum += np.sum(rhosq,0)
            ntot += rhosq.shape[0]
        
    navg = nsum/ntot

    return(navg, ntot)

#==============================================================================


def probability4D(imode, fbase, istart=1, iend=-1, binsize=1, numbins=100, numbinsPhi=32, fbaseout=''):
#
# 4D probability for modes [i-2,i-1,i]
    
    import numpy as np
    import time
    import os

    time_read = 0
    time_phase = 0
    time_pdf = 0
    time_other = 0
    

    tic = time.perf_counter()
    
    pi = np.pi

    nmodes, nsave, sets = datasets(fbase)

    navg, ntot = avgOccupation(nmodes, sets, istart, iend)

    time_read += time.perf_counter() - tic

    tic = time.perf_counter()

    i = imode - 1
    navg1 = navg[i-2]
    navg2 = navg[i-1]
    navg3 = navg[i]
    
    edg1 = np.arange(numbins+1)*binsize*navg1;  edg1[-1] = 1.e16
    edg2 = np.arange(numbins+1)*binsize*navg2;  edg2[-1] = 1.e16 
    edg3 = np.arange(numbins+1)*binsize*navg3;  edg3[-1] = 1.e16
    edgPhi =  (np.arange(-numbinsPhi/2-1, numbinsPhi/2+1) + 0.5 )/numbinsPhi*2*pi
    
    edges123 = (edg1, edg2, edg3, edgPhi)
    edges12  = (edg1, edg2, edgPhi)
    edges23  = (edg2, edg3, edgPhi)
    edges31  = (edg3, edg1, edgPhi)

    P123 = np.zeros((numbins, numbins, numbins, numbinsPhi+1), 'int32')
    P12  = np.zeros((numbins, numbins, numbinsPhi+1), 'int32')
    P23  = np.zeros((numbins, numbins, numbinsPhi+1), 'int32')
    P31  = np.zeros((numbins, numbins, numbinsPhi+1), 'int32')

    time_other = time.perf_counter() - tic
     
    for s in sets:

        ifile = istart;
        
        while True:

            tic = time.perf_counter()
           
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
    
            phi = np.angle(dat[:,ire] + 1j*dat[:,iim])
            rhosq = dat[:,ire]**2 + dat[:,iim]**2

            time_read += time.perf_counter() - tic

            tic = time.perf_counter()
        
            i = imode - 1
        
            n1 = rhosq[:, i-2]
            n2 = rhosq[:, i-1]
            n3 = rhosq[:, i]
        
            f123 = -phi[:,i-2] - phi[:,i-1] + phi[:,i] 
            f12  =             - phi[:,i-1] + phi[:,i] 
            f23  = -phi[:,i-2] + phi[:,i-1] 
            f31  =  phi[:,i-2]              - phi[:,i] 
        
            ind = (f123 >  pi);  f123[ind] += -2*pi;
            ind = (f123 >  pi);  f123[ind] += -2*pi;
            ind = (f123 < -pi);  f123[ind] +=  2*pi;
            ind = (f123 < -pi);  f123[ind] +=  2*pi;

            ind = (f12 >  pi);   f12[ind] += -2*pi;
            ind = (f12 < -pi);   f12[ind] +=  2*pi;
        
            ind = (f23 >  pi);   f23[ind] += -2*pi;
            ind = (f23 < -pi);   f23[ind] +=  2*pi;
            
            ind = (f31 >  pi);   f31[ind] += -2*pi;
            ind = (f31 < -pi);   f31[ind] +=  2*pi;

            time_phase +=  time.perf_counter() - tic
        
            tic = time.perf_counter()
        
            u = np.array([n1, n2, n3, f123]).T      
            (h, edges) = np.histogramdd(u, edges123)
            P123 += h.astype('int32')

            u = np.array([n1, n2, f12]).T      
            (h, edges) = np.histogramdd(u, edges12)
            P12 += h.astype('int32')
            
            u = np.array([n2, n3, f23]).T      
            (h, edges) = np.histogramdd(u, edges23)
            P23 += h.astype('int32')

            u = np.array([n3, n1, f31]).T      
            (h, edges) = np.histogramdd(u, edges31)
            P31 += h.astype('int32')

            time_pdf += time.perf_counter() - tic

    #print("ntot = ", ntot,  ",   ncount = ",  np.sum(P123) ) 

    tic += time.perf_counter()
    
    P123[:,:,:,0] +=  P123[:,:,:,-1]
    P12[:,:,0]  +=  P12[:,:,-1]
    P23[:,:,0]  +=  P23[:,:,-1]
    P31[:,:,0]  +=  P31[:,:,-1]

    P123 = P123[:,:,:,:-1]
    P12  = P12[:,:,:-1]
    P23  = P23[:,:,:-1]
    P31  = P31[:,:,:-1]
    
    ncount = np.sum(P123)
    
    print("ntot = ", ntot,  ",   ncount = ",  ncount )

    print(f"Time read:  {time_read:0.4f} seconds")
    print(f"Time phase: {time_phase:0.4f} seconds")
    print(f"Time pdf:   {time_pdf:0.4f} seconds")


    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_mi.param"
       
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
        f.write("numbinsPhi  " + str(numbinsPhi) + "\n")

        f.write("\n")

        f.write("ntot        " + str(ntot)   + "\n")
        f.write("ncount      " + str(ncount) + "\n")
        f.write("navg1       " + str(navg1)  + "\n")
        f.write("navg2       " + str(navg2)  + "\n")
        f.write("navg3       " + str(navg3)  + "\n")

        f.write("\n")

        f.close()

        P123.tofile(fbaseout + "_P123.dat")
        P12.tofile(fbaseout + "_P12.dat")
        P23.tofile(fbaseout + "_P23.dat")
        P31.tofile(fbaseout + "_P31.dat")


        time_other += time.perf_counter() - tic

#================================================================

def MI(fbase):
    
    import numpy as np

    pi = np.pi
    
    imode, numbins, numbinsPhi, binsize, (navg1, navg2, navg3), ncount =  ReadParamMI(fbase)

    dn1  = navg1*binsize
    dn2  = navg2*binsize
    dn3  = navg3*binsize
    dphi = 2*pi/numbinsPhi
    
    P123 = np.fromfile(fbase + "_P123.dat", 'int32')
    P12  = np.fromfile(fbase + "_P12.dat",  'int32')
    P23  = np.fromfile(fbase + "_P23.dat",  'int32')
    P31  = np.fromfile(fbase + "_P31.dat",  'int32')

    #-- normalize and reshape probabilities --
    
    ntot = P123.sum()

    # print (ntot, ncount, P12.sum(), P23.sum(), P31.sum() )
    
    P123 = P123 / (ncount *dn1*dn2*dn3*dphi)
    P12  = P12  / (ncount *dn1*dn2    *dphi)
    P23  = P23  / (ncount     *dn2*dn3*dphi)
    P31  = P31  / (ncount *dn1    *dn3*dphi)

    P123 = P123.reshape(numbins, numbins, numbins, numbinsPhi)
    P12  = P12.reshape(numbins, numbins, numbinsPhi)
    P23  = P23.reshape(numbins, numbins, numbinsPhi)
    P31  = P31.reshape(numbins, numbins, numbinsPhi)
 

    #-- integrate 4D probability to obtain 1D probalities of amplitudes -- 
        
    P1   = P123.sum((1,2,3))    *dn2*dn3*dphi
    P2   = P123.sum((0,2,3))*dn1    *dn3*dphi
    P3   = P123.sum((0,1,3))*dn1*dn2    *dphi

    
    #--  averages for staircase distribution --

    n = (np.arange(0,numbins)+0.5)*dn1
    navg1s = sum(P1*n)/sum(P1)

    n = (np.arange(0,numbins)+0.5)*dn2
    navg2s = sum(P2*n)/sum(P2)
    
    n = (np.arange(0,numbins)+0.5)*dn3
    navg3s = sum(P3*n)/sum(P3)

    
    #-- integrate 4D probability to obtain 1D probabiliy of Theta123 --
    
    Pa123 = P123.sum((0,1,2))*dn1*dn2*dn3

    
    #-- integrate 3D probabilities to obtain 1D probabilities of Theta12, Theta23, and Theta31 --
    
    Pa12  = P12.sum((0,1))   *dn1*dn2
    Pa23  = P23.sum((0,1))       *dn2*dn3
    Pa31  = P31.sum((0,1))   *dn1    *dn3


    #-- Entropy of 4D distribution (n1,n2,n3, theta123) --
    
    ind = (P123>0);  S123 = -np.sum(P123[ind]* np.log2(P123[ind])) *dn1*dn2*dn3*dphi

    
    #-- Entropy of 3D distribution (n1,n2, theta12) and like --

    ind = (P12>0);   S12  = -np.sum(P12[ind] * np.log2(P12[ind])) *dn1*dn2    *dphi
    ind = (P23>0);   S23  = -np.sum(P23[ind] * np.log2(P23[ind]))     *dn2*dn3*dphi
    ind = (P31>0);   S31  = -np.sum(P31[ind] * np.log2(P31[ind])) *dn1    *dn3*dphi


    #-- Entropy of 1D distribution for amplitudes --

    ind = (P1>0);    S1  = -np.sum(P1[ind] * np.log2(P1[ind])) * dn1
    ind = (P2>0);    S2  = -np.sum(P2[ind] * np.log2(P2[ind])) * dn2
    ind = (P3>0);    S3  = -np.sum(P3[ind] * np.log2(P3[ind])) * dn3


    #-- Entropy of 1D distribution for angles --
    
    ind = (Pa123>0); Sa123  = -np.sum(Pa123[ind] * np.log2(Pa123[ind])) * dphi
    ind = (Pa12>0);  Sa12   = -np.sum(Pa12[ind]  * np.log2(Pa12[ind]))  * dphi
    ind = (Pa23>0);  Sa23   = -np.sum(Pa23[ind]  * np.log2(Pa23[ind]))  * dphi
    ind = (Pa31>0);  Sa31   = -np.sum(Pa31[ind]  * np.log2(Pa31[ind]))  * dphi

    
    #-- Equilibrium, assuming exponential fit with <n_i> and equidistribution for theta --

    S1o = 1/np.log(2) * (1 + np.log(navg1));
    S2o = 1/np.log(2) * (1 + np.log(navg2));
    S3o = 1/np.log(2) * (1 + np.log(navg3));

    Sao = np.log2(2*pi);

    S1s = 1/np.log(2) * (1 + np.log(navg1s));
    S2s = 1/np.log(2) * (1 + np.log(navg2s));
    S3s = 1/np.log(2) * (1 + np.log(navg3s));
    

    #-- Mutual information -- 

    I123 = S1 + S2 + S3 + Sao - S123
    
    I12  = S1 + S2 + Sao - S12 
    I23  = S2 + S3 + Sao - S23
    I31  = S3 + S1 + Sao - S31 

    II   = I12 + I23 + I31 - I123


    #-- with respect to equilibrium (old approach) --
    
    if False:

        S4D = S123  - (S1o + S2o + S3o + Sa0) ;

        S1 = S1 - S1o;
        S2 = S2 - S2o;
        S3 = S3 - S3o;
        Sa123 = Sa123 - Sao;

        I4D =  S1 + S2 + S3 - S4D;


    #-- Output to text file --
 
 
    f = open(fbase + "_mi.txt", 'w')
 
    f.write("\n") 

    f.write("# 1.i  2.S123   3.S12  4.S23  5.S31  6.S1  7.S2  8.S3\n")  
    f.write("# 9.S1o  10.S2o  11.S3o  12.Sao  13.Sa123  14.Sa12  15.Sa23  16.Sa31\n")  
    f.write("# 17.S1s  18.S2s  19.S3s \n\n")  

    dat = (imode,  S123, S12, S23, S31, S1, S2, S3,
           S1o, S2o, S3o, Sao, Sa123, Sa12, Sa23, Sa31, S1s, S2s, S3s)

    f.write(str(dat) + "\n\n")

    f.write("# 1.i  2.I123   3.I12  4.I23  5.I31  6.II\n\n")  

    dat = (imode, I123, I12, I23, I31, II)  
    
    f.write(str(dat) + "\n\n")
        
    f.close()


    
    #-- plot 1D distributions of probabilities ---

    
    if False:
        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8,4))

        plt.subplot(121)
        plt.yscale('log')

        x = np.arange(numbins) + 0.5
        
        plt.plot(x*dn1,P1, 'bo-');
        plt.plot(x*dn2,P2, 'go-');
        plt.plot(x*dn3,P3, 'ro-');
        plt.grid(True)

        plt.subplot(122)
 
        x = (np.arange(numbinsPhi) + 0.5)/numbinsPhi*2 - 1
        
        plt.plot(x,Pa123*pi, 'ko-');
        plt.plot(x,Pa12*pi,  'ro-');
        plt.plot(x,Pa23*pi,  'bo-');
        plt.plot(x,Pa31*pi,  'go-');

        
        plt.grid(True)
        
        plt.show()

    return
    

#================================================================

def probability1D(fbase):
    
    import numpy as np

    pi = np.pi
    
    imode, numbins, numbinsPhi, binsize, (navg1, navg2, navg3), ncount =  ReadParamMI(fbase)

    dn1  = navg1*binsize
    dn2  = navg2*binsize
    dn3  = navg3*binsize
    dphi = 2*pi/numbinsPhi
    
    P123 = np.fromfile(fbase + "_P123.dat", 'int32')
  
    #-- normalize and reshape probabilities --
    
    ntot = P123.sum()

    # print (ntot, ncount, P12.sum(), P23.sum(), P31.sum() )
    
    P123 = P123 / (ncount *dn1*dn2*dn3*dphi)
    P123 = P123.reshape(numbins, numbins, numbins, numbinsPhi)
 

    #-- integrate 4D probability to obtain 1D probalities of amplitudes -- 
        
    P1   = P123.sum((1,2,3))    *dn2*dn3*dphi
    P2   = P123.sum((0,2,3))*dn1    *dn3*dphi
    P3   = P123.sum((0,1,3))*dn1*dn2    *dphi
    P4   = P123.sum((0,1,2))*dn1*dn2*dn3

    
    #--  averages for staircase distribution --

    bin = np.arange(numbins) + 0.5
    
    n1 = bin*dn1
    navg1s = sum(P1*n1)/sum(P1)

    n2 = bin*dn2
    navg2s = sum(P2*n2)/sum(P2)
    
    n3 = bin*dn3
    navg3s = sum(P3*n3)/sum(P3)


    #-- Output to text file --
            
    dataout = np.vstack((bin*binsize, P1*navg1, P2*navg2, P3*navg3)).T
    
    h = "Run \"" + fbase + "\":  " + str(ntot) + " elements;  "
    h += "modes ({}, {}, {})\n".format(imode-2, imode-1, imode)  
    h += "True averages (1,2,3):        {:e}    {:e}    {:e}\n".format(navg1, navg2, navg3)
    h += "Staircase averages (1,2,3):   {:e}    {:e}    {:e}\n".format(navg1s, navg2s, navg3s) 
    h += "1.n/n_avg   2.P1    3.P2   4.P3\n\n"
        
    fnameout = fbase + "_prob.txt"
    np.savetxt(fnameout, dataout, header = h)
        


    binphi = ( np.arange(numbinsPhi) - numbinsPhi/2 )*dphi
    dataout = np.vstack((binphi, P4)).T
    dataout = np.vstack((dataout, dataout[0,:]))
    dataout[-1,0] = -dataout[-1,0]
    
    h = "Run \"" + fbase + "\":  " + str(ntot) + " elements;  "
    h += "modes ({}, {}, {})\n".format(imode-2, imode-1, imode)  
    h += "1.phi   2.P4\n\n"
        
    fnameout = fbase + "_probphi.txt"
    np.savetxt(fnameout, dataout, header = h)
        

    
    #-- plot 1D distributions of probabilities ---

    
    if False:
        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        import matplotlib.pyplot as plt

        plt.figure(figsize=(4,4))

        plt.subplot(111)
        plt.yscale('log')
        
        plt.plot(n1,P1, 'bo-');
        plt.plot(n2,P2, 'go-');
        plt.plot(n3,P3, 'ro-');
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
        
        if ( s2[0].count('name: nsave;') > 0):
            nsave = int(s2[1])
        if ( s2[0].count('name: m;') > 0):
            m = int(s2[1])

                      
    return([m, nsave])


  
#================================================================
    
def ReadParamMI(fbase):
  
    fname = fbase + '_mi.param'

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
            if (s2[0] == 'numbinsPhi'):
                numbinsPhi = int(s2[1])
            if (s2[0] == 'navg1'):
                navg1 = float(s2[1])
            if (s2[0] == 'navg2'):
                navg2 = float(s2[1])
            if (s2[0] == 'navg3'):
                navg3 = float(s2[1])
            if (s2[0] == 'binsize'):
                binsize = float(s2[1])


    return(imode, numbins, numbinsPhi, binsize, [navg1, navg2, navg3], ncount)
    
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

