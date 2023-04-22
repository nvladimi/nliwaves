
#import importlib ; import fiboKMI
#importlib.reload(fiboKMI ; fiboKMI.Test('q5_v0m30', iend=4, fbaseout='tmp')

def Test(fbase, istart=1, iend=-1, fbaseout=''):

    i = 25

    probability3D(i, fbase, istart, iend, binsize=0.2, numbins=11, fbaseout=fbaseout)

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


def probability3D(imode, fbase, istart=1, iend=-1, binsize=0.1, numbins=200, fbaseout=''):
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
    P   = np.zeros((numbins, numbins, numbins), 'int32')

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
            rho = rho[:, imode-2:imode+2]
            
            s1 = np.log(rho[:,1]/rho[:,0])
            s2 = np.log(rho[:,2]/rho[:,1])
            s3 = np.log(rho[:,3]/rho[:,2])
           
            ntot += rho.shape[0]
 
            u = np.array([s1, s2, s3]).T      
            (h, edges) = np.histogramdd(u, (edg, edg, edg))
            
            P += h.astype('int32')

    ncount = np.sum(P)
    
    #print("ntot = ", ntot,  ",   ncount = ",  ncount )


    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_kmi.param"
       
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
 
        f.write("\n")

        f.close()

        P.tofile(fbaseout + "_kP123.dat")
 
  

#================================================================

def MI(fbase):
    
    import numpy as np

    pi = np.pi
    
    imode, numbins, bs, ntot =  ReadParamMI(fbase)
   
    P123  = np.fromfile(fbase + "_kP123.dat",  'int32')

    #-- normalize and reshape probabilities --
    
    ncount = P123.sum()

    #print (ntot, ncount)
    
    P123  = P123 / (ncount * bs*bs*bs)
    P123  = P123.reshape(numbins, numbins, numbins)
    
    P12   = P123.sum(2)*bs
    P23   = P123.sum(0)*bs
    P31   = P123.sum(1)*bs
       
    P1    = P123.sum((1,2))*bs*bs
    P2    = P123.sum((0,2))*bs*bs
    P3    = P123.sum((0,1))*bs*bs


    #-- Entropies --
    
    ind = (P123>0);  S123 = -np.sum(P123[ind]* np.log2(P123[ind])) *bs*bs*bs

    ind = (P12>0);   S12  = -np.sum(P12[ind] * np.log2(P12[ind])) *bs*bs
    ind = (P23>0);   S23  = -np.sum(P23[ind] * np.log2(P23[ind])) *bs*bs
    ind = (P31>0);   S31  = -np.sum(P31[ind] * np.log2(P31[ind])) *bs*bs

    ind = (P1>0);    S1  = -np.sum(P1[ind] * np.log2(P1[ind])) * bs
    ind = (P2>0);    S2  = -np.sum(P2[ind] * np.log2(P2[ind])) * bs
    ind = (P3>0);    S3  = -np.sum(P3[ind] * np.log2(P3[ind])) * bs

    
    #-- Mutual information -- 

    I123 = S1 + S2 + S3 - S123
    
    I12  = S1 + S2 - S12 
    I23  = S2 + S3 - S23
    I31  = S3 + S1 - S31 

    II   = I12 + I23 + I31 - I123


    #-- Output to text file --
 
 
    f = open(fbase + "_kmi.txt", 'w')
 
    f.write("\n") 

    f.write("# 1.i  2.S123   3.S12  4.S23  5.S31  6.S1  7.S2  8.S3\n\n")  

    dat = (imode,  S123, S12, S23, S31, S1, S2, S3)

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
  
    fname = fbase + '_kmi.param'

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

