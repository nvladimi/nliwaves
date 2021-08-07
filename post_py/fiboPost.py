
#import importlib ; import fiboPost
#importlib.reload(fiboPost) ; fiboPost.Test('q5_v0m30', iend=3)

def Test(fbase, alpha, istart=1, iend=-1, fbaseout=''):


    CorrPumping(fbase, alpha, istart=istart, iend=iend, fnameout=fbaseout)
    #Moments(fbase, nmom=6, istart=istart, iend=iend, fbaseout=fbaseout)

    #avgSpectra("POSTspc/q7", fnameout='')

#==============================================================================
    
def datasets(fbase):
    
    import numpy as np
    import os
    
    sets = list()
    
    for file in os.listdir("."):
        if file.startswith(fbase):
            sets.append(file)
            #print(os.path.join("/mydir", file))

    return(sets)

#==============================================================================

def AtoB(m,v,P):
   import numpy as np

   p = np.arange(m) + 1
   p = (p*(1+v) - 1 - 2*v)/3
   #c = ( (1 + np.sqrt(5))/2 )**p * P**(-1/3)
   phi = (1 + np.sqrt(5))/2
   c = phi**p * (P/2)**(-1/3) * 5**(-v/6)
       
   return(c)

#==============================================================================

def avgSpectra(dir, fnameout=''):

    import numpy as np
    import os

    count = 0
    
    for subdir in os.listdir(dir):
        for file in os.listdir(dir + '/' + subdir):
            
            fname = dir + '/' + subdir + '/' + file 
            dat = np.loadtxt(fname)

            if 'alldat' in locals():
                alldat += dat
            else:
                alldat = dat

            count += 1

    dataout = alldat/count

    h = "Averaging " + str(count) + " files in " + dir + " by fiboPost.py\n\n"
    
    np.savetxt(fnameout, dataout, header = h)


    
#==============================================================================

def Moments(fbase, nmom=12, istart=1, iend=-1, convert=False, fbaseout=''):

    import numpy as np
    import os

    sets = datasets(fbase)
    nmodes, nsave, v, P = ReadParam(sets[0])
    a2b = AtoB(nmodes, v, P)

    momsum = np.zeros((nmodes, nmom))
     
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
    
            ak = np.sqrt(dat[:,ire]**2 + dat[:,iim]**2)

            n = ak.shape[0]
            if convert:
                for i in range(n):
                    ak[i,:] = ak[i,:] * a2b
            
            for k in range(1,nmom+1):
    
                momsum[:,k-1] += np.average(ak**k,0)
                
            ntot += 1
        
    momavg = momsum/ntot

    #-- Output to text file --

    Fi = Fibonacci(nmodes)
    
    imodes = np.arange(nmodes)+1

    if (len(fbaseout) > 0):
        
        dataout = np.vstack((imodes, Fi, momavg.T)).T
       
        h = "Run \"" + fbase + "\":  " + str(ntot*nsave) + " elements\n\n"
        h = h + "1.i  2.Fi  3.m1  4.m2 ...  \n\n"
        
        fnameout = fbaseout + "_mts.txt"
        np.savetxt(fnameout, dataout, header = h)
  
   
    if False:
        import matplotlib.pyplot as plt

        a2  =  momavg[:,1]
        a4  =  momavg[:,3]
        a6  =  momavg[:,5]
              
        plt.figure(figsize=(12,6))

        plt.subplot(121)
        plt.yscale('log')
        
        plt.plot(imodes, 22*Fi**(1/3), 'bx-');
        plt.plot(imodes, Fi*a2, 'ro-');
        plt.grid(True)

        plt.subplot(122)
        plt.yscale('log')
        
        plt.plot(imodes, a4 / (2*a2**2), 'bo-');
        plt.plot(imodes, a6 / (6*a2**3), 'ro-');
        plt.grid(True)

        plt.show()


#==============================================================================
#
# Compute "n_k" and multmode correlator Q[k] = <a^*[k] a[k+j1] a[k+j2] ... >
# where positive or negative "j" are taken from an input array "jj" 
# for flux J use jj = [-2, -1]
# for K use jj = [-4, -3, -1]
# for H use jj = [-6, -5, -3, -1] 

def modeCorr(fbase, jj, istart=1, iend=-1, convert=True, fbaseout=''):
    
    import numpy as np
    import os

    sets = datasets(fbase)
    m, nsave, v, P = ReadParam(sets[0])
    
    a2b = AtoB(m, v, P)

    #-- avearaging --

    Q  = np.zeros(m) + 1j *  np.zeros(m)
    nk = np.zeros(m)
    ntot = 0
    
    for s in sets:
        
        ifile = istart
        
        while True:

            #-- read data --
            
            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'       

            print(fname)
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*m+1, round(len(dat)/(2*m+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break

            ire = np.arange(m)*2 + 1
            iim = ire + 1;
    
            ak = dat[:,ire] + 1j*dat[:,iim]
            n = ak.shape[0]
            if convert:
                for i in range(n):
                    ak[i,:] = ak[i,:] * a2b

            #-- correlator --

            q   = np.real(ak * np.conjugate(ak))
            nk += np.sum(q, 0)
                    
            q = np.conjugate(ak)
            for j in jj:
                if (j<0):
                    q = q * np.hstack((np.zeros([n,-j]), ak[:,:j]))
                elif (j>0):
                    q = q * np.hstack(ak[j:,:], np.zeros([n,j]))
                else:
                    q = q * ak
                    
            Q  += np.sum(q,0)

            ntot += n

    nk  = nk/ntot
    Q   = Q/ntot

    #-- normalization factor --

    q = nk
    for j in jj:
         
        if (j<0):
            q = q * np.hstack(( np.zeros(-j), nk[:j] ))
        elif (j>0):
            q = q * np.hstack(( nk[j:], np.zeros(j) ))
        else:
            q = q * nk

   
    #-- Output to text file --
    
    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_mcr.txt"

        i = np.arange(m)+1
        
        dataout = np.vstack( (i, nk, np.sqrt(q), np.real(Q), np.imag(Q) ) ).T
        
        
        h = 'Run \"'+ fbase + '\':  ' + str(ntot) + ' snapshots\n'
        h = h + "j = {}\n".format(jj)
        h = h + "1.i  2.nk  3.sqrt|nj|   4.Re(Q)  5.Im(Q) \n\n"
        
        np.savetxt(fnameout, dataout, header = h)

        
        
#================================================================

def Fibonacci(n):

    import numpy as np
     
    Fi = np.zeros(n, 'uint64')
     
    Fi[0]=1
    Fi[1]=1

    for i in range(2,n):
        
        Fi[i] = Fi[i-1] + Fi[i-2]
        
    return(Fi)
  
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

        if ( s2[0].count('name: v') > 0):
            v = float(s2[1])

        if ( s2[0].count('name: P') > 0):
            P = float(s2[1])

            
                      
    return([m, nsave, v, P])


 
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

