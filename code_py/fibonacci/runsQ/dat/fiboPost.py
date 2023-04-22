
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

    [m, nsave, v] = ReadParam(sets[0])

    return(m, nsave, sets, v)
    

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

def Moments(fbase, nmom=12, istart=1, iend=-1, fbaseout=''):

    import numpy as np
    import os

    nmodes, nsave, sets = datasets(fbase)
    
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

def Flux(fbase, istart=1, iend=-1, fnameout=''):

    import numpy as np
    import os

    m, nsave, sets, v = datasets(fbase)

    Fi = Fibonacci(m)
    
    if v == 0:
        p  = m-11
        r  = p-3
    elif v == 1:
        p  = 9
        r  = p+3
    else:
        print("alpha = {}".format(v))
        
 
    #-- avearaging --

    Q  = np.zeros(m) + 1j *  np.zeros(m)
    Qp = np.zeros(m) + 1j *  np.zeros(m)
    Qr = np.zeros(m) + 1j *  np.zeros(m)
    nk = np.zeros(m)

    ntot = 0
    
    for s in sets:
        
        ifile = istart
        
        while True:
           
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

            q   = np.real(ak * np.conjugate(ak))
            nk += np.average(q, 0)

            q   = ak[:,:-2] * ak[:,1:-1] * np.conjugate(ak[:,2:])
            Q  += np.hstack(([0,0], np.average(q, 0)))

            if v == 0:
                
                a  = np.tile(ak[:,p], (ak.shape[1]-2, 1)).T
                q  = ak[:,:-2] * ak[:,1:-1] * np.conjugate(a)
                Qp +=  np.hstack(([0,0], np.average(q, 0)))
                
                a  = np.tile(ak[:,r], (ak.shape[1]-2, 1)).T
                q  = ak[:,:-2] * ak[:,1:-1] * np.conjugate(a)
                Qr +=  np.hstack(([0,0], np.average(q, 0)))

            elif v == 1:

                a  = np.tile(ak[:,p], (ak.shape[1]-2, 1)).T
                q  = a * ak[:,1:-1] * np.conjugate(ak[:,2:])
                Qp +=  np.hstack(([0,0], np.average(q, 0)))

                a  = np.tile(ak[:,r], (ak.shape[1]-2, 1)).T
                q  = a * ak[:,1:-1] * np.conjugate(ak[:,2:])
                Qr +=  np.hstack(([0,0], np.average(q, 0)))
                             
            ntot += 1

    nk  = nk/ntot
    Q   = Q/ntot
    Qp  = Qp/ntot
    Qr  = Qr/ntot

    #-- fluxes --

    J = np.imag(Q)
    K = np.real(Q)
    V = Fi**v
    Vm = np.hstack((0, V[:-2])) 
    
    PI = 2*np.hstack((Fi[1:] * (Vm * J[1:]), 0)) + 2*np.hstack((Fi[:-2] * (V[:-2] * J[2:]), [0,0]) )

    #-- correlators with pumping --

    Q = Q / np.hstack(( [1,1], np.sqrt(nk[:-2]*nk[1:-1]*nk[2:]) ))

    ap = np.sqrt(nk[p])
    ar = np.sqrt(nk[r])
    
    if v == 0:
        Qp  = Qp/np.hstack(( [1,1],  ap * np.sqrt(nk[:-2] * nk[1:-1])  ))
        Qr  = Qr/np.hstack(( [1,1],  ar * np.sqrt(nk[:-2] * nk[1:-1])  ))
        
    elif v == 1:
        
        Qp  = Qp/np.hstack(( [1,1],  ap * np.sqrt(nk[1:-1] *nk[2:])  ))
        Qr  = Qr/np.hstack(( [1,1],  ar * np.sqrt(nk[1:-1] *nk[2:])  ))
  
   
    #-- Output to text file --
    
    if (len(fnameout) > 0):

        i = np.arange(m)+1
        
        dataout = np.vstack( (i, Fi, nk, PI, J, K,
                              np.real(Q), np.imag(Q),
                              np.real(Qp), np.imag(Qp),
                              np.real(Qr), np.imag(Qr),
        ) ).T
        
        
        h = 'Run \"'+ fbase + '\':  ' + str(ntot) + ' snapshots\n'
        h = h + "p = {},   r = {}\n".format(p, r)
        h = h + "1.i  2.Fi  3.nk  4.PI  5.J  6.K "
        h = h + "7.Re(Q) 8.Im(Q) 9.Re(Qp) 10.Im(Qp) 11.Re(Qr) 12.Im(Qr)\n\n"
        
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
            v = int(s2[1])

                      
    return([m, nsave, v])


 
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

