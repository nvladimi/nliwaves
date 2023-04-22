
#import importlib ; import fiboPost
#importlib.reload(fiboPost) ; fiboPost.Test()

def Test():     #fbase,  fbaseout=''):

    f="s5m10"
    istart=1000
    iend=-1
    
    #modeCorr(f, [2,1], [0], istart=istart, convert=False, fbaseout='testout/'+f+'_21-0');
    #modeCorr(f, [2,2,1,1], [0,0], istart=istart, convert=False, fbaseout='testout/'+f+'_2211-00');
    #modeCorr(f, [2,1,0], [2,1,0], istart=istart, convert=False, fbaseout='testout/'+f+'_210-210');

    pdfFlux(f, binsize=1, numbins=2000, istart=istart, convert=False, fbaseout='testout/'+f);
        
    #Moments(fbase, nmom=6, istart=istart, iend=iend, convert=True, fbaseout=fbaseout)
    #Evol(fbase, istart=istart, iend=iend, convert=True, fbaseout=fbaseout)


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

def AtoB(m,alpha,PI):
   import numpy as np

   p = np.arange(m) + 1
   p = (p*(1+alpha) - 1 - 2*alpha)/3
   #c = ( (1 + np.sqrt(5))/2 )**p * P**(-1/3)
   phi = (1 + np.sqrt(5))/2
   c = phi**p * (PI/2)**(-1/3) * 5**(-alpha/6)
       
   return(c)



#==============================================================================

def Evol(fbase, istart=1, iend=-1, convert=False, fbaseout=''):

    import numpy as np
    import os

    sets = datasets(fbase)
    nmodes, nsave, alpha, PI = ReadParam(sets[0])
    a2b = AtoB(nmodes, alpha, PI)
    a2bsq = a2b**2

     
    ntot = 0
    
    for s in sets:

        ifile = istart
        
        while True:
           
            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'       
            fout  = fbaseout + '/'+ s + '.' + str(ifile).zfill(4) + '.aatot'       
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break
            if np.isnan(np.sum(dat)):
                print("NAN in " + fname)
                continue

            ire = np.arange(nmodes)*2 + 1
            iim = ire + 1;
    
            aak = dat[:,ire]**2 + dat[:,iim]**2

            n = aak.shape[0]
            if convert:
                for i in range(n):
                    aak[i,:] = aak[i,:] * a2bsq
            
            aatot = np.sum(aak, 1)
            
            aatot.tofile(fout)
            #P123 = np.fromfile(fout)
    
#==============================================================================

def Moments(fbase, nmom=12, istart=1, iend=-1, convert=False, fbaseout=''):

    import numpy as np
    import os

    sets = datasets(fbase)
    nmodes, nsave, alpha, PI = ReadParam(sets[0])
    a2b = AtoB(nmodes, alpha, PI)

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
            if np.isnan(np.sum(dat)):
                print("NAN in " + fname)
                continue

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
        h = h + "0.i  1.Fi  2.m1  3.m2 ...  \n\n"
        
        fnameout = fbaseout + "_mts.txt"
        np.savetxt(fnameout, dataout, header = h)
  

#==============================================================================
#
# Compute "n_k", multmode correlator
# Q[i] = < a[i-jRe] ... a[i] ... conj(a[i-jIm])  ... conj(a[i])>
# and nj[i] = < n[i-jRe] ... n[i] ... n[i-jIm]  ... n[i])>
# where "jRe" and "jIm" are taken from an input arrays.
# for flux J use jRe = [2, 1, 0], jIm = [0, 0, 1]
# Output is text file with columns:  0.i  1.nk  2.sqrt(nj)   3.Re(Q)  4.Im(Q)


def modeCorr(fbase, jRe, jIm, istart=1, iend=-1, convert=False, fbaseout=''):
    
    import numpy as np
    import os

    sets = datasets(fbase)

    print(sets)
    nmodes, nsave, alpha, PI = ReadParam(sets[0])
    
    a2b = AtoB(nmodes, alpha, PI)

    #-- avearaging --

    Q  = np.zeros(nmodes) + 1j *  np.zeros(nmodes)
    R = np.zeros(nmodes)
    nk = np.zeros(nmodes)
    ntot = 0
    
    for s in sets:
        
        ifile = istart
        
        while True:

            #-- read data --
            
            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'       

            #print(fname)
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break

            ire = np.arange(nmodes)*2 + 1
            iim = ire + 1;
    
            ak = dat[:,ire] + 1j*dat[:,iim]
            n = ak.shape[0]
            if convert:
                for i in range(n):
                    ak[i,:] = ak[i,:] * a2b

            #-- occupation numbers --

            q   = np.real(ak * np.conjugate(ak))
            nk += np.sum(q, 0)

            #-- correlator --
            
            q = np.ones([n, nmodes])
        
            for j in jRe:
                if (j>0):
                    q = q * np.hstack(( np.zeros([n,j]), ak[:,:-j] ))
                else:
                    q = q * ak
                    
            for j in jIm:
                if (j>0):
                    q = q * np.conj(np.hstack(( np.zeros([n,j]), ak[:,:-j] )))
                else:
                    q = q * np.conj(ak)
                    
            Q  += np.sum(q,0)
            R  += np.sum(np.abs(q*np.conj(q)), 0)

            ntot += n

    nk  = nk/ntot
    R   = R/ntot
    Q   = Q/ntot
    Q   = Q * (1j)**(len(jRe) - len(jIm))

    #-- normalization factor --
  
    nj = np.ones(nmodes)
    
    for j in jRe:     
        if (j>0):
            nj = nj * np.hstack(( np.zeros(j), nk[:-j] ))
        else:
            nj = nj * nk
    for j in jIm:     
        if (j>0):
            nj = nj * np.hstack(( np.zeros(j), nk[:-j] ))
        else:
            nj = nj * nk

    #-- Output to text file --
    
    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_mcr.txt"

        i = np.arange(nmodes)+1
        
        dataout = np.vstack( (i, nk, np.sqrt(nj), np.real(Q), np.imag(Q),  np.sqrt(R) ) ).T
        
        
        h = 'Run \"'+ fbase + '\':  ' + str(ntot) + ' snapshots\n'
        h = h + "jRe = {}\n".format(jRe)  + "jIm = {}\n".format(jIm) 
        h = h + "0.i  1.nk  2.sqrt|nj|   3.Re(Q)  4.Im(Q)   5.R \n\n"
        
        np.savetxt(fnameout, dataout, header = h)

#==============================================================================
#
# Compute pdf of flux J = a[i-2] a[i-1] conj(a)

def pdfFlux(fbase, binsize=0.1, numbins=100, istart=1, iend=-1, convert=False, fbaseout=''):
    
    import numpy as np
    import os

    sets = datasets(fbase)

    print(sets)
    nmodes, nsave, alpha, PI = ReadParam(sets[0])
    
    a2b = AtoB(nmodes, alpha, PI)

    #-- bins --

    #phi = (1 + np.sqrt(5))/2
    #Javg = phi**(-3/2)

    edg = (np.arange(-numbins/2,numbins/2+1)) * binsize;
    edg[0] = -1.e16; edg[-1] = 1.e16

    Pre = np.zeros((numbins, nmodes), 'int32')
    Pim = np.zeros((numbins, nmodes), 'int32')

    
    #-- count data --
    
    ntot = 0
    
    for s in sets:
        
        ifile = istart
        
        while True:

            #-- read data --
            
            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'       

            #print(fname)
            
            if (os.path.isfile(fname) &  ((iend < 0) | (ifile < iend)) ) :
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, round(len(dat)/(2*nmodes+1))).transpose()
                ifile += 1
            else:
                #print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
                break

            ire = np.arange(nmodes)*2 + 1
            iim = ire + 1;
    
            ak = dat[:,ire] + 1j*dat[:,iim]
            n = ak.shape[0]
            if convert:
                for i in range(n):
                    ak[i,:] = ak[i,:] * a2b

            #-- flux --

            J = 1j*np.conj(ak)
            j=1;  J = J * np.hstack(( np.zeros([n,j]), ak[:,:-j] ))
            j=2;  J = J * np.hstack(( np.zeros([n,j]), ak[:,:-j] ))
           
            #-- histogram ---

            for i in range(0,nmodes):
                (h, edges) = np.histogram(np.real(J[:,i]), edg)                
                Pre[:,i] += h.astype('int32')
                (h, edges) = np.histogram(np.imag(J[:,i]), edg)
                Pim[:,i] += h.astype('int32')
 
            ntot += n
            
    
    #-- Output to text file --
    
    if (len(fbaseout) > 0):

        fnameout = fbaseout + "_Jpdf.param"
       
        f = open(fnameout, 'w')
        
        f.write("\n")

        f.write("fbase       " + fbase + "\n")
        f.write("istart      " + str(istart) + "\n")
        f.write("iend        " + str(iend) + "\n")
        f.write("ntot        " + str(ntot)   + "\n")
      
        f.write("\n")

        f.write("nmodes      " + str(nmodes)     + "\n")
        f.write("binsize     " + str(binsize)    + "\n")
        f.write("numbins     " + str(numbins)    + "\n")

        f.write("\n")

        f.close()

        Pim.tofile(fbaseout + "_JPim.dat")
        Pre.tofile(fbaseout + "_JPre.dat")
 
        
              
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
        
        if ( s2[0].count('name: nsave;') > 0):
            nsave = int(s2[1])
        if ( s2[0].count('name: m;') > 0):
            m = int(s2[1])

        if ( s2[0].count('name: alpha;') > 0):
            alpha = float(s2[1])

        if ( s2[0].count('name: PI;') > 0):
            PI = float(s2[1])

            
                      
    return([m, nsave, alpha, PI])


 
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


#==============================================================================

def avgAll(dir, fnameout=''):

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

