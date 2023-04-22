
#import importlib ; import fiboPost
#importlib.reload(fiboPost) ; fiboPost.Test()

def Test():     #fbase,  fbaseout=''):

    f="s5m10"
    istart=1000
    iend=-1
    
    #modeCorr(f, [2,1], [0], istart=istart, convert=False, fbaseout='testout/'+f+'_21-0');
    #modeCorr(f, [2,2,1,1], [0,0], istart=istart, convert=False, fbaseout='testout/'+f+'_2211-00');
    #modeCorr(f, [2,1,0], [2,1,0], istart=istart, convert=False, fbaseout='testout/'+f+'_210-210');

    #pdfFlux(f, binsize=1, numbins=2000, istart=istart, convert=False, fbaseout='testout/'+f);
        
    #Moments(fbase, nmom=6, istart=istart, iend=iend, convert=True, fbaseout=fbaseout)
    #Evol(fbase, istart=istart, iend=iend, convert=True, fbaseout=fbaseout)


    #-- corr: text "21-0",  "322-0", "11-30", "222-40", "44333-0", "63333-0", "4111-00"

    tcorr2D(f, "21-0", 30, nt1=16, nt2=10, n=10000, frange=(1000,1010), fbaseout='testout/'+f+'-21-0');
    #tcorr2D(f, "11-30", 30, nt1=10, nt2=10, n=10000, frange=(1000,1010), fbaseout='testout/'+f+'-11-30');

    
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


#--  end modeCorr
   
#==============================================================================

def tcorr2D(fbase, corr, j, nt1=10, nt2=10, n=10000, frange=(0,10), fbaseout=''):
    #
    # corr - text: "21-0",  "322-0", "11-30", "222-40", "44333-0", "63333-0", "4111-00"
    #
    
    import numpy as np
    import os

    sets = datasets(fbase)

    print(sets)
    nmodes, nsave, alpha, PI = ReadParam(sets[0])
    
    #-- avearaging --

    Q  = np.zeros((2*nt1, 2*nt2), dtype=np.complex128)
    ntot = 0
    
    for s in sets:

        A0 = np.zeros( (frange[1]-frange[0])*n,  dtype=np.complex64)
        A1 = np.zeros( (frange[1]-frange[0])*n,  dtype=np.complex64)
        A2 = np.zeros( (frange[1]-frange[0])*n,  dtype=np.complex64)

        #-- collect data --

        for ifile in range(frange[0], frange[1]):
                    
            fname = s + '/'+ s + '.' + str(ifile).zfill(4) + '.ak'       

            #print(fname)
            
            try:
                dat = np.fromfile(fname, 'float64')
                dat = dat.reshape(2*nmodes+1, n).transpose()
            except:
                print("Problem with file read:  " + fname )
                break

            ire = np.arange(nmodes)*2 + 1
            iim = ire + 1;
    
            ak = dat[:,ire] + 1j*dat[:,iim]
            
            i0 = (ifile-frange[0])*n
            A0[i0:i0+n] = ak[:,j]

            if corr == "21-0":
                A1[i0:i0+n] = ak[:,j-1]
                A2[i0:i0+n] = ak[:,j-2]
              
            elif corr == "322-0":
                A1[i0:i0+n] = ak[:,j-2]
                A2[i0:i0+n] = ak[:,j-3]
 
            elif corr == "11-30":
                A1[i0:i0+n] = ak[:,j-1]
                A2[i0:i0+n] = ak[:,j-3]
               
            elif corr == "222-40":
                A1[i0:i0+n] = ak[:,j-2]
                A2[i0:i0+n] = ak[:,j-4]
                
            elif corr == "44333-0":
                A1[i0:i0+n] = ak[:,j-3]
                A2[i0:i0+n] = ak[:,j-4]

            elif corr == "63333-0":
                A1[i0:i0+n] = ak[:,j-3]
                A2[i0:i0+n] = ak[:,j-6]
               
            elif corr =="4111-00":
                A1[i0:i0+n] = ak[:,j-1]
                A2[i0:i0+n] = ak[:,j-4]
                 
            else:
                print("Unknown correlator:  " + corr )
                break

        #-- correlator --

        nt0 = np.maximum(nt1,nt2)
                  
        for n1 in np.arange(-nt1,nt1):
            for n2 in np.arange(-nt2,nt2):

                a0 = A0[ nt0    : -nt0    ]
                a1 = A1[ nt0+n1 : -nt0+n1 ]
                a2 = A2[ nt0+n2 : -nt0+n2 ]

                if corr == "21-0":
                    q = a2 * a1 * np.conjugate(a0)

                elif corr == "322-0":
                    q = a2 * a1 * a1 * np.conjugate(a0)
                   
                elif corr == "11-30":
                    q = a1 * a1 * np.conjugate(a2*a0)
              
                elif corr == "222-40":
                    q = a1 * a1 * a1 * np.conjugate(a2*a0)
                
                elif corr == "44333-0":
                    q = a2 * a2 * a1 * a1 * a1 * np.conjugate(a0)

                elif corr == "63333-0":
                    q = a2 * a1 * a1 * a1 * a1 * np.conjugate(a0)
               
                elif corr =="4111-00":
                    q = a2 * a1 * a1 * a1 * np.conjugate(a0*a0)

                else:
                    print("Unknown correlator:  " + corr )
                    break

                Q[nt1+n1,nt2+n2]  += np.sum(q)
                # for a single set compute average and deviation

        ntot += len(q)

    Q  = Q/ntot
 
    Q.tofile( fbaseout + "_j" + str(j).zfill(3) + ".dat" )
 
    #print(Q.dtype, R.dtype,  Q.shape, R.shape)
    
    #P123 = np.fromfile(fout)
   
#==============================================================================

def Cmodel(corr, j, d):
#
# returns value of non-normalized irreducible cumulant predicted using data
# from Table in MMSupplement.
#
# corr - text such as "21-0"
# j    - mode number or array
# d    - damped mode, for inverse cascade d=1
# 
# corr - text: "21-0",  "322-0", "11-30", "222-40", "44333-0", "63333-0", "4111-00"


    import numpy as np

    jd = j - d
    
    if corr == "21-0":
        if d == 1:  #inverse
            E =  0.4358
        else: 
            E = -0.4241
        m  = 3
        mu = 1
        cj = (jd - 2)*(jd - 1)*jd
        
        
    elif corr == "322-0":
        if d == 1:  #inverse
            E =  0.1707
        else:
            E = -0.1634
        m  = 4
        mu = 2
        cj = (jd - 3)*(jd - 2)*(jd - 2)*jd


    elif corr == "431-0":
        if d == 1:  #inverse
            E =  0.0869
        else:
            E = -0.0879
        m  = 4
        mu = 1
        cj = (jd - 4)*(jd - 3)*(jd - 1)*jd


    elif corr == "11-30":
        if d == 1:  #inverse
            E = -0.4368
        else:
            E =  0.4381
        m  = 4
        mu = 2
        cj = (jd - 1)*(jd - 1)*(jd - 3)*jd

        
    elif corr == "222-40":
        if d == 1:  #inverse
            E =  0.0354 
        else:
            E = -0.0351
        m  = 5
        mu = 6
        cj = (jd - 2)**3 * (jd - 4)*jd

        
    elif corr == "44333-0":
        if d == 1:  #inverse
            E =  0.0264
        else:
            E = -0.0260
        m  = 6
        mu = 12
        cj = (jd - 4)**2 * (jd - 3)**3 *jd

        
    elif corr == "63333-0":
        if d == 1:  #inverse
            E = -0.0098
        else:
            E =  0.0094
        m  = 6
        mu = 24
        cj = (jd - 6) * (jd - 3)**4 *jd

        
    elif corr == "4111-00":
        if d == 1:  #inverse
            E = -0.0617
        else:
            E =  0.0629
        m  = 6
        mu = 12
        cj = (jd - 4) * (jd - 1)**3 * jd**2


        
    else:
        print("Unknown correlator: " + corr)
        return(0)

        
    Ctilde =  1.13**(m/3) * np.sqrt(mu) * (np.abs(cj))**(1/3)
    
    C = E * Ctilde / np.abs(jd)

    return(C)

    
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

