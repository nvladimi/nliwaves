
#import importlib ; import twomode as tm;
#importlib.reload(tm) ; tm.Test()


def Test():
    
    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt


    #   plotICDC(datdir = "POSTMI/1", startswith="ig", bins=(100,100,32))
    #   plotICDC(datdir = "POSTMI/3", startswith="ig", bins=(200,40,32))
    #   plotICDC(datdir = "POSTMI/3", startswith="dg", bins=(20,30,32))

    miICDC()
    
#   Probability3D("dg001_p1e2_dt1", istart=1, iend=-1, binsize=1, numbins=(100, 50, 32), fbaseout='test')
#   MI('test')
     
#   ProbabilityN(dat, numbins=10, binsize=nu1*10, fbaseout=fbase,  showplot=False)
#   ProbabilityPhi(dat, numbins=32, fbaseout=fbase,  showplot=True)
#   Bursts(dat, fbaseout="", showplot=True)

    return



#================================================================
    
#================================================================


def miICDC(datdir = "POSTMI/5", startswith="y", outpkl = "./POSTMI/mi5_ic.pkl"):
    
    import numpy as np
    import pickle
    import os

    dat = list()

    #-- read all files in directory, compute MI, store in"dat" --

    for fname in os.listdir(datdir):
        if fname.endswith("B_mi.pkl") & fname.startswith(startswith):
            with open(datdir +"/"+ fname, 'rb') as handle:
                dat0 = pickle.load(handle)
                dat.append(dat0)


    #-- sort by "chi" -- 
                
    keys = dat0.keys()
    print(keys)

    def myFunc(e):
        return e['chi']

    dat.sort(key=myFunc)
    
    for d in dat:
        print("{:10.4f} {:6f} {:6f} {:6f} {:d}".format(d['chi'], d['I12'], d['n1'], d['n2'], d['ntot']))


    #-- create new dictionary of arrays --

    newdict   = dict.fromkeys(keys, 0)

    for key in keys: 
        newvalues = list()
        for d in dat:
            newvalues.append(d[key])
        newdict[key] = np.array(newvalues)


    print(newdict["chi"])
    print(newdict["I12"])
    
            
    #-- save dictionary of arrays as pickle --

    with open(outpkl, 'wb') as handle:
        pickle.dump(newdict, handle, protocol=pickle.HIGHEST_PROTOCOL)


#================================================================

def plotICDC(datdir = "POSTMI/1", startswith="ig", bins=(100,100,32)):
    
    import numpy as np
    import pickle
    import os
    import matplotlib.pyplot as plt

    dat = list()

    #-- canvas --


    fig, ax = plt.subplots(ncols=2, nrows=1, tight_layout=True, figsize=(12,6))
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[0].set(xlabel = "rho1", ylabel = "bincount")  
    ax[1].set(xlabel = "rho2", ylabel = "bincount")  
    
    n1 = np.arange(bins[0])
    n2 = np.arange(bins[1])
    

    #-- read all files in directory, compute MI, store in"dat" --

    for fname in os.listdir(datdir):
        if fname.endswith("_P123.dat") & fname.startswith(startswith):
            P123 = np.fromfile(datdir + "/" + fname, dtype="uint32")
            P123 = np.reshape(P123, bins) 
            P1  = P123.sum((1,2))
            P2  = P123.sum((0,2))
            ax[0].plot(n1, P1, label= fname[:-8].replace("_", " ") )
            ax[1].plot(n2, P2, label= fname[:-8].replace("_", " ") )

    import operator
    handles, labels = ax[0].get_legend_handles_labels()
    hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)
    ax[0].legend(handles2, labels2)


    import operator
    handles, labels = ax[1].get_legend_handles_labels()
    hl = sorted(zip(handles, labels), key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)
    ax[1].legend(handles2, labels2)

    plt.show()


#================================================================


def Probability3D(fbase, istart=1, iend=-1, binsize=1, numbins=(32, 32, 32), fbaseout=''):
    
    import numpy as np
    import pickle

     
    pi = np.pi

    dat = ReadData(fbase, istart, iend)

    
    #-- setup bins --
    
    n1avg, n2avg = averageN(dat)
   
    edg1 = np.arange(numbins[0]+1)*binsize*n1avg;  edg1[-1] = 1.e16
    edg2 = np.arange(numbins[1]+1)*binsize*n2avg;  edg2[-1] = 1.e16 
    edgPhi =  (np.arange(-numbins[2]/2-1, numbins[2]/2+1) + 0.5 )/numbins[2]*2*pi
    
    edge123  = (edg1, edg2, edgPhi)

    P123 = np.zeros((numbins[0], numbins[1], numbins[2]+1), 'int32')

    
    #-- collect propability --
    
    n1 = dat[:,1]*dat[:,1] + dat[:,2]*dat[:,2]
    n2 = dat[:,3]*dat[:,3] + dat[:,4]*dat[:,4]
   
    phi1 = np.angle(dat[:,1] + 1j*dat[:,2])
    phi2 = np.angle(dat[:,3] + 1j*dat[:,4])
                 
    dphi = 2*phi1 - phi2 
        
    ind = (dphi >  pi);  dphi[ind] += -2*pi;
    ind = (dphi < -pi);  dphi[ind] +=  2*pi;
 
    u = np.array([n1, n2, dphi]).T
    
    (h, edges) = np.histogramdd(u, edge123)
    P123 += h.astype('uint32')
                
    P123[:,:,0] +=  P123[:,:,-1]
    P123 = P123[:,:,:-1]

    ncount = np.sum(P123)
    
    #-- save data --

    [p1, p2,  g1, g2, dt] = ReadParam(fbase)

    T = (p1 + 2*p2)/2/(g1+g2)
    chi = (g1 + g2)**2 /2/T

    try:
        dT = p1/g1 - 2*p2/g2
        dTT = dT/T
    except:
        if g1 == 0:
            dTT =  999
        else:
            dTT = -999

    dictout = {"chi":      chi,
              "dTT":      dTT,
              "p1":       p1,
              "p2":       p2,
              "g1":       g1,
              "g2":       g2,
              "n1":       n1avg,
              "n2":       n2avg,
              "ntot":     ncount,
              "binsize":  binsize,
              "numbins":  numbins}
    

    if (len(fbaseout) > 0):

        with open(fbaseout + "_P123param.pkl", 'wb') as handle:
            pickle.dump(dictout, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
        P123.tofile(fbaseout + "_P123.dat")
 
        
        fnameout = fbaseout + "_P123param.txt"
       
        f = open(fnameout, 'w')
        
        f.write("\n")

        f.write("fbase:       " + fbase + "\n")
        f.write("istart:      " + str(istart) + "\n")
        f.write("iend:        " + str(iend)   + "\n")
        
        f.write("\n")
        
        f.write("ncount:      " + str(ncount) + "\n")
        f.write("n1avg:       " + str(n1avg)  + "\n")
        f.write("n2avg:       " + str(n2avg)  + "\n")

        f.write("\n")

        f.write("binsize     " + str(binsize)    + "\n")
        f.write("numbinsN1   " + str(numbins[0]) + "\n")
        f.write("numbinsN2   " + str(numbins[1]) + "\n")
        f.write("numbinsNphi " + str(numbins[2]) + "\n")

        f.write("\n")

        f.close()

#================================================================

def MI(fbase):

    import numpy as np
    import pickle
    
    pi = np.pi

    
    with open(fbase + "_P123param.pkl", 'rb') as handle:
        param = pickle.load(handle)

    P123 = np.fromfile(fbase + "_P123.dat", dtype="uint32")
    
    numbins = param["numbins"]
    nN1=numbins[0]
    nN2=numbins[1]
    nPhi=numbins[2]
    
    dn1  = param["n1"] * param["binsize"]
    dn2  = param["n2"] * param["binsize"]
    dphi = 2*pi/nPhi

    P123 = np.reshape(P123, (nN1, nN2, nPhi))
    
    #-- normalize and integrate probabilities --
    
    ncount = P123.sum()

    if ncount != param["ntot"]:
         print('ERROR: Mispatch in size of data, ntot = {}, nread={}'.format(param['ntot'], ncount))
    
    P123 = P123 / (ncount *dn1*dn2*dphi)
        
    P1   = P123.sum((1,2))    *dn2*dphi
    P2   = P123.sum((0,2))*dn1    *dphi
    Pphi = P123.sum((0,1))*dn1*dn2

    
    #--  averages for staircase distribution --

    n = (np.arange(0,nN1)+0.5)*dn1
    navg1s = np.sum(P1*n)/np.sum(P1)

    n = (np.arange(0,nN2)+0.5)*dn2
    navg2s = np.sum(P2*n)/np.sum(P2)
    
 
    #-- Entropy of 3D distribution (n1,n2,theta) --
    
    ind = (P123>0);  S12 = -np.sum(P123[ind]* np.log2(P123[ind])) *dn1*dn2*dphi

 
    #-- Entropy of 1D distributions --

    ind = (P1>0);    S1  = -np.sum(P1[ind] * np.log2(P1[ind])) * dn1
    ind = (P2>0);    S2  = -np.sum(P2[ind] * np.log2(P2[ind])) * dn2
    ind = (Pphi>0);  Sphi  = -np.sum(Pphi[ind] * np.log2(Pphi[ind])) * dphi

      
    #-- Equilibrium, assuming exponential fit with <n_i> and equidistribution for theta --
    
    
    S1s = 1/np.log(2) * (1 + np.log(navg1s));
    S2s = 1/np.log(2) * (1 + np.log(navg2s));
    Sph0 = np.log2(2*pi);
    

    #-- Mutual information -- 

    
    I12  = S1 + S2 + Sph0 - S12 
 

    #-- Output to dictionary --

    dictout = {"chi":  param["chi"],
              "dTT":  param["dTT"],
              "S12":  S12,
              "S1":   S1,
              "S2":   S2,
              "Sphi": Sphi,
              "I12":  I12,
              "n1s":  navg1s,
              "n2s":  navg2s,
              "n1":   param["n1"],
              "n2":   param["n2"],
              "ntot": ncount
              }
    
    with open(fbase + "_mi.pkl", 'wb') as handle:
        pickle.dump(dictout, handle, protocol=pickle.HIGHEST_PROTOCOL)
 

#================================================================
    
def averageN(dat):

    import numpy as np
     
    n1 = dat[:,1]*dat[:,1] + dat[:,2]*dat[:,2]
    n2 = dat[:,3]*dat[:,3] + dat[:,4]*dat[:,4]

    n1avg = np.average(n1)
    n2avg = np.average(n2)
 
    return n1avg, n2avg


#================================================================
    
def ProbabilityN(dat, numbins, binsize, fbaseout = "",  showplot = False):

    import numpy as np
     
    n1 = dat[:,1]*dat[:,1] + dat[:,2]*dat[:,2]
    n2 = dat[:,3]*dat[:,3] + dat[:,4]*dat[:,4]

    n1avg = np.average(n1)
    n2avg = np.average(n2)
    
    binedges = np.arange(numbins)*binsize

    (n1prob, x) = np.histogram(n1, binedges, density=True)
    (n2prob, x) = np.histogram(n2, binedges, density=True)

    x = 0.5 * (binedges[:-1] + binedges[1:])

    if (len(fbaseout) > 0):
        dataout = np.vstack((x, n1prob, n2prob)).transpose()
        h = "Run \"" + fbaseout + "\":  " + str(len(n1)) + " elements"
        h = h + " with n1avg = " + str(n1avg) + ",  n2avg = " + str(n2avg) + "\n"
        h = h + "1.n1,n2   2.prob(n1)   3.prob(n2)\n"
        fnameout = fbaseout + "_probN.txt"
        np.savetxt(fnameout, dataout, header = h)
    
    if showplot:
        
        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots();
        plt.yscale('log')
        ax.plot(x,n1prob, 'ro');
        ax.plot(x,n2prob, 'bo');

        ax.set(xlabel='$n_1$,  $n_2$', ylabel='probability', title = fbaseout)
        plt.show()

    return


#================================================================

def ProbabilityPhi(dat, numbins, fbaseout = "",  showplot = False):

    import numpy as np

    pi = np.pi
    
    a1=np.angle(dat[:,1] + 1j*dat[:,2]);
    a2=np.angle(dat[:,3] + 1j*dat[:,4]);
    phi = 2*a1 - a2;

    del a1, a2
    
    ind = (phi>pi);  phi[ind] = phi[ind] - 2*pi;
    ind = (phi<-pi); phi[ind] = phi[ind] + 2*pi;
    
    binedges =  (np.arange(-numbins/2-1, numbins/2+1) + 0.5 )/numbins*2*pi 
    
    (prob, x) = np.histogram(phi, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])

    prob[0] = prob[0] + prob[-1]

    prob = prob[:-1]
    x = x[:-1]

    x = x/pi
    prob = prob*pi
    
    if (len(fbaseout) > 0):
        dataout = np.vstack((x, prob)).transpose()
        h = "Run \"" + fbaseout + "\":  " + str(len(phi)) + " elements\n"
        h = h + "1.phi   2.prob(phi)\n"
        fnameout = fbaseout + "_probPhi.txt"
        np.savetxt(fnameout, dataout, header = h)

    if showplot:

        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots();
        ax.plot(x,prob, '-ro');

        ax.set(xlabel='$\phi / \pi$', ylabel='probability', title = fbaseout)

        plt.grid(True)
        plt.show()

        
    return


#================================================================
    
def Bursts(dat, spacing = 0.1, fbaseout = "",  showplot = False):

    import numpy as np

    #-- extract time dependence of occupation number n1(t) from data --
    
    t  = dat[:,0]
    n1 = dat[:,1]*dat[:,1] + dat[:,2]*dat[:,2]
    
    n1avg = np.average(n1)

    #-- find all local maxima in n1 array --
    
    f0 = n1[1:-1]
    fm = n1[0:-2]
    fp = n1[2:]
    
    q = ( (f0>fm) & (f0>fp) & (f0>2*n1avg) )
    q = np.hstack((np.full((1,), False), q))
    
    indmax = np.nonzero(q)

    del f0, fm, fp, q

    n1max = n1[indmax];
    tmax  = t[indmax];
    
    ismax = np.full(n1max.shape, True)  # array indicating which maxima are isolated

    #-- find false maxima that are too close to others, set corrsponding ismax to false --

    n1maxsort = np.flip(np.argsort(n1max)) 
     
    for i in n1maxsort:
        if (ismax[i]):
             t0 = tmax[i]
             n0 = n1max[i]
             dtmax = abs(tmax - t0)
             dtThreshold = spacing/n1max[i]
             ind =  (( dtmax < dtThreshold ) & (n1max < n0))
             ind[i] = False
             ismax[ind] = False
             
    #-- create array of selected maxima --
    
    n1Burst = n1max[ismax]
    tBurst  = tmax[ismax]

    if (len(fbaseout) > 0):
        dataout = np.vstack((tBurst, n1Burst)).transpose()
        h = "Run \"" + fbaseout + "\":  " + str(len(tBurst)) + " bursts"
        h = h + " with spacing > " + str(spacing) + "\n"
        h = h + "1.t_burst   2.n1_burst\n"
        fnameout = fbaseout + "_bursts.txt"
        np.savetxt(fnameout, dataout, header = h)

    #-- plot maxima superimposed to time dependence --

    if showplot:

        import matplotlib
        matplotlib.rcParams['text.usetex'] = True
        import matplotlib.pyplot as plt
        
        print("Number of bursts:  " + str(len(tBurst)) )   

        #plt.ion()    
        fig, ax = plt.subplots();
        ax.plot(t,n1, 'r');
        ax.plot(tBurst,n1Burst, 'bo');
        ax.grid();
        ax.set(xlabel='t', ylabel='$N_1$', title = fbaseout)
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
        
        if ( s2[0].count('isave') > 0):
            isave = int(s2[1])
        if ( s2[0].count('nsave') > 0):
            nsave = int(s2[1])
        if ( s2[0].count('dt') > 0):
            dt = float(s2[1])
        if ( s2[0].count('Gamma') > 0):
            s4 = s2[1].strip().split(" ")
            gamma = [float(s0) for s0 in s4]
        if ( s2[0].count('Rflux') > 0):
            s4 = s2[1].strip().split(" ")
            rflux = [float(s0) for s0 in s4]       

    p1 =  rflux[0]
    p2 =  rflux[2]
    g1 = -gamma[0]
    g2 = -gamma[2]
            
    return([p1, p2,  g1, g2, dt])

#================================================================


def ReadData(fbase, istart = 1,  iend = -1):

    import os
    import numpy as np
 
    i = istart;
    dat = np.empty( shape=(0, 5) )
    
    while True:

        fname = fbase + '/'+ fbase + '.' + str(i).zfill(4) + '.a1a2'
        
        if (os.path.isfile(fname) &  ((iend < 0) | (i < iend)) ) :
            dat0 = np.fromfile(fname, 'float64')
            dat0 = dat0.reshape(5, round(len(dat0)/5)).transpose()
            dat = np.vstack((dat, dat0))
            i = i+1
        else:
            print("Last file read:  " + fname + ";    datapoints:  " +  str(len(dat)) )
            break
        
    return(dat)    


#================================================================

