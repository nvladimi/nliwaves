
#import importlib ; import twomodeMI
#importlib.reload(twomodeMI) ; twomodeMI.allMI()
   
        
#================================================================


def allMI(datdir = "./POST/post05", datpkl = "./POST/post05_mi.pkl"):
    
    import numpy as np
    import pickle
    import os

    dat = list()

    #-- read all files in directory, compute MI, store in"dat" --

    for file in os.listdir(datdir):
        if file.endswith(".P123"):
            fbase = file[:-5]
            dat0 = MI(datdir + "/" + fbase)
            dat.append(dat0)

    keys = dat0.keys()

    #-- find and sort all "chi" --
    
    allchi = set()    
    for d in dat:
        allchi.add(d['chi'])
    allchi = list(allchi)
    allchi.sort()

    
    #-- for each "chi" sort by "dTT" --
    
    def myFunc(e):
        return e['dTT']

    alldict = list()
    
    for chi in allchi:

        thisdict = dict.fromkeys(keys, 0)
        
        dat0 = list()
        for d in dat:
            if d['chi'] == chi:
                dat0.append(d)
        dat0.sort(key=myFunc)

        for key in keys:
            values = list()
            for d in dat0:
                values.append(d[key])
            thisdict[key] = np.array(values)

        thisdict['chi'] = round(chi, 6)

        alldict.append(thisdict)
        
    #-- save list of dictionaries as pickle --

    with open(datpkl, 'wb') as handle:
        pickle.dump(alldict, handle, protocol=pickle.HIGHEST_PROTOCOL)


#================================================================

def MI(fbase):

    import numpy as np
    
    pi = np.pi
    
    param = ReadParam(fbase + ".param")
    P123  = ReadData(fbase + ".P123", param)  #P123[Nphi, nN2, nN1]

    nN1=param["nN1"]
    nN2=param["nN2"]
    nPhi=param["nPhi"]
    
    dn1  = param["max_N1"]/nN1
    dn2  = param["max_N2"]/nN2
    dphi = 2*pi/nPhi

    
    #-- normalize and integrate probabilities --
    
    ncount = P123.sum()
    
    P123 = P123 / (ncount *dn1*dn2*dphi)
        
    P1   = P123.sum((0,1))    *dn2*dphi
    P2   = P123.sum((0,2))*dn1    *dphi
    Pphi = P123.sum((1,2))*dn1*dn2

    
    #--  averages for staircase distribution --

    n = (np.arange(0,nN1)+0.5)*dn1
    navg1s = sum(P1*n)/sum(P1)

    n = (np.arange(0,nN2)+0.5)*dn2
    navg2s = sum(P2*n)/sum(P2)
    
 
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

    g1 = param["g1"]
    g2 = param["g2"]
    p1 = param["p1"]
    p2 = param["p2"]
    
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


    datout = {"chi":  chi,
              "dTT":  dTT,
              "S12":  S12,
              "S1":   S1,
              "S2":   S2,
              "Sphi": Sphi,
              "I12":  I12,
              "n1s":  navg1s,
              "n2s":  navg2s,
              "ntot": ncount
              }

    return(datout)
  
#================================================================
    
def ReadParam(fname):
  
    param = {
        "ntot":   0,
        "nPhi":   0,
        "nN1":    0,
        "nN2":    0,
        "max_N1": 0.0,
        "max_N2": 0.0,
        "g1":     0.0,
        "g2":     0.0,
        "p1":     0.0,
        "p2":     0.0
    }


    f = open(fname, 'r')
    s = f.read()
    f.close()
    
    ss = s.split("\n\n")

    for s in ss:
        
        s1 = s.replace("\n#", ";")
        s2 = s1.split("\n")        

        for i in param:

            if ( s2[0].count('name: '+i) > 0):
                param[i] = type(param[i])(s2[1])
                      
    return(param)

  
#================================================================

def ReadData(fname, param):

    import os
    import numpy as np

    dat = np.fromfile(fname, 'uint32')
    dat = dat.reshape(param['nPhi'], param['nN2'], param['nN1'])

    if param['ntot'] !=  np.sum(dat):
        print('ERROR: Mispatch in size of data, ntot = {}, nread={}'.format(param['ntot'], np.sum(dat)))
 
    return(dat)    



#================================================================
