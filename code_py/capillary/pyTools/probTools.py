#import importlib ; import gfProbFlat
#importlib.reload(gfProbFlat) ; gfProbFlat.Test(iz=1)

def Test():
   
 
       
#=================================================================================================


#=================================================================================================

def probTrim(cI, dIo):
    
    import numpy as np

    q = np.nonzero(cI>0)[0][-1] + 1
    cIr = cI[:q]
    
    Ir = (np.arange(len(cIr)) + 0.5) * dIo
    dIr = (np.ones(len(cIr))) * dIo

    
    pIr  = cIr  / np.sum(cIr)  / dIo
   
    return Ir, pIr, dIr


#=================================================================================================

def propRebin(cI, dIo, ranges, debug=False):
    
    import numpy as np
    
    cIr = list({})
    dIr = list({})
    Ir  = list({})
    
    i1 = 0
    
    for r in ranges:
        
        i2 = r[1]/dIo
        n = max(0, int(r[0]/dIo))
        
        while (i1 < i2):
            c  = np.sum(cI[i1:i1+n])
            center = dIo * np.sqrt(i1*(i1+n))
            width = dIo * n
            cIr.append(c)            
            Ir.append(center)
            dIr.append(width)

            if debug:
                print("first bin = {}, bins = {}, interval = [{}, {}], count {}".format(
                    i1, n, i1*dIo, (i1+n)*dIo, c))
            i1 += n
                  
    cIr = np.asarray(cIr)
    Ir  = np.asarray(Ir)
    dIr = np.asarray(dIr)

    q = np.nonzero(cIr>0)[0][-1] + 1
    cIr = cIr[:q]
    Ir  = Ir[:q]
    dIr = dIr[:q]
    pIr = cIr / np.sum(cIr) / dIr


    if debug:

        I = (np.arange(len(cI)) + 0.5) * dIo
        pI  = cI  / np.sum(cI)  / dIo

        Iavg  = np.sum(I  * pI * dIo) 
        Iravg = np.sum(Ir * pIr * dIr) 
    
        print("Totals: ", np.sum(cI), np.sum(cIr)),  
 
        print("Averages: ", Iavg,  Iravg) 

        #-- show plot --
        
        import matplotlib.pyplot as plt

        fig, ax1 = plt.subplots(1, 1, figsize=(6,4))

        ax1.grid()
 
        ax1.set_yscale('log')

        q = np.nonzero(cI>0)[0][-1] + 1

        ax1.set_xlim(0,q*dIo)
        ax1.set(xlabel="$I$", ylabel="PDF")
    
        ax1.plot(I, cI/np.sum(cI)/dIo, 'ob',  markersize=2, linewidth=1,  mfc='none')
        ax1.plot(Ir, pIr, 'or',  markersize=4, linewidth=1,  mfc='none')
        
        plt.tight_layout()
        
        plt.show()
    
    return Ir, pIr, dIr


#=================================================================================================


def probBinEqual(f,  numbins=2000,  bincoef=0.01):

    import numpy as np

    favg = np.average(f)

    binsize = bincoef*favg    
    binedges = np.arange(numbins)*binsize
    
    (p, x) = np.histogram(f, binedges, density=True)
    
    x  = 0.5 * (binedges[:-1] + binedges[1:])
    dx = np.vstack((x-binedges[:-1],  binedges[1:]-x))

    return(p, x, dx)

#=================================================================================================

def probBinStretch(f,  numbins=50,  coef0=1e-2, coeflast=1e2):

    import numpy as np

    favg = np.average(f)

    binsize0 = coef0*favg
    lastedge = coeflast*favg
    c = np.log(1+ lastedge/binsize0)/numbins
    binedges = (np.exp(c * np.arange(numbins+1)) - 1)*binsize0
    
    (p, x) = np.histogram(f, binedges, density=True)
    
    x = np.sqrt( binedges[:-1] * binedges[1:] )
    dx = np.vstack((x-binedges[:-1],  binedges[1:]-x))

    ind = (p>0)
    x = x[ind]
    p = p[ind]
    
    return (p, x, dx)

#=================================================================================================

