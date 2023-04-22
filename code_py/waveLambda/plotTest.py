# this is a comment

import kinwave as kw
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

#-------------------------------------------



q = [+1.5, +1.5];  mr=16;  dr=1/4;   fbase="fldFtest"


fname = fbase + ".npy"


kw.testNfld(q, mr, dr, fname)
   
#--------------------------------------------

def plot_N1fld(fbase):

    fname = fbase + ".npy"
    #outfile = "2023-04-13_fld1/" + fbase + ".pdf"
    outfile = "plotTest.pdf" 
    
    
    with open(fname, 'rb') as f:
        p   = np.load(f)
        p1  = np.load(f)
        q   = np.load(f)
        (mr, dr, S) = np.load(f) 
        NF = np.load(f)
        m = int ( (NF.shape[0]-1)/2 )
    

        
    print(m,p,p1,q, np.shape(NF))
    print(mr, dr, S)
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,4)) 

    Lrange = 100
    Frange = 50
    
    #Lrange = np.max( [np.abs(np.max(NL)), np.abs(np.min(NL))] )
    #print("Lrange = ", Lrange)
    Frange = np.max( [np.abs(np.max(NF)), np.abs(np.min(NF))] )
    print("Frange = ", Frange)

    extent=((-m-0.5)*dr, (m+0.5)*dr, (-m-0.5)*dr, (m+0.5)*dr)

    
    #ax[0].imshow(np.transpose(NL), cmap="jet", origin = 'lower', vmin=-Lrange, vmax=Lrange, extent=extent)
    ax[1].imshow(np.transpose(NF), cmap="jet", origin = 'lower', vmin=-Frange, vmax=Frange, extent=extent)

    ax[0].set_title( "q = [{}, {}], \ \  p = [1,0], \ \  p\' = [0,1]".format(q[0], q[1]))
    ax[1].set_title( "$m_r = {}, \quad \Delta r = 1/{}$".format(int(mr), int(1/dr)) )

    
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.03)


    plt.savefig(outfile) #pad_inches=0

    plt.close()


    return() 

#plot_N1fld(fbase)

print("Done.")




#===========================================



#=========================================================
