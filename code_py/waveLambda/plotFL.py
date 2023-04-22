# this is a comment

import kinwave as kw
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

#-------------------------------------------


#q = [-1.5, +1.5];  mr=24;  dr=1/8;   fbase="fld_pi2-1.5+1.5_r24o08"
#q = [+0.0, +1.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.0+1.5_r24o08"
#q = [+1.5, +1.5];  mr=24;  dr=1/8;   fbase="fld_pi2+1.5+1.5_r24o08"
#q = [+1.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi2+1.5+0.0_r24o08"
#q = [+1.5, -1.5];  mr=24;  dr=1/8;   fbase="fld_pi2+1.5-1.5_r24o08"
#q = [+0.0, -1.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.0-1.5_r24o08"
#q = [-1.5, -1.5];  mr=24;  dr=1/8;   fbase="fld_pi2-1.5-1.5_r24o08"
#q = [-1.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi2-1.5+0.0_r24o08"

#q = [-0.5, +0.5];  mr=24;  dr=1/8;   fbase="fld_pi2-0.5+0.5_r24o08"
#q = [+0.0, +0.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.0+0.5_r24o08"
#q = [+0.5, +0.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.5+0.5_r24o08"
#q = [+0.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi2+0.5+0.0_r24o08"
#q = [+0.5, -0.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.5-0.5_r24o08"
#q = [+0.0, -0.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.0-0.5_r24o08"
#q = [-0.5, -0.5];  mr=24;  dr=1/8;   fbase="fld_pi2-0.5-0.5_r24o08"
#q = [-0.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi2-0.5+0.0_r24o08"

#q = [-2.5, +2.5];  mr=24;  dr=1/8;   fbase="fld_pi2-2.5+2.5_r24o08"
#q = [+0.0, +2.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.0+2.5_r24o08"
#q = [+2.5, +2.5];  mr=24;  dr=1/8;   fbase="fld_pi2+2.5+2.5_r24o08"
#q = [+2.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi2+2.5+0.0_r24o08"
#q = [+2.5, -2.5];  mr=24;  dr=1/8;   fbase="fld_pi2+2.5-2.5_r24o08"
#q = [+0.0, -2.5];  mr=24;  dr=1/8;   fbase="fld_pi2+0.0-2.5_r24o08"
#q = [-2.5, -2.5];  mr=24;  dr=1/8;   fbase="fld_pi2-2.5-2.5_r24o08"
#q = [-2.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi2-2.5+0.0_r24o08"


#q = [-1.5, +1.5];  mr=24;  dr=1/8;   fbase="fld_pi4-1.5+1.5_r24o08"
#q = [+0.0, +1.5];  mr=24;  dr=1/8;   fbase="fld_pi4+0.0+1.5_r24o08"
#q = [+1.5, +1.5];  mr=24;  dr=1/8;   fbase="fld_pi4+1.5+1.5_r24o08"
#q = [+1.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi4+1.5+0.0_r24o08"
#q = [+1.5, -1.5];  mr=24;  dr=1/8;   fbase="fld_pi4+1.5-1.5_r24o08"
#q = [+0.0, -1.5];  mr=24;  dr=1/8;   fbase="fld_pi4+0.0-1.5_r24o08"
#q = [-1.5, -1.5];  mr=24;  dr=1/8;   fbase="fld_pi4-1.5-1.5_r24o08"
#q = [-1.5, +0.0];  mr=24;  dr=1/8;   fbase="fld_pi4-1.5+0.0_r24o08"

fname = "2023-04-12_fld/" + fbase + ".npy"


#kw.N2Dfld_pi2(q, mr, dr, fname)
   
#--------------------------------------------

def plot_Nfld(fbase):

    fname = "2023-04-12_fld/" + fbase + ".npy"
    outfile = "2023-04-12_fld/" + fbase + ".pdf"
    
    
    with open(fname, 'rb') as f:
        p   = np.load(f)
        p1  = np.load(f)
        q   = np.load(f)
        (mr, dr, S) = np.load(f) 
        NF = np.load(f)
        NL = np.load(f)
        m = int ( (NF.shape[0]-1)/2 )
    

        
    print(m,p,p1,q, np.shape(NF))
    print(mr, dr, S)
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,4)) 

    Lrange = 100
    Frange = 50
    
    #Lrange = np.max( [np.abs(np.max(NL)), np.abs(np.min(NL))] )
    #print("Lrange = ", Lrange)
    #Frange = np.max( [np.abs(np.max(NF)), np.abs(np.min(NF))] )/10
    #print("Frange = ", Frange)

    extent=((-m-0.5)*dr, (m+0.5)*dr, (-m-0.5)*dr, (m+0.5)*dr)

    
    ax[0].imshow(np.transpose(NL), cmap="jet", origin = 'lower', vmin=-Lrange, vmax=Lrange, extent=extent)
    ax[1].imshow(np.transpose(NF), cmap="jet", origin = 'lower', vmin=-Frange, vmax=Frange, extent=extent)

    ax[0].set_title( "q = [{}, {}], \ \  p = [1,0], \ \  p\' = [1,1]/sqrt(2)".format(q[0], q[1]))
    ax[1].set_title( "$m_r = {}, \quad \Delta r = 1/{}$".format(int(mr), int(1/dr)) )

    
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.03)


    plt.savefig(outfile) #pad_inches=0

    plt.close()


    return() 

plot_Nfld(fbase)

print("Done.")




#===========================================



#=========================================================
