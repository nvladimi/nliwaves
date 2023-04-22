
#import importlib ; import fig_Kmult2D
#importlib.reload(fig_Kmult2D) ; fig_Kmult2D.plot()

#================================================================

def plot():
    
    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    fnameout = "fig_Kmult2D.pdf"
    
    import matplotlib.pyplot as plt
       
    pi = np.pi

    phi = 1.61803398875

    equil = {0: -0.1604039416865561,    # -1/3 ln (phi)
             1: -0.3208078833731123,    # -2/3 ln (phi)
             -1: 0}


    #-- read data, alpha = 0 --

    fbase = "KM/q7_v0dt1_i20"
    
    imode, numbins, binsize, ncount =  ReadParamMI(fbase)
   
    P12  = np.fromfile(fbase + "_Pkmult.dat",  'int32')
    P12  = P12 / (ncount * binsize* binsize)
    P12  = P12.reshape(numbins, numbins)
    logPalpha0 =np.log10(P12).transpose()
 
    

    
    #-- read data, alpha = 1 --
    
    fbase = "KM/q7b_g3500dt5_i40"

    imode, numbins, binsize, ncount =  ReadParamMI(fbase)
   
    P12  = np.fromfile(fbase + "_Pkmult.dat",  'int32')
    P12  = P12 / (ncount * binsize* binsize)
    P12  = P12.reshape(numbins, numbins)
    logPalpha1 =np.log10(P12).transpose()


    #-- analytical probability --

    s = np.linspace(-10,10,401)
    s1,s2 = np.meshgrid(s,s)

    
    P = 8*np.exp(4*s1 + 2*s2) * (1 + np.exp(2*s1) + np.exp(2*s1 + 2*s2) )**(-3)

    logP = np.log10(P)

    #plt.imshow(logP, origin="lower", cmap="jet")

    #plt.show()

    #return()
        
    #-- plot 2D distributions of probabilities ---

    levels = np.hstack((np.arange(-5,0,1), -0.55))
    level0 = (-0.55,)
    
 
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4,4),constrained_layout=True)
    
    plt.xlim(-7,7)
    plt.ylim(-7,7)
    plt.xlabel("$\sigma_i$")
    plt.ylabel("$\sigma_{i+1}$")
    plt.axvline(c='grey', lw=0.5)
    plt.axhline(c='grey', lw=0.5)

    #plt.grid(True)
        
    eq = equil[0]
    extent0 = (-10 - eq, 10 - eq, -10 - eq, 10 - eq)
    
    eq = equil[1]
    extent1 = (-10 - eq, 10 - eq, -10 - eq, 10 - eq)

    
    plt.contour(logPalpha0, levels,
                colors="red", linestyles='solid', linewidths=1,
                extent=extent0)
    plt.contour(logPalpha1, levels,
                colors="#0080FF", linestyles='dashed', linewidths=1,
                extent=extent1) 
    plt.contour(logP, levels,
                colors="black", linestyles='solid', linewidths=0.5,
                extent=(-10,10,-10,10))

    #ax = fig.add_subplot(333)


    # inset axes....
    axins = ax.inset_axes([0.68, 0.68, 0.305, 0.305])
    x1, x2, y1, y2 = -0.5, 0.5, -0.5, 0.5
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    #axins.set_xticklabels('')
    #axins.set_yticklabels('')
    axins.set_xticks(())
    axins.set_yticks(())
    ax.indicate_inset_zoom(axins)


    axins.contour(logPalpha0, level0,
                colors="red", linestyles='solid', linewidths=1,
                extent=extent0)
    axins.contour(logPalpha1, level0,
                colors="#0080FF", linestyles='dashed', linewidths=1,
                extent=extent1) 
    axins.contour(logP, level0,
                colors="black", linestyles='solid', linewidths=0.5,
                extent=(-10,10,-10,10))

    axins.axvline(c='grey', lw=0.5)
    axins.axhline(c='grey', lw=0.5)


    
     
    plt.savefig(fnameout, pad_inches=0)
        
    #plt.show()
        
    plt.close()
        
    return
    


  
#================================================================
    
def ReadParamMI(fbase):
  
    fname = fbase + '_Kmult.param'

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

