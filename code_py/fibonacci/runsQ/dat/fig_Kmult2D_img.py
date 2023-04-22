
#import importlib ; import fig_Kmult2D_img
#importlib.reload(fig_Kmult2D_img) ; fig_Kmult2D_img.plot()

#================================================================

def plot():
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cm

    #fnameout = "fig_Kmult2D_img7.pdf"
    fnameout = "fig_Kmult2D_img7.png"
    
    import matplotlib.pyplot as plt
       
    pi = np.pi

    phi = 1.61803398875

    equil = {0: -0.1604039416865561,    # -1/3 ln (phi)
             1: -0.3208078833731123,    # -2/3 ln (phi)
             -1: 0}

    eq = equil[0]
    extent0 = (-10 - eq, 10 - eq, -10 - eq, 10 - eq)

    lowlim = -7.5
    highlim =0
    
    #-- read data, alpha = 0 --

    fbase = "KM/q7_v0dt1_i20"
    
    imode, numbins, binsize, ncount =  ReadParamMI(fbase)
   
    P12  = np.fromfile(fbase + "_Pkmult.dat",  'int32')
    P12  = P12 / (ncount * binsize* binsize)
    P12  = P12.reshape(numbins, numbins)
    logPalpha0 =np.log10(P12).transpose()

    
    ind = (logPalpha0 < lowlim);  logPalpha0[ind] = lowlim;

    print(np.min(logPalpha0))

    #-- analytical probability --

    s = np.linspace(-10,10,401)
    s1,s2 = np.meshgrid(s,s)
    
    P = 8*np.exp(4*s1 + 2*s2) * (1 + np.exp(2*s1) + np.exp(2*s1 + 2*s2) )**(-3)

    logP = np.log10(P12)
   


    #-- plot 2D distributions of probabilities ---

    #levels = np.hstack((np.arange(-5,0,1), -0.55, 0, 1));     print(levels)
    #cmap="inferno_r"

    levels = np.hstack((lowlim,  np.arange(-5,0,1), -0.55, highlim));     print(levels)
    cmap="inferno"

    #levels = np.hstack((lowlim,  np.arange(-5,0,1), -0.55, 1));     print(levels)
    #cmap="RdYlBu_r"

    
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4,4))
    
    plt.xlim(-6.5,6.5)
    plt.ylim(-6.5,6.5)
 
    sc = plt.contourf(logPalpha0, levels, cmap=cmap, linestyles='solid', extent=extent0)
    #sc.cmap.set_under = "w"
    
    plt.contour(logPalpha0, levels[1:], colors="k", linestyles='solid', linewidths=1, extent=extent0)
    


    #norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
    #cset1 = axs.contourf(X, Y, Z, levels, norm=norm, cmap=cm.get_cmap(cmap, len(levels) - 1))
    

    #plt.imshow(logP, origin="lower", cmap="jet")


    ax.margins(x=0)

    plt.axis("off")
    plt.axis("off")
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    
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

