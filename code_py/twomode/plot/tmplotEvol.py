
#import importlib ; import tmplotEvol
#importlib.reload(tmplotEvol) ; tmplotEvol.Evol()



#================================================================

def Evol():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D

    outfile = 'PLOT/avgNevol.pdf'

    prefix   =  "POST/"
    suffix1  = "_dt1_tavg.txt";
    suffix2  = "_dt2_tavg.txt";

    #-------------------------------------------
    
    fbases1 = ("dg001_p1e4", "dg001_p1e3", "dg001_p1e2")
    pp    = np.asarray([1e-4, 1e-3, 1e-2])
    gg    = np.asarray([0.01, 0.01, 0.01])
    chis   = gg*gg*gg / pp
    nus1   = pp/(4*gg)

    fbases2 = ("dg100_p1e2", "dg100_p4e3", "dg100_p2e3")
    pp    = np.asarray([1e-2, 4e-3, 2e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    chis   = gg*gg*gg / pp
    nus2   = pp/(4*gg)


    fbases3 = ("ig001_p5e5", "ig001_p5e4", "ig001_p5e3")
    pp    = np.asarray([5e-5, 5e-4, 5e-3])
    gg    = np.asarray([0.01, 0.01, 0.01])
    chis   = gg*gg*gg / (2*pp)
    nus3   = pp/gg

    fbases4 = ("ig100_p5e3", "ig100_p2e3", "ig100_p1e3")
    pp    = np.asarray([5e-3, 2e-3, 1e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    chis   = gg*gg*gg / (2*pp)
    nus4   = pp/gg


    
    
    #-------------------------------------------
    
    fig, f1_axes = plt.subplots(ncols=2, nrows=4, constrained_layout=True, figsize=(8,10))

    ax1a = f1_axes[0,0]
    ax1b = f1_axes[0,1]
    ax2a = f1_axes[1,0]
    ax2b = f1_axes[1,1]
    ax3a = f1_axes[2,0]
    ax3b = f1_axes[2,1]
    ax4a = f1_axes[3,0]
    ax4b = f1_axes[3,1]


    ax1a.set(ylabel='$\\rho_1^2 / \\nu_2$', title = "direct cascade, small $\chi$")
    ax1b.set(ylabel='$\\rho_2^2 / \\nu_2$', title = "direct cascade, small $\chi$")
 
    ax2a.set(ylabel='$(\\rho_1^2 / \\nu_2 ) \chi^{-0.45} $', title = "direct cascade, large $\chi$")
    ax2b.set(ylabel='$\\rho_2^2 / \\nu_2$', title = "direct cascade, large $\chi$")
 
    ax3a.set(ylabel='$\\rho_1^2 / \\nu_1$', title = "inverse cascade, small $\chi$")
    ax3b.set(ylabel='$\\rho_2^2 / \\nu_1$', title = "inverse cascade, small $\chi$")
 
    ax4a.set(ylabel='$\\rho_1^2 / \\nu_1$', title = "inverse cascade, large $\chi$")
    ax4b.set(ylabel='$\\rho_2^2 / \\nu_1$', title = "inverse cascade, large $\chi$")
 

    #ax1.set_yscale('log')
    #ax1.set_xlim(0,35)
    #ax1.set_ylim(1e-6,1)

    LW = 0.5
    
    #-------------------------------------------
    
    suffix = suffix1    
    fnamein = prefix + fbases1[0] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus1[0]
   
    ax1a.plot(dat[:,0], dat[:,1]/nu, 'r-', mfc='none', linewidth = LW)
    ax1b.plot(dat[:,0], dat[:,2]/nu, 'r-', mfc='none', linewidth = LW)
    
    fnamein = prefix + fbases1[1] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus1[1]
 
    ax1a.plot(dat[:,0], dat[:,1]/nu, 'g-', mfc='none', linewidth = LW)
    ax1b.plot(dat[:,0], dat[:,2]/nu, 'g-', mfc='none', linewidth = LW)
   
    fnamein = prefix + fbases1[2] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus1[2]

    ax1a.plot(dat[:,0], dat[:,1]/nu, 'b-', mfc='none', linewidth = LW)
    ax1b.plot(dat[:,0], dat[:,2]/nu, 'b-', mfc='none', linewidth = LW)


    #-------------------------------------------
    
    suffix = suffix1    
    fnamein = prefix + fbases2[0] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus2[0]

    p=-0.45
    c=100**p
    ax2a.plot(dat[:,0], dat[:,1]/nu*c, 'r-', mfc='none', linewidth = LW)
    ax2b.plot(dat[:,0], dat[:,2]/nu, 'r-', mfc='none', linewidth = LW)
    
    fnamein = prefix + fbases2[1] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus2[1]

    c=250**p
    ax2a.plot(dat[:,0], dat[:,1]/nu*c, 'g-', mfc='none', linewidth = LW)
    ax2b.plot(dat[:,0], dat[:,2]/nu, 'g-', mfc='none', linewidth = LW)
   
    fnamein = prefix + fbases2[2] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus2[2]

    c=500**p
    ax2a.plot(dat[:,0], dat[:,1]/nu*c, 'b-', mfc='none', linewidth = LW)
    ax2b.plot(dat[:,0], dat[:,2]/nu, 'b-', mfc='none', linewidth = LW)


    #-------------------------------------------
    
    suffix = suffix1    
    fnamein = prefix + fbases3[0] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus3[0]
   
    ax3a.plot(dat[:,0], dat[:,1]/nu, 'r-', mfc='none', linewidth = LW)
    ax3b.plot(dat[:,0], dat[:,2]/nu, 'r-', mfc='none', linewidth = LW)
    
    fnamein = prefix + fbases3[1] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus3[1]
 
    ax3a.plot(dat[:,0], dat[:,1]/nu, 'g-', mfc='none', linewidth = LW)
    ax3b.plot(dat[:,0], dat[:,2]/nu, 'g-', mfc='none', linewidth = LW)
   
    fnamein = prefix + fbases3[2] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus3[2]

    ax3a.plot(dat[:,0], dat[:,1]/nu, 'b-', mfc='none', linewidth = LW)
    ax3b.plot(dat[:,0], dat[:,2]/nu, 'b-', mfc='none', linewidth = LW)

    
    #-------------------------------------------
    
    suffix = suffix2    
    fnamein = prefix + fbases4[0] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus4[0]
   
    ax4a.plot(dat[:,0], dat[:,1]/nu, 'r-', mfc='none', linewidth = LW)
    ax4b.plot(dat[:,0], dat[:,2]/nu, 'r-', mfc='none', linewidth = LW)
    
    fnamein = prefix + fbases4[1] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus4[1]
 
    ax4a.plot(dat[:,0], dat[:,1]/nu, 'g-', mfc='none', linewidth = LW)
    ax4b.plot(dat[:,0], dat[:,2]/nu, 'g-', mfc='none', linewidth = LW)
   
    fnamein = prefix + fbases4[2] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus4[2]

    ax4a.plot(dat[:,0], dat[:,1]/nu, 'b-', mfc='none', linewidth = LW)
    ax4b.plot(dat[:,0], dat[:,2]/nu, 'b-', mfc='none', linewidth = LW)

  
   
    
    #-------------------------------------------

    #plt.show() #debug

    #plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================
 



#================================================================

