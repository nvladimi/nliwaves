
#import importlib ; import tmplotBursts
#importlib.reload(tmplotBursts) ; tmplotBursts.Bursts()



#================================================================

def Bursts():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D

    outfile = 'PLOT/pdfBursts.pdf'

    prefix   =  "POST/"
    suffix1  = "_dt1_factor2_bursts.txt"
    suffix2  = "_dt2_factor2_bursts.txt"


    
    #-------------------------------------------
    
    fbases = ("ig100_p5e3", "ig100_p2e3", "ig100_p1e3")
    pp    = np.asarray([5e-3, 2e-3, 1e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    chis  = gg*gg*gg / (2*pp)
    nus   = pp/gg

    tbinsize = np.asarray([50, 100, 500])
    tnumbins = np.asarray([20, 20, 10])

    rbinsize = np.asarray([0.1, 0.1, 0.2])
    rnumbins = np.asarray([20, 20, 10])

    #-------------------------------------------
    
    fig, f_axes = plt.subplots(ncols=2, nrows=3, constrained_layout=True, figsize=(8,10))

    for i in (0,1,2):
        f_axes[i,0].set(xlabel='$\Delta t_b$', ylabel='probability', yscale='log')
 
    for i in (0,1,2):
        f_axes[i,1].set(xlabel='$\\rho_b$', ylabel='probability', yscale='log')

    
    #ax1.set_xlim(0,35)
    #ax1.set_ylim(1e-6,1)

    LW1 = 1.0
    LW2 = 0.5
    MS1 = 6
    MS2 = 3


    ax1a = f_axes[0,0] ; ax1a.set_xlim(0,1000);  ax1a.set_ylim(1e-5,1e-2)
    ax1b = f_axes[0,1] ; ax1b.set_xlim(0,2);     ax1b.set_ylim(1e-2,2)
    ax2a = f_axes[1,0] ; ax2a.set_xlim(0,2000);  ax2a.set_ylim(1e-5,1e-2)
    ax2b = f_axes[1,1] ; ax2b.set_xlim(0,2);     ax2b.set_ylim(1e-2, 2)
    ax3a = f_axes[2,0] ; ax3a.set_xlim(0,5000);  ax3a.set_ylim(1e-5,1e-2)
    ax3b = f_axes[2,1] ; ax3b.set_xlim(0,2);     ax3b.set_ylim(1e-2, 2)


    
#================================================================

    irun = 0
    suffix = suffix1

    
    fnamein = prefix + fbases[irun] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus[irun]
 

    tb = dat[1:,0] - dat[:-1,0]
    rb = dat[1:,1]

    binedges = np.arange(tnumbins[irun])*tbinsize[irun]

    (p, x) = np.histogram(tb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
    ax1a.plot(x, p, 'ro-', mfc='none', linewidth = LW1, markersize = MS1)
     
    binedges = np.arange(rnumbins[irun])*rbinsize[irun]

    (p, x) = np.histogram(rb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
   
    ax1b.plot(x, p, 'ro-', mfc='none', linewidth = LW1, markersize = MS1)
    
    #-------------------------------------------
    
    irun = 0
    suffix = suffix2

    
    fnamein = prefix + fbases[irun] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus[irun]
 

    tb = dat[1:,0] - dat[:-1,0]
    rb = dat[1:,1]

    binedges = np.arange(tnumbins[irun])*tbinsize[irun]

    (p, x) = np.histogram(tb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
    ax1a.plot(x, p, 'bo-', mfc='none', linewidth = LW2, markersize = MS2)

    
    binedges = np.arange(rnumbins[irun])*rbinsize[irun]

    (p, x) = np.histogram(rb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
   
    ax1b.plot(x, p, 'bo-', mfc='none', linewidth = LW2, markersize = MS2)
    

#================================================================

    irun = 1
    suffix = suffix1

    
    fnamein = prefix + fbases[irun] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus[irun]
 

    tb = dat[1:,0] - dat[:-1,0]
    rb = dat[1:,1]

    binedges = np.arange(tnumbins[irun])*tbinsize[irun]

    (p, x) = np.histogram(tb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
    
    ax2a.plot(x, p, 'ro-', mfc='none', linewidth = LW1, markersize = MS1)

    
    binedges = np.arange(rnumbins[irun])*rbinsize[irun]

    (p, x) = np.histogram(rb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
   
    ax2b.plot(x, p, 'ro-', mfc='none', linewidth = LW1, markersize = MS1)
    
    #-------------------------------------------
    
    irun = 1
    suffix = suffix2

    
    fnamein = prefix + fbases[irun] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus[irun]
 

    tb = dat[1:,0] - dat[:-1,0]
    rb = dat[1:,1]

    binedges = np.arange(tnumbins[irun])*tbinsize[irun]

    (p, x) = np.histogram(tb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
    
    ax2a.plot(x, p, 'bo-', mfc='none', linewidth = LW2, markersize = MS2)

    
    binedges = np.arange(rnumbins[irun])*rbinsize[irun]

    (p, x) = np.histogram(rb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
   
    ax2b.plot(x, p, 'bo-', mfc='none', linewidth = LW2, markersize = MS2)
    
 #================================================================

    irun = 2
    suffix = suffix2

    
    fnamein = prefix + fbases[irun] + suffix
    dat = np.loadtxt(fnamein)
    nu = nus[irun]
 

    tb = dat[1:,0] - dat[:-1,0]
    rb = dat[1:,1]

    binedges = np.arange(tnumbins[irun])*tbinsize[irun]

    (p, x) = np.histogram(tb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
    
    ax3a.plot(x, p, 'bo-', mfc='none', linewidth = LW2, markersize = MS2)

    
    binedges = np.arange(rnumbins[irun])*rbinsize[irun]

    (p, x) = np.histogram(rb, binedges, density=True)
    x = 0.5 * (binedges[:-1] + binedges[1:])
   
    ax3b.plot(x, p, 'bo-', mfc='none', linewidth = LW2, markersize = MS2)

  #================================================================


    ax1a.text( 800, 6e-3, '$\chi = 100$')
    ax1b.text(1.6,  0.5, ' $\chi = 100$')
    ax2a.text(1600, 6e-3, '$\chi = 200$')
    ax2b.text(1.6,  1.3,  '$\chi = 200$')
    ax3a.text(4000, 6e-3, '$\chi = 500$')
    ax3b.text(1.6,  1.3,  '$\chi = 500$')


    ax1a.text(  400, 2.0e-4, '$P(t_b) = a e^{-ax}, \quad a = 0.0055$', rotation = -35)
    ax2a.text(  600, 2.5e-4, '$P(t_b) = a e^{-ax}, \quad a = 0.0022$', rotation = -27)
    ax3a.text( 1500, 6.0e-5, '$P(t_b) = a e^{-ax}, \quad a = 0.0011$', rotation = -30)


  
    x = np.arange(0,1100,100)
    a = 0.0055
    ax1a.plot(x, a*np.exp(-a*x), '--k', linewidth = LW2)

    x = np.arange(0,2200,200)
    a = 0.0022
    ax2a.plot(x, a*np.exp(-a*x), '--k', linewidth = LW2)

    x = np.arange(0,5500,500)
    a = 0.0011
    ax3a.plot(x, a*np.exp(-a*x), '--k', linewidth = LW2)

    
    x = np.arange(0,2.1,0.01)
    
    a = 3.5
    lbl = '$a^2 x e^{-ax}, \quad a = 3.5$'
    
    ax1b.plot(x, a*a*x*np.exp(-a*x), '-k', linewidth = LW2, label=lbl)
    ax2b.plot(x, a*a*x*np.exp(-a*x), '-k', linewidth = LW2)
    ax3b.plot(x, a*a*x*np.exp(-a*x), '-k', linewidth = LW2)


    a=2
    lbl = '$2 a x e^{-ax^2}, \quad a = 2$'

    ax1b.plot(x, 2*a*x*np.exp(-a*x*x), '--k', linewidth = LW2, label=lbl)
    ax2b.plot(x, 2*a*x*np.exp(-a*x*x), '--k', linewidth = LW2)
    ax3b.plot(x, 2*a*x*np.exp(-a*x*x), '--k', linewidth = LW2)
 
 
    ax1b.legend(frameon=False)
    ax1b.legend(frameon=False)
    ax1b.legend(frameon=False)


  #================================================================
   

    #plt.show() #debug

    #plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================
 



#================================================================

