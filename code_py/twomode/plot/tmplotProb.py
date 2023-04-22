
#import importlib ; import tmplotProb
#importlib.reload(tmplotProb) ; tmplot.probAll()



def probAll():

    probDCsmallChi()
    probICsmallChi()
    probDClargeChi()
    probIClargeChi()


#================================================================

def probDCsmallChi():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D

    outfile = 'PLOT/probDCsmallChi.pdf'

    prefix    =  "POST/"
    suffixN   = "_dt1_probN.txt";
    suffixPhi = "_dt1_probPhi.txt";
    
    fbases = ("dg001_p1e4", "dg001_p1e3", "dg001_p1e2")
    pp    = np.asarray([1e-4, 1e-3, 1e-2])
    gg    = np.asarray([0.01, 0.01, 0.01])
    
    chis   = gg*gg*gg / pp
    nus   = pp/(4*gg)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    MS = 4
    
    #-------------------------------------------
    
    ax1.set_yscale('log')
    ax1.set_xlim(0,35)
    ax1.set_ylim(1e-6,1)

    ax1.set(xlabel='$\\rho_1^2 / \\nu_2, \quad \\rho_2^2 / \\nu_2$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[0])

    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu, 'r^', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[1])
    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu,  'go', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[2])

    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );

    
    x = np.arange(0,35)    
    ax1.plot(x, 0.5*np.exp(-0.5*x), 'k--', linewidth = 1)
    ax1.plot(x, np.exp(-x), 'k--', linewidth = 1)

    ax1.legend()
    ax1.legend(frameon=False)
 
    ax1.text(13.5, 2e-6, '$P(\\rho_1^2) \sim \\frac{1}{2} e^{-\\rho_1^2 / (2\\nu_2)}$', rotation=-55)
    ax1.text( 4.0, 1e-5, '$P(\\rho_2^2) \sim e^{-\\rho_2^2 / \\nu_2}$', rotation=-75)
 
    
    #-------------------------------------------

    ax2.set_xlim(-1, 1)
    ax2.set(xlabel='$\phi / \pi$', ylabel='probability')
    ax2.grid(True)

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    ax2.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax2.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax2.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1);

    
    #-------------------------------------------

    #plt.show() #debug

    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================
    


def probICsmallChi():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D


    outfile = 'PLOT/probICsmallChi.pdf'
    
    prefix    =  "POST/"
    suffixN   = "_dt1_probN.txt";
    suffixPhi = "_dt1_probPhi.txt";
    
    fbases = ("ig001_p5e5", "ig001_p5e4", "ig001_p5e3")
    pp    = np.asarray([5e-5, 5e-4, 5e-3])
    gg    = np.asarray([0.01, 0.01, 0.01])
    
    chis   = gg*gg*gg / (2*pp)
    nus   = pp/gg
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    MS = 4
    
    #-------------------------------------------
    
    ax1.set_yscale('log')
    ax1.set_xlim(0,20)
    ax1.set_ylim(1e-6,1)

    
    ax1.set(xlabel='$\\rho_1^2 / \\nu_1, \quad \\rho_2^2 / \\nu_1$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[0])

    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu, 'r^', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[1])
    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu,  'go', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[2])

    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );

    
    x = np.arange(0,20)    
    ax1.plot(x, 2*np.exp(-2*x), 'k--', linewidth = 1)
    ax1.plot(x, np.exp(-x), 'k--', linewidth = 1)

    ax1.legend()
    ax1.legend(frameon=False)
 
    ax1.text(12, 2e-5, '$P(\\rho_1^2) >  e^{-\\rho_1^2 / \\nu_2}$', rotation=-40)
    ax1.text( 3, 3e-6, '$P(\\rho_2^2) > 2 e^{-\\rho_2^2 / (2\\nu_2)}$', rotation=-75)
 
    
    #-------------------------------------------

    ax2.set_xlim(-1, 1)
    ax2.set(xlabel='$\phi / \pi$', ylabel='probability')
    ax2.grid(True)

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    ax2.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax2.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax2.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1);

    
    #-------------------------------------------

    #plt.show() #debug

    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================



#================================================================

def probDClargeChi():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D

    outfile = 'PLOT/probDClargeChi.pdf'

    prefix    =  "POST/"
    suffixN   = "_dt1_probN.txt";
    suffixNx4 = "_dt1x4_probN.txt";
    suffixPhi = "_dt1_probPhi.txt";
    
    fbases = ("dg100_p1e2", "dg100_p4e3", "dg100_p2e3")
    
    pp    = np.asarray([1e-2, 4e-3, 2e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    
    chis   = gg*gg*gg / pp
    nus   = pp/(4*gg)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    ax3 = inset_axes(ax1, width="47%", height="47%", borderpad=0.5)


    MS = 4
    
    #-------------------------------------------
    
    ax1.set_yscale('log')
    ax1.set_xlim(0,300)
    ax1.set_ylim(1e-8,1)
    
    ax3.set_yscale('log')
    ax3.set_xlim(0,20)
    ax3.set_ylim(1e-6,1)

    
    
    ax1.set(xlabel='$\\rho_1^2 / \\nu_2, \quad \\rho_2^2 / \\nu_2$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein);
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    fnamein = prefix + fbases[0] + suffixNx4
    dat4 = np.loadtxt(fnamein);
    
    
    ax1.plot(dat4[:,0]/nu, dat4[:,1]*nu, 'r^', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );
    ax3.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS);


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])

    fnamein = prefix + fbases[1] + suffixNx4
    dat4 = np.loadtxt(fnamein);
    
    ax1.plot(dat4[:,0]/nu, dat4[:,1]*nu,  'go', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );
    ax3.plot(dat[:,0]/nu, dat[:,2]*nu, 'go', mfc='none', markersize = MS);

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    fnamein = prefix + fbases[2] + suffixNx4
    dat4 = np.loadtxt(fnamein);
    
    ax1.plot(dat4[:,0]/nu, dat4[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );
    ax3.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs', mfc='none', markersize = MS);

    
    x = np.arange(0.1,35,0.1)
    f = 0.5*np.exp(-x*0.5)/np.sqrt(2*np.pi*x)
    ax1.plot(x, f, 'k--', linewidth = 1)
    ax3.plot(x, f, 'k-', linewidth = 1)

    
    #ax1.legend(frameon=False, loc="upper left")
 
    ax1.text(250, 1e-6, '$P(\\rho_1^2)$', rotation = 0)
    ax1.text( 35, 8e-8, '$P(\\rho_2^2)$', rotation = 0)
    ax3.text( 8, 1e-1, '$P(\\rho_2^2)$', rotation = 0)
    ax3.text( 8, 1e-2, '$\\frac{1}{\sqrt{2 \pi \\tilde{n}_2}} e^{- \\tilde{n}_2/2}$')
    ax3.text( 3, 1e-5, '$\\tilde{n}_2 \equiv \\rho_2^2 / \\nu_2$')


    x = np.arange(1,150)
    c = 1/100.
    f = (1 - c*x) * np.exp( -c/8 *x*x + c*c/12 *x*x*x - 0*c*c*c/16 *x*x*x*x  +
                            0*c*c*c*c/20 *x*x*x*x*x + 0*c*c*c*c*c/24 *x*x*x*x*x*x)
    f = f/sum(f)
    ax1.plot(x, f, 'r--', linewidth = 1)


    x = np.arange(1,250)
    c = 1/250.
    f = (1 - c*x) * np.exp( -c/8 *x*x + c*c/12 *x*x*x  - 0*c*c*c/16 *x*x*x*x  +
                            0*c*c*c*c/20 *x*x*x*x*x + 0*c*c*c*c*c/24 *x*x*x*x*x*x) 
    f = f/sum(f)
    ax1.plot(x, f, 'g--', linewidth = 1)

    x = np.arange(1,300)
    c = 1/500.
    f = (1 - c*x) * np.exp( -c/8 *x*x + c*c/12 *x*x*x  - 0*c*c*c/16 *x*x*x*x  +
                            0*c*c*c*c/20 *x*x*x*x*x + 0*c*c*c*c*c/24 *x*x*x*x*x*x) 
    f = f/sum(f)
    ax1.plot(x, f, 'b--', linewidth = 1)


    
    
    #-------------------------------------------

    ax2.set_xlim(-1, 1)
    ax2.set(xlabel='$\phi / \pi$', ylabel='probability')
    ax2.grid(True)

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    
    ax2.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1, label=lbl);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])

    
    ax2.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1, label=lbl);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    
    ax2.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1, label=lbl);

    ax2.legend(frameon=True)

    #-------------------------------------------



    
    
    #-------------------------------------------

    #plt.show() #debug

    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================
    


def probIClargeChi():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D


    outfile = 'PLOT/probIClargeChi.pdf'
    
    prefix    =  "POST/"
    suffixN   = "_dt2_probN.txt";
    suffixPhi = "_dt2_probPhi.txt";
    
    fbases = ("ig100_p5e3", "ig100_p2e3", "ig100_p1e3")
    pp    = np.asarray([5e-3, 2e-3, 1e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    
    chis   = gg*gg*gg / (2*pp)
    nus   = pp/gg
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    MS = 4
    
    #-------------------------------------------
    
    ax1.set_yscale('log')
    ax1.set_xlim(0,1000)
    ax1.set_ylim(1e-8,1)

    ax1.set(xlabel='$\\rho_1^2 / \\nu_1, \quad \\rho_2^2 / \\nu_1$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu, 'r^', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])
    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu,  'go', mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    
    ax1.plot(dat[:,0]/nu, dat[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    ax1.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );

    
    x = np.arange(0,1000,10)    
    ax1.plot(x, 0.0030*np.exp(-0.0030*x), 'k--', linewidth = 1)
    ax1.plot(x, 0.0060*np.exp(-0.0060*x), 'k--', linewidth = 1)
    ax1.plot(x, 0.0130*np.exp(-0.0130*x), 'k--', linewidth = 1)

    ax1.legend()
    ax1.legend(frameon=False)
 
    ax1.text(90, 7e-3, '$P(\\rho_1^2)$')
    ax1.text(90, 2e-6, '$P(\\rho_2^2)$')

    ax1.text(800, 5e-4, '$e^{-0.003 x}$')
    ax1.text(660, 2e-5, '$e^{-0.006 x}$')
    ax1.text(650, 4e-6, '$e^{-0.013 x}$')

    
    #-------------------------------------------

    ax2.set_xlim(-1, 1)
    ax2.set(xlabel='$\phi / \pi$', ylabel='probability')
    ax2.grid(True)

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    ax2.plot(dat[:,0], dat[:,1], 'r-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax2.plot(dat[:,0], dat[:,1], 'g-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax2.plot(dat[:,0], dat[:,1], 'b-', mfc='none', markersize = MS, linewidth = 1);


    #-------------------------------------------

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    ax3 = inset_axes(ax2, width="50%", height="50%", borderpad=0.5)

    ax3.set_xlim(-0.56, -0.44)
    #ax3.set(xlabel='$\phi / \pi$', ylabel='probability')
    ax3.grid(True)

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    ax3.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax3.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax3.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1);


    
    #plt.show() #debug

    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================






#================================================================




#================================================================

