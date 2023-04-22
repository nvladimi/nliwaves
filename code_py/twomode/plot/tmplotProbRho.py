
#import importlib ; import tmplotProbRho
#importlib.reload(tmplotProbRho) ; tmplotProbRho.probAll()



def probAll():

    #probRhoICsmallChi()
    probRhoIClargeChi()


#================================================================

    
#================================================================
    


def probRhoICsmallChi():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D


    outfile = 'PLOT/probRhoICsmallChi.pdf'
    
    prefix    =  "POST/"
    suffix1   = "_dt1_probRho1.txt";
    suffix2   = "_dt1_probRho2.txt";
    
    fbases = ("ig001_p5e5", "ig001_p5e4", "ig001_p5e3")
    pp    = np.asarray([5e-5, 5e-4, 5e-3])
    gg    = np.asarray([0.01, 0.01, 0.01])
    
    chis   = gg*gg*gg / (2*pp)
    nus   = pp/gg
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    MS = 4
    LW1 = 1.5
    
    #-------------------------------------------
    
    ax1.set_yscale('log')
    ax1.set_xlim(0,4)
    ax1.set_ylim(1e-6,100)

    ax1.set(xlabel='$\\rho_1$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffix1
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[0])

    
    ax1.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, label = lbl );


    fnamein = prefix + fbases[1] + suffix1
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[1])
    
    ax1.plot(dat[:,0], dat[:,1],  'go-', mfc='none', markersize = MS, label = lbl );
 
    
    fnamein = prefix + fbases[2] + suffix1
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[2])
    
    ax1.plot(dat[:,0], dat[:,1], 'bs-',  mfc='none', markersize = MS, label = lbl );


    if True:
 
        dx=0.01;
        x = np.arange(dx, 4, dx)

        g = gg[0]
        ggg = g*g*g
    
        p = -1 + 4*pp[0]/ggg
        print(p)
        y = x**p * np.exp( -4/g * x )
        y = y / np.sum(y)/dx    
        ax1.plot(x, y, 'r--', linewidth = 1)
    

        p = -1 + 4*pp[1]/ggg
        print(p)
        y = x**p * np.exp( -4/g * x )
        y = y / np.sum(y)/dx    
        ax1.plot(x, y, 'g--', linewidth = 1)

        p = -1 + 4*pp[2]/ggg
        print(p)
        y = x**p * np.exp( -4/g * x )
        y = y / np.sum(y)/dx    
        ax1.plot(x, y, 'b--', linewidth = 1)


        
    ax1.legend(frameon=False)
 
    # ax1.text(12, 2e-5, '$P(\\rho_1^2) >  e^{-\\rho_1^2 / \\nu_2}$', rotation=-40)
    # ax1.text( 3, 3e-6, '$P(\\rho_2^2) > 2 e^{-\\rho_2^2 / (2\\nu_2)}$', rotation=-75)
 
    
    #-------------------------------------------

    ax2.set_yscale('log')
    ax2.set_xlim(0,2.5)
    ax2.set_ylim(1e-6,10)

    ax2.set(xlabel='$\\rho_2$', ylabel='probability')

    fnamein = prefix + fbases[0] + suffix2
    dat = np.loadtxt(fnamein)
    
    ax2.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = LW1);

    fnamein = prefix + fbases[1] + suffix2
    dat = np.loadtxt(fnamein)
    
    ax2.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = LW1);

    fnamein = prefix + fbases[2] + suffix2
    dat = np.loadtxt(fnamein)
    
    ax2.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = LW1);


    dx=0.01;
    x = np.arange(0, 4, dx)
   
    y = x * np.exp( -4/gg[0] * x )
    y = y / np.sum(y)/dx    
    ax2.plot(x, y, 'k--', linewidth = 1)

    
    #ax2.legend(frameon=False)

    
    #-------------------------------------------

    
    #plt.show() #debug

    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================



    
#================================================================
    


def probRhoIClargeChi():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    #from matplotlib.lines import Line2D


    outfile = 'PLOT/probRhoIClargeChi.pdf'
    
    prefix    =  "POST/"
    suffix1   = "_dt2_probRho1.txt";
    suffix2 = "_dt2_probRho2.txt";
    
    fbases = ("ig100_p5e3", "ig100_p2e3", "ig100_p1e3")
    pp    = np.asarray([5e-3, 2e-3, 1e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    
    chis   = gg*gg*gg / (2*pp)
    nus   = pp/gg
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    MS = 4
    LW1 = 1.5
    
    #-------------------------------------------

    ax1.set_yscale('log')
    ax1.set_xlim(0,2.5)
    ax1.set_ylim(1e-6,10)

    ax1.set(xlabel='$\\rho_1$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffix1
    nu = 1
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    
    ax1.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffix1
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])
    
    ax1.plot(dat[:,0], dat[:,1],  'go-', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffix1
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    ax1.plot(dat[:,0], dat[:,1],  'bs-', mfc='none', markersize = MS,  label = lbl );


    
    #-------------------------------------

    dx = 0.01
    x = np.arange(dx,3, dx)

    g = gg[0]
    ggg = g*g*g

    lbl = "$ \\rho^{-1 + 2/\\chi} \exp ( - \\frac{2V^2}{\gamma^2} \\rho^2 )$" 
     
    c=1
    
    p = -1 + 2/chis[0]*c
    y = x**p * np.exp( -2/g/g * x*x)
    y = y / np.sum(y)/dx    
    ax1.plot(x, y, 'r--', linewidth = 1)
    
    p = -1 + 2/chis[1]*c
    y = x**p * np.exp( -2/g/g * x*x)
    y = y / np.sum(y)/dx    
    ax1.plot(x, y, 'g--', linewidth = 1)

    p = -1 + 2/chis[2]*c
    y = x**p * np.exp( -2/g/g * x*x)
    y = y / np.sum(y)/dx    
    ax1.plot(x, y, 'b--', linewidth = 1, label = lbl)


    #lbl = "$ \exp ( - \\frac{4V}{\gamma} \\rho )$" 
    
    #y = np.exp( -4/g * x)
    #y = y / np.sum(y)/dx  * 0.1   
    #ax1.plot(x, y, 'k-', linewidth = 0.5, label = lbl)

    lbl = "$ \\rho^{-1} \exp ( - 2 \\rho^2 )$" 

    y = 1/x * np.exp( -2*x*x)
    y = y / np.sum(y)/dx    
    ax1.plot(x, y, 'k-', linewidth = 0.5, label = lbl)

    
    
    ax1.legend(frameon=False)


    #======================================================================

    ax2.set_yscale('log')
    ax2.set_xlim(0,2.5)
    ax2.set_ylim(1e-4,5)

    ax2.set(xlabel='$\\rho_2$', ylabel='probability')

    fnamein = prefix + fbases[0] + suffix2
    dat = np.loadtxt(fnamein)

    ax2.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = LW1);

    fnamein = prefix + fbases[1] + suffix2
    dat = np.loadtxt(fnamein)
    
    ax2.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = LW1);

    fnamein = prefix + fbases[2] + suffix2
    dat = np.loadtxt(fnamein)

    ax2.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = LW1);
    


    #-------------------------------------------

    
    lbl = "$ \\rho \exp (- \\frac{4V}{\gamma} \\rho )$" 

    dx=0.01;
    x = np.arange(0,3, dx)

    g=gg[0]
    
    y = x * np.exp( -4 * x/g)
    y = y / np.sum(y)/dx    
    ax2.plot(x, y, 'k-', linewidth = 0.5,  label = lbl)

    
    lbl = "$ \\rho^{3/2} \exp (- 4 (\\frac{V}{\gamma} \\rho)^{3/2} )$" 

    dx=0.01;
    x = np.arange(0,3, dx)

    q=1.5
    y = x**q * np.exp( -4 * (x/g)**q )
    y = y / np.sum(y)/dx    
    ax2.plot(x, y, 'k--', linewidth = 1,  label = lbl)

    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])


    ax2.legend(frameon=False)


   #-------------------------------------------

    
    #plt.show() #debug

    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================



