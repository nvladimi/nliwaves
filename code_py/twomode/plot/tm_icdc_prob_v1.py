
#import importlib ; import tm_icdc_prob_v1
#importlib.reload(tm_icdc_prob_v1) ; tm_icdc_prob_v1.all()



def all():


    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    
    outfile = 'tm_icdc_prob_v1.pdf'

    prefix    =  "../POST/"

    
    fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(10,5.1))
    #fig, ax = plt.subplots(ncols=4, nrows=2, tight_layout=True, figsize=(13,6.5))


    MS = 4
    LP = 3
 
    
    #== INVERSE CASCADE, SMALL CHI ==============================================================


    suffixN   = "_dt1_probN.txt";
    suffixPhi = "_dt1_probPhi.txt";
    
    fbases = ("ig001_p5e5", "ig001_p5e4", "ig001_p5e3")
    pp    = np.asarray([5e-5, 5e-4, 5e-3])
    gg    = np.asarray([0.01, 0.01, 0.01])
    
    chis   = gg*gg*gg / (2*pp)
    nus   = pp/gg
    
    
    #-------------------------------------------

    iax = ax[0,0]
    iax.text( 2.5, 0.3, '(a)')
    iax.yaxis.labelpad = LP  


    
    iax.set_yscale('log')
    iax.set_xlim(0,20)
    iax.set_ylim(1e-6,1)


    ticks = (1e-6, 1e-4, 1e-2, 1)
    iax.set_yticks(ticks)
    iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))
    
    iax.set(xlabel='$\\rho_1^2 / \\nu_1, \quad \\rho_2^2 / \\nu_1$', ylabel='probability')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[0])

    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu, 'r^', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[1])
    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu,  'go', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[2])

    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );

    
    x = np.arange(0,20)    
    iax.plot(x, 2*np.exp(-2*x), 'k--', linewidth = 1)
    iax.plot(x, np.exp(-x), 'k--', linewidth = 1)

    iax.legend()
    iax.legend(frameon=False,  handletextpad=0)
 
    #iax.text(12, 2e-5, '$P(\\rho_1^2) >  e^{-\\rho_1^2 / \\nu_2}$', rotation=-40)
    #iax.text( 2, 3e-6, '$P(\\rho_2^2) > 2 e^{-\\rho_2^2 / (2\\nu_2)}$', rotation=-75)
    
    iax.text(15, 3e-5, '$P(\\rho_1^2)$')
    iax.text( 2, 3e-6, '$P(\\rho_2^2)$')


        
    
    #-------------------------------------------

    iax = ax[0,1]
    iax.text( -0.82, 0.94, '(b)')
    iax.yaxis.labelpad = LP  

    
    iax.set_xlim(-1, 1)
    iax.set_ylim(0.35, 1.0)
    iax.set(xlabel='$\phi / \pi$')
    iax.grid(axis='x')
    ticks = (0.4, 0.6, 0.8, 1)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])

    
    

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    iax.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    iax.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    iax.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1);
 




    #== DIRECT CASCADE, SMALL CHI ==============================================================

    
    suffixN   = "_dt1_probN.txt";
    suffixPhi = "_dt1_probPhi.txt";
    
    fbases = ("dg001_p1e4", "dg001_p1e3", "dg001_p1e2")
    pp    = np.asarray([1e-4, 1e-3, 1e-2])
    gg    = np.asarray([0.01, 0.01, 0.01])
    
    chis   = gg*gg*gg / pp
    nus   = pp/(4*gg)

    
    #-------------------------------------------

    iax = ax[0,2]
    iax.text( 4, 0.3, '(c)')
    iax.yaxis.labelpad = LP  


    
    iax.set_yscale('log')
    iax.set_xlim(0,35)
    iax.set_ylim(1e-6,1)

    ticks = (1e-6, 1e-4, 1e-2, 1)
    iax.set_yticks(ticks)
    iax.set_yticklabels(("$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$1$"))

    
    iax.set(xlabel='$\\rho_1^2 / \\nu_2, \quad \\rho_2^2 / \\nu_2$')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[0])

    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu, 'r^', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[1])
    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu,  'go', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0e}'.format(chis[2])

    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );

    
    x = np.arange(0,35)    
    iax.plot(x, 0.5*np.exp(-0.5*x), 'k--', linewidth = 1)
    iax.plot(x, np.exp(-x), 'k--', linewidth = 1)

    #iax.legend(frameon=False)
 
    #iax.text(11.0, 2e-6, '$P(\\rho_1^2) \sim \\frac{1}{2} e^{-\\rho_1^2 / (2\\nu_2)}$', rotation=-55)
    #iax.text( 3.08, 1e-5, '$P(\\rho_2^2) \sim e^{-\\rho_2^2 / \\nu_2}$', rotation=-75)

    iax.text(23.0, 1e-4, '$P(\\rho_1^2)$')
    iax.text( 3.5, 2e-6, '$P(\\rho_2^2)$')

    
    #-------------------------------------------

    iax = ax[0,3]
    iax.text( -0.8, 0.62, '(d)')
    iax.yaxis.labelpad = LP  

    
    iax.set_xlim(-1, 1)
    iax.set(xlabel='$\phi / \pi$')
    iax.grid(axis="x")

    iax.set_ylim(0.3, 0.65)
    ticks = (0.3, 0.4, 0.5, 0.6)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])

    

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    iax.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    iax.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    iax.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1);


    
#== INVERSE CASCADE, LARGE CHI ==============================================================
        
    suffixN   = "_dt2_probN.txt";
    suffixPhi = "_dt2_probPhi.txt";
    
    fbases = ("ig100_p5e3", "ig100_p2e3", "ig100_p1e3")
    pp    = np.asarray([5e-3, 2e-3, 1e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    
    chis   = gg*gg*gg / (2*pp)
    nus   = pp/gg
    

    iax = ax[1,0]
    iax.text(100, 0.20, '(e)')
    iax.yaxis.labelpad = LP  

    
    #-------------------------------------------
    
    iax.set_yscale('log')
    iax.set_xlim(0,1000)
    iax.set_ylim(1e-8,1)

    iax.set(xlabel='$\\rho_1^2 / \\nu_1, \quad \\rho_2^2 / \\nu_1$', ylabel="probability")
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu, 'r^', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])
    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu,  'go', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    
    iax.plot(dat[:,0]/nu, dat[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );

    
    x = np.arange(0,1000,10)    
    iax.plot(x, 0.0030*np.exp(-0.0030*x), 'k--', linewidth = 1)
    iax.plot(x, 0.0060*np.exp(-0.0060*x), 'k--', linewidth = 1)
    iax.plot(x, 0.0130*np.exp(-0.0130*x), 'k--', linewidth = 1)

    iax.legend(frameon=False, handletextpad=0)
 
    iax.text(200, 2.5e-3, '$P(\\rho_1^2)$')
    iax.text(80, 5e-7, '$P(\\rho_2^2)$')

    #iax.text(800, 5e-4, '$e^{-0.003 x}$')
    #iax.text(660, 2e-5, '$e^{-0.006 x}$')
    #iax.text(650, 4e-6, '$e^{-0.013 x}$')

    
    #-------------------------------------------

    iax = ax[1,1]
    iax.text(-0.82, 54, '(f)')
    iax.axvline(x=-0.5, color='gray', lw = 0.5)
    iax.axhline(color='gray', lw = 0.5)
    iax.yaxis.labelpad = LP  

    iax.set_ylim(-5, 60)
    ticks = (0, 20, 40, 60)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])

    
    
    iax.set_xlim(-1, 1)
    iax.set(xlabel='$\phi / \pi$')
    #iax.grid(True)

    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    iax.plot(dat[:,0], dat[:,1], 'r-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    iax.plot(dat[:,0], dat[:,1], 'g-', mfc='none', markersize = MS, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    iax.plot(dat[:,0], dat[:,1], 'b-', mfc='none', markersize = MS, linewidth = 1);


    #-------------------------------------------

    ax3 = inset_axes(iax, width="47%", height="47%", borderpad=0.5)

    ax3.set_xlim(-0.56, -0.44)
    ax3.axvline(x=-0.5, color='gray', lw = 0.5)
    ax3.axhline(color='gray', lw = 0.5)

    ax3.set_ylim(-5, 60)
    ticks = (0, 20, 40, 60)
    ax3.set_yticks(ticks)
    ax3.set_yticklabels(["$%g$"%y for y in ticks], fontsize=8)

    ticks = (-0.55, -0.5, -0.45)
    ax3.set_xticks(ticks)
    ax3.set_xticklabels(("$-0.55$", "$-0.5$", "$-0.45$"), fontsize=8)
   

    
    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1

    
    ax3.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS/2, linewidth = 1);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax3.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS/2, linewidth = 1);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    
    ax3.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS/2, linewidth = 1);


    
#== DIRECT CASCADE, LARGE CHI ==============================================================


    suffixN   = "_dt1_probN.txt";
    suffixNx4 = "_dt1x4_probN.txt";
    suffixPhi = "_dt1_probPhi.txt";
    
    fbases = ("dg100_p1e2", "dg100_p4e3", "dg100_p2e3")
    
    pp    = np.asarray([1e-2, 4e-3, 2e-3])
    gg    = np.asarray([1.00, 1.00, 1.00])
    
    chis   = gg*gg*gg / pp
    nus   = pp/(4*gg)
    
    ax3 = inset_axes(ax[1,2], width="47%", height="47%", borderpad=0.5)

    
    #-------------------------------------------

    iax = ax[1,2]
    iax.text( 30, 0.2, '(g)')
    iax.yaxis.labelpad = LP  

    
    iax.set_yscale('log')
    iax.set_xlim(0,300)
    iax.set_ylim(1e-8,1)
    
    ax3.set_yscale('log')
    ax3.set_xlim(0,20)
    ax3.set_ylim(1e-6,1)

    
    
    iax.set(xlabel='$\\rho_1^2 / \\nu_2, \quad \\rho_2^2 / \\nu_2$')
 
        
    fnamein = prefix + fbases[0] + suffixN
    nu = nus[0]
    dat = np.loadtxt(fnamein);
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    fnamein = prefix + fbases[0] + suffixNx4
    dat4 = np.loadtxt(fnamein);
    
    
    iax.plot(dat4[:,0]/nu, dat4[:,1]*nu, 'r^', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS,  label = lbl );
    ax3.plot(dat[:,0]/nu, dat[:,2]*nu, 'r^', mfc='none', markersize = MS/2, lw=0.5);


    fnamein = prefix + fbases[1] + suffixN
    nu = nus[1]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])

    fnamein = prefix + fbases[1] + suffixNx4
    dat4 = np.loadtxt(fnamein);
    
    iax.plot(dat4[:,0]/nu, dat4[:,1]*nu,  'go', mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu,  'go', mfc='none', markersize = MS,  label = lbl );
    ax3.plot(dat[:,0]/nu, dat[:,2]*nu, 'go', mfc='none', markersize = MS/2, lw=0.5);

    
    fnamein = prefix + fbases[2] + suffixN
    nu = nus[2]
    dat = np.loadtxt(fnamein)
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    fnamein = prefix + fbases[2] + suffixNx4
    dat4 = np.loadtxt(fnamein);
    
    iax.plot(dat4[:,0]/nu, dat4[:,1]*nu, 'bs',  mfc='none', markersize = MS);
    iax.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs',  mfc='none', markersize = MS, label = lbl );
    ax3.plot(dat[:,0]/nu, dat[:,2]*nu, 'bs', mfc='none', markersize = MS/2, lw=0.5);

    
    x = np.arange(0.1,35,0.1)
    f = 0.5*np.exp(-x*0.5)/np.sqrt(2*np.pi*x)
    iax.plot(x, f, 'k--', linewidth = 1)
    ax3.plot(x, f, 'k-', linewidth = 1)

    
    #iax.legend(frameon=False, loc="upper left")
 
    iax.text(240, 7e-7, '$P(\\rho_1^2)$', rotation = 0)
    iax.text( 35, 8e-8, '$P(\\rho_2^2)$', rotation = 0)
    ax3.text( 8, 5e-2, '$P(\\rho_2^2)$', rotation = 0)
   # ax3.text( 8, 1e-2, '$\\frac{1}{\sqrt{2 \pi \\tilde{n}_2}} e^{- \\tilde{n}_2/2}$')
   # ax3.text( 3, 1e-5, '$\\tilde{n}_2 \equiv \\rho_2^2 / \\nu_2$')

    ax3.set_xlim(0, 20)
    ticks = (5, 10, 15, 20)
    ax3.set_xticks(ticks)
    ax3.set_xticklabels(("", "$10$", "$15$", "$20$"), fontsize=8)

    ticks = (1e-4, 1e-2, 1)
    ax3.set_yticks(ticks)
    ax3.set_yticklabels(("", "$10^{-2}$", "$1$"), fontsize=8)

   

    x = np.arange(1,150)
    c = 1/100.
    f = (1 - c*x) * np.exp( -c/8 *x*x + c*c/12 *x*x*x - 0*c*c*c/16 *x*x*x*x  +
                            0*c*c*c*c/20 *x*x*x*x*x + 0*c*c*c*c*c/24 *x*x*x*x*x*x)
    f = f/sum(f)
    iax.plot(x, f, 'r--', linewidth = 1)


    x = np.arange(1,250)
    c = 1/250.
    f = (1 - c*x) * np.exp( -c/8 *x*x + c*c/12 *x*x*x  - 0*c*c*c/16 *x*x*x*x  +
                            0*c*c*c*c/20 *x*x*x*x*x + 0*c*c*c*c*c/24 *x*x*x*x*x*x) 
    f = f/sum(f)
    iax.plot(x, f, 'g--', linewidth = 1)

    x = np.arange(1,300)
    c = 1/500.
    f = (1 - c*x) * np.exp( -c/8 *x*x + c*c/12 *x*x*x  - 0*c*c*c/16 *x*x*x*x  +
                            0*c*c*c*c/20 *x*x*x*x*x + 0*c*c*c*c*c/24 *x*x*x*x*x*x) 
    f = f/sum(f)
    iax.plot(x, f, 'b--', linewidth = 1)


    
    
    #-------------------------------------------

    iax = ax[1,3]
    iax.text( -0.8, 3.6, '(h)')
    iax.yaxis.labelpad = LP  

    
    iax.set_xlim(-1, 1)
    iax.set_ylim(-0.3, 4)

    
    iax.set(xlabel='$\phi / \pi$')
    #iax.grid(True)
    iax.axvline(x=0.5, color='gray', lw = 0.5)
    iax.axhline(color='gray', lw = 0.5)
    
    fnamein = prefix + fbases[0] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[0])

    
    iax.plot(dat[:,0], dat[:,1], 'r^-', mfc='none', markersize = MS, linewidth = 1, label=lbl);

    fnamein = prefix + fbases[1] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[1])

    
    iax.plot(dat[:,0], dat[:,1], 'go-', mfc='none', markersize = MS, linewidth = 1, label=lbl);

    fnamein = prefix + fbases[2] + suffixPhi
    dat = np.loadtxt(fnamein)
    dat = np.vstack((dat, dat[0,:]))
    dat[-1,0] = 1
    lbl = "$\chi = $" + ' {:1.0f}'.format(chis[2])

    
    iax.plot(dat[:,0], dat[:,1], 'bs-', mfc='none', markersize = MS, linewidth = 1, label=lbl);

    #iax.legend(frameon=True)






    
#================================================================

    plt.subplots_adjust(left=0.06, right=0.99, top=0.98, bottom=0.09, hspace=0.3, wspace=0.25)

    #plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
    
#================================================================
