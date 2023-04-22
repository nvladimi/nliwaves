#import importlib ; import plotCap512anisotropy
#importlib.reload(plotCap512anisotropy) ; plotCap512anisotropy.All()


import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np
import nlwTools


#=================================================================

def All(every=10):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import nlwTools

    fnameout = 'a1n512_anisotropy.pdf'
    prefix = "../wexac/a1n512_20210123/n512evol_"
    suffix = ".txt"
    
    seeds = ('s50', 's51', 's52', 's53')
    seedcolor = {'s50': 'C0',
                 's51': 'C1',
                 's52': 'C2',
                 's53': 'C3'}
    LW = 0.5
    
    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(12,8),  constrained_layout=True)
    fig.suptitle('Anisotropy of capillary runs with forcing 0.001 on 512x512 grid', fontsize=16)

    setcanvas(ax)

    
    for s in seeds:

        fname = prefix + s + suffix
        
        dat = np.loadtxt(fname)  # (t, nk, nq, pkx, pky, qkx, qky)
                
        t     = dat[:,0]
        nk    = dat[:,1]
        nq    = dat[:,2]
        pkx   = dat[:,3]
        pky   = dat[:,4]
        qkx   = dat[:,5]
        qky   = dat[:,6]


        pmod = np.sqrt(pkx**2 + pky**2)
        pphi = np.arctan2(pky, pkx)/np.pi
        qmod = np.sqrt(qkx**2 + qky**2)
        qphi = np.arctan2(qky, qkx)/np.pi
        
        color = seedcolor[s]
        
        ax[0,0].plot(t, nk,   '-', linewidth = LW, color = color)
        ax[1,0].plot(t, nq,   '-', linewidth = LW, color = color)
        ax[0,1].plot(t, pmod, '-', linewidth = LW, color = color)
        ax[1,1].plot(t, pphi, '-', linewidth = LW, color = color)
        ax[0,2].plot(t, qmod, '-', linewidth = LW, color = color)
        ax[1,2].plot(t, qphi, '-', linewidth = LW, color = color)

        
    #add_legends(ax, loc)

    #plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()


#=================================================================

def setcanvas(ax):


    tmax = 1700
    
    iax = ax[0,0]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\Sigma n_k$')
    iax.grid()
    iax.set_ylim((0,0.03))
    iax.set_xlim((0,tmax))

    iax = ax[1,0]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\Sigma n_k k^4$')
    iax.grid()
    iax.set_ylim((0,25))
    iax.set_xlim((0,tmax))

    
    iax = ax[0,1]
    iax.set_xlabel('$t$')
    iax.set_ylabel('momentum')
    iax.grid()
    iax.set_ylim((0,0.02))
    iax.set_xlim((0,tmax))

    
    iax = ax[1,1]
    iax.set_xlabel('$t$')
    iax.set_ylabel('momentum angle / $\pi$')
    iax.grid()
    iax.set_ylim((-1,1))
    iax.set_xlim((0,tmax))


    iax = ax[0,2]
    iax.set_xlabel('$t$')
    iax.set_ylabel('anisotropy')
    iax.grid()
    iax.set_ylim((0,2))
    iax.set_xlim((0,tmax))

    iax = ax[1,2]
    iax.set_xlabel('$t$')
    iax.set_ylabel('anisotropy angle / $\pi$')
    iax.grid()
    iax.set_ylim((-1,1))
    iax.set_xlim((0,tmax))

#    ax[1,0].set_ylim(0,1.5)
#    ax[1,1].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 

            
#=================================================================


def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

#=================================================================

