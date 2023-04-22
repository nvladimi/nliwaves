#import importlib ; import plotCap256anisotropy
#importlib.reload(plotCap256anisotropy) ; plotCap256anisotropy.All()


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

    fnameout = 'a1n256_anisotropy.pdf'
    prefix = "../2020-01-06_klo/postdata/"
    suffix = "n256_evol.txt"
    
    seeds = ('a1', 'a2', 'a5')
    seedcolor = {'a1': 'green',
                 'a2': 'blue',
                 'a5': 'red'}
    label     = {'a1': 'forcing = 0.001',
                 'a2': 'forcing = 0.002',
                 'a5': 'forcing = 0.005'}

    
    LW = 1
    
    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(12,8),  constrained_layout=True)
    fig.suptitle('Old data - anisotropy of capillary runs on 256x256 grid', fontsize=16)

    setcanvas(ax)

    
    for s in seeds:

        fname = prefix + s + suffix
        
        dat = np.loadtxt(fname)  # (t, nk, nq, pkx, pky, qkx, qky)

        ind = np.arange(0, dat.shape[0], every)
        dat = dat[ind,:]
        
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
        
        ax[0,0].plot(t, nk,   '-', linewidth = LW, color = color, label=label[s])
        ax[1,0].plot(t, nq,   '-', linewidth = LW, color = color)
        ax[0,1].plot(t, pmod, '-', linewidth = LW, color = color)
        ax[1,1].plot(t, pphi, '.', linewidth = LW/2, color = color, markersize=1)
        ax[0,2].plot(t, qmod, '-', linewidth = LW, color = color)
        ax[1,2].plot(t, qphi, '.', linewidth = LW/2, color = color, markersize=1)

        
    #add_legends(ax, loc)

    ax[0,0].legend(frameon=False)

    
    #plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()


#=================================================================

def setcanvas(ax):


    iax = ax[0,0]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\Sigma n_k$')
    iax.grid()
    iax.set_ylim(0,0.14)

    iax = ax[1,0]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\Sigma n_k k^4$')
    iax.grid()
    iax.set_ylim(0,40)


    iax = ax[0,1]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$|P|$')
    iax.grid()
    iax.set_ylim(0,0.12)

    iax = ax[1,1]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\\theta_P / \pi$')
    iax.grid()
    iax.set_ylim(-1,1)

    iax = ax[0,2]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$|Q|$')
    iax.grid()
    iax.set_ylim(0,5)

    iax = ax[1,2]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\\theta_Q / \pi$')
    iax.grid()
    iax.set_ylim(-1,1)

    
#    ax[1,0].set_ylim(0,1.5)
#    ax[1,1].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 

            
#=================================================================


def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

#=================================================================

