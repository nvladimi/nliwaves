#import importlib ; import figVsTimeLong
#importlib.reload(figVsTimeLong) ; figVsTimeLong.All()


import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np
import nlwTools


#=================================================================

def All(every=20):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import nlwTools

    fnameout = 'figVsTimeLong.pdf'
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
    
    fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(6,3),  constrained_layout=True)

    iax = ax[0]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$|P|$')
    iax.grid()
    iax.set_ylim(0,0.01)

    iax = ax[1]
    iax.set_xlabel('$t$')
    iax.set_ylabel('$\\theta  / \pi$')
    iax.grid()
    iax.set_ylim(0,1)
    ticks = np.arange(0,1.01,0.25)
    iax.set_yticks(ticks)
    iax.set_yticklabels(["$%g$"%y for y in ticks])
    
    s = "a1"

    fname = prefix + s + suffix

    dat = np.loadtxt(fname)  # (t, nk, nq, pkx, pky, qkx, qky)

    ind = np.arange(0, dat.shape[0], every)
    dat = dat[ind,:]

    t     = dat[:,0]
    pkx   = dat[:,3]
    pky   = dat[:,4]


    pmod = np.sqrt(pkx**2 + pky**2) / 2     # <-- to take into account sqrt(2) definition of a_k 
    pphi = np.arctan2(pky, pkx)/np.pi

    
    color = seedcolor[s]

    ax[0].plot(t, pmod, '-', linewidth = LW, color = color)
    ax[1].plot(t, pphi, '.', linewidth = LW/2, color = color, markersize=1)

        
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()

