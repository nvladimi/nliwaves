#import importlib ; import plotCap512vstime
#importlib.reload(plotCap512vstime) ; plotCap512vstime.All()


#=================================================================

def All(every=20):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import nlwTools

    fnameout = 'a1n512_vstime.pdf'
    prefix = "../wexac/a1n512_20210123/a1"
    suffix = ".dat"
    
    seeds = ('s50', 's51', 's52', 's53')
    seedcolor = {'s50': 'C0',
                 's51': 'C1',
                 's52': 'C2',
                 's53': 'C3'}
    LW = 0.5
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(8,8),  constrained_layout=True)
    fig.suptitle('Progress of capillary runs on 512x512 grid', fontsize=16)

    setcanvas(ax)

    v = nlwTools.data_index()
   
    for s in seeds:

        fname = prefix + s + suffix
        
        dat = np.loadtxt(fname, comments="%")
        
        t   = dat[:,v["time"]]
        ind = nlwTools.CleanTime(t)       
        dat = dat[ind,:]
        ind = np.arange(0,len(dat), every)
        dat = dat[ind,:]
        
        t     = dat[:,v["time"]]
        
        hmin = dat[:,v["eta_min"]]
        hmax = dat[:,v["eta_max"]]
        hrms = dat[:,v["eta_rms"]]
        
        H    = dat[:,v["E_pot"]] + dat[:,v["E_kin"]] + dat[:,v["E_nl"]]
        NL   = dat[:,v["E_nl"]]/H

        color = seedcolor[s]
        
        ax[0,0].plot(t, hmin, '-', linewidth = LW, color = color)
        ax[0,0].plot(t, hmax, '-', linewidth = LW, color = color)
        ax[0,1].plot(t, hrms, '-', linewidth = LW, color = color)
        ax[1,0].plot(t, H,    '-', linewidth = LW, color = color)
        ax[1,1].plot(t, NL,   '-', linewidth = LW, color = color)

        
    #add_legends(ax, loc)

    #plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()


#=================================================================

def setcanvas(ax):
    import matplotlib.pyplot as plt

    ax[0,0].set_xlabel('$t$')
    ax[0,0].set_ylabel('$\eta_{min}$,  $\eta_{max}$')
    ax[0,0].grid()

    ax[0,1].set_xlabel('$t$')
    ax[0,1].set_ylabel('$\eta_{rms}$')
    ax[0,1].grid()
    
    ax[1,0].set_xlabel('$t$')
    ax[1,0].set_ylabel('$H$')
    ax[1,0].grid()
    
    ax[1,1].set_xlabel('$t$')
    ax[1,1].set_ylabel('$H_{NL} / H$')
    ax[1,1].grid()

#    ax[1,0].set_ylim(0,1.5)
#    ax[1,1].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 

            
#=================================================================


def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

#=================================================================

