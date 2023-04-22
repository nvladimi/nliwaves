#import importlib ; import plot_klo
#importlib.reload(plot_klo) ; plot_klo.Plot()


#=================================================================

def Plot(ampname="a2e3", fnameout='test.pdf'):
    import matplotlib.pyplot as plt
    import numpy as np

    kmax = (2, 3, 4, 6, 9, 10)
       
    fig, ax = plt.subplots(ncols=5, nrows=6, figsize=(10,12), constrained_layout = True)

    setcanvas(ax)

    for i in range(6):
        plotKloRow(ax, i, ampname, kmax[i])
                 
    #add_legends(ax, loc)

    #plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    #plt.show()
    plt.close()

    
#=================================================================

def plotKloRow(ax, row, amp, kmax):
    import matplotlib.pyplot as plt
    import numpy as np
    import resonant
    
    fname = "runsA/dat/klo_m" + str(2*kmax+1).zfill(2) + amp + ".npy" 

    klo = np.load(fname)
    klo = np.fft.fftshift(klo)
    klo = klo[1:,1:]
    logklo = np.log(klo)
    #print(np.min(logklo), np.max(logklo))
    
    k = np.arange(-kmax,kmax+1) 
    kx, ky = np.meshgrid(k,k)
    k4 = (kx*kx + ky*ky)**2
    k4klo = k4*klo
    #print(np.min(k4klo), np.max(k4klo))
    
    side, diag = resonant.Modes(kmax)
    cmax = np.max(side+diag)+1
    
    ax[row, 0].imshow(logklo, cmap='jet', origin = 'lower', vmin=-20, vmax=-4);
    ax[row, 1].imshow(k4klo,  cmap='jet', origin = 'lower', vmin=0,   vmax=0.005);

    ax[row, 2].imshow(side+diag, cmap='jet', origin = 'lower', vmin=0, vmax=cmax);
    ax[row, 3].imshow(diag,      cmap='jet', origin = 'lower', vmin=0, vmax=cmax);
    ax[row, 4].imshow(side,      cmap='jet', origin = 'lower', vmin=0, vmax=cmax);

#=================================================================

def setcanvas(ax):
    import matplotlib.pyplot as plt

    
    for i in range(ax.shape[0]):
        for j in range(ax.shape[1]):
            ax[i,j].set_xticks(())
            ax[i,j].set_yticks(())

"""                    
            ax[i,j].set_title(modes[i][j])
            ax[i,j].set_xlabel('$t$')
            ax[i,j].set_ylabel(ylabel)
            ax[i,j].set_ylim(y1,y2)
            ax[i,j].grid()
            ax[i,j].set_xlim(0,1e6)
            ax[i,j].ticklabel_format(axis='x', style='sci', scilimits =(-3,0)) 
            if i == 2 and varname == "NL":
                ax[i,j].set_ylim(-0.12,0.02)
"""               
            
            
#=================================================================

#=================================================================

    
#=================================================================

def add_legends(ax, loc):
    import matplotlib.pyplot as plt
    
    for i in (0,1,2,3):
        for j in (0,1):
 
            ax[i,j].legend(frameon=False, loc=loc)

   
#=================================================================


#=================================================================

#=================================================================

