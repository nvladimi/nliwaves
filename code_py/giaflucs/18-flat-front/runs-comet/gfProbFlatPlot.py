#import importlib ; import gfProbFlatPlot
#importlib.reload(gfProbFlatPlot) ; gfProbFlatPlot.Plot()


def Plot():
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatIpdf.pdf'
    
    fbase = ( "dx1n04_dz1nz100", "dx1n08_dz1nz100", "dx2n08_dz1nz100", "dx1n04_dz2nz200" )

    ML  =   ("base", "2Lx2L", "dx/2", "dz/2")
    MT  =   ("or", "sg", "vb", "^c")
    MFC = 'none'                        # color - marker face color
    MS  = 2

    iz=(1, 5, 20, 100)
    dz=1.062944882793539
   
    xIhi = {1: 7.0, 5: 60, 20: 50.0, 100:30.0} 
     

    #-- set canvas --

    fig, ax = plt.subplots(ncols=4, nrows=2, constrained_layout=True, figsize=(12,6))

    for i in range (0,4):

        ax[0,i].set(title='z = {:.2f} km'.format(iz[i]*dz))               

        ax[0,i].set(xlabel = "$I$")
        ax[0,i].set(ylabel = "PDF")
        ax[0,i].set_yscale('log')
        ax[0,i].set_xlim(0, xIhi[iz[i]])
        ax[0,i].set_ylim(1e-10, 10)
        ax[0,i].grid()

        ax[1,i].set(xlabel = "$(\ln I  - \langle \ln I \\rangle ) /\sigma$")
        ax[1,i].set(ylabel = "PDF")
        ax[1,i].set_yscale('log')
        ax[1,i].set_xlim(-5, 4)
        ax[1,i].set_ylim(1e-8, 10)
        ax[1,i].grid()

    #-- add data --

    for f in (0, 1, 2, 3):
        for i in range (0,4):

            iiz = iz[i]
            if f==3: iiz = 2*iiz
            
            fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
            dat = np.loadtxt(fname)
            ax[0,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])

            fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_problogI.dat"    
            dat = np.loadtxt(fname)
            ax[1,i].plot(dat[:,0], dat[:,1], MT[f],  markersize=MS,  mfc=MFC, label=ML[f])

            
    #-- add model curves --

    for i in range (0,4):

        x,y,lbl = model(i)
        ax[0,i].plot(x, y, '-k', linewidth=1, label=lbl)
    
    #-- show and export --

    for i in range (0,4):

        ax[0,i].legend()
    
    plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    #plt.show()


#---------------------------------------------------------------------

def model(i):
    import numpy as np
 
    a = (3.3, 1.8, 1.8, 0.95)
    q = (1, 0.62,  0.65, 1)
   
    lbl = ('$\exp(-3.3 \, I)$',
           '$\exp(-1.8 \, I^{0.62})$',
           '$\exp(-1.8 \, I^{0.65})$',
           '$\exp(-0.95 \, I)$')
    
    dx = 0.1
    x  = np.arange(0, 100, dx)
    y  = np.exp( -a[i]*x**q[i])
    y  = y / sum(y)/dx
    
    return x, y, lbl[i]

    
#---------------------------------------------------------------------


    

#=================================================================================================
