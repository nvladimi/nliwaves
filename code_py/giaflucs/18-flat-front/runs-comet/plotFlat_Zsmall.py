#import importlib ; import plotFlat_Zsmall
#importlib.reload(plotFlat_Zsmall) ; plotFlat_Zsmall.Plot()


#--------------------------------------------------------------------


def Plot(c1=9, c2=11):
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatI_smallZ.pdf'
    
    fbase = ( "DATA/flat_dx1n04_dz8nz040", 
              "DATA/flat_dx1n04_dzAnz040",
              "DATA_trim/flat_dx1n04_dz1nz100a",
              "DATA_trim/flat_dx1n04_dzCnz400a", 
              "DATA/flat_dx1n04_dz1nz100" 
    )

    ML  =   ('dz/8', 'dz/16', '3Lx3L', 'dz/64', 'base')
    MT  =   ('sb', 'xb', '-r', '-b', 'xr')
    MFC = 'none'                        # color - marker face color
    MS  = 4

    dz=1.062944882793539/64

    iz= (( 4,  8, 12, 16), (32, 64, 96, 128), (192, 256, 320, 384))

    xIhi = ((2, 2, 2, 2), (2, 2, 2, 2), (2, 2, 2, 2)) 
    yPlo = ((1e-10, 1e-10, 1e-10, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10)) 
    yPhi = ((15, 9, 6, 5), (2.5, 1.4, 1.2, 1), (1, 1, 1.2, 1.4)) 
     

    #-- set canvas --

    fig, ax = plt.subplots(ncols=4, nrows=3, constrained_layout=True, figsize=(12,9))


    for i in range (0,4):        
        for j in range (0,3):

            ax[j,i].set(title='z = {:.2f} km'.format(iz[j][i]*dz))               
 
            ax[j,i].set(xlabel = "$I$")
            ax[j,i].set(ylabel = "PDF")
            #ax[j,i].set_yscale('log')
            ax[j,i].set_xlim(1e-4, 2)
            ax[j,i].set_ylim(1e-2, yPhi[j][i])
            #ax[j,i].set_ylim(1e-2, 120)
            #ax[j,i].grid()
  
 
    #-- add data --

    for f in (3,2,1,0,4):
        for i in (0, 1, 2, 3):
            for j in range (0,3):

                try:
                    iiz = iz[j][i]
                    
                    if f == 0:
                        if iiz%8 == 0:
                            iiz = int(iiz/8)
                        else:
                            iiz = -1
                    elif f == 1:
                        if iiz%4 == 0:
                            iiz = int(iiz/4)
                        else:
                            iiz = -1
                    elif f == 2 or f == 4:
                        if iiz%64 == 0:
                            iiz = int(iiz/64)
                        else:
                            iiz = -1
                        
                    fname = fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"                   
                    dat = np.loadtxt(fname)
                    #print(fname)
                    ax[j,i].plot(dat[:,0], dat[:,1], MT[f], label=ML[f],
                                 markersize=MS,  mfc=MFC, linewidth=0.3)
                except:
                    pass


            
    #-- add model curves --

    # x, y = modelLinear(1.0, 0, 2)

    # for i in range (0,4):
    #    for j in range (0,3):
    #        ax[j,i].plot(x, y, '--k', linewidth=0.5, label = 'exp(-I)')

    
    
    #for i in range (0,4):

    #    x,y,lbl = model(i)
    #    ax[0,i].plot(x, y, '-k', linewidth=1, label=lbl)
    
    #-- show and export --

    for i in range (0,4):

        ax[0,i].legend()
        ax[1,i].legend()
        ax[2,i].legend()

    #ax[1,1].legend()

    
    #plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    #plt.show()



#---------------------------------------------------------------------

def modelCubic(a, x1, x2, b=1e6):
    import numpy as np
    
    dx = (x2-x1)/100
    x  = np.arange(x1, x2, dx)
    y  = np.exp( -((x-1)/a)**2 + ((x-1)/b)**3 )
    #y  = np.exp( -((x-1)/a)**2 )
    y  = y / sum(y)/dx
    
    return x,y


#---------------------------------------------------------------------

def modelLinear(a, x1, x2):
    import numpy as np
    
    dx = (x2-x1)/1000
    x  = np.arange(x1, x2, dx)
    y  = np.exp( - x/a )
    #y  = y / sum(y)/dx
    
    return x,y



#---------------------------------------------------------------------

def modelTvelve(z, x1, x2, c1, c2):
    import numpy as np

    c0 = 1
    #c1 = 10
    #c2 = 30
    
    C1 = (z/c1)**(11/5.)
    C2 = (z/c2)**(-11/5.)

    
    dx = (x2-x1)/100
    x  = np.arange(x1, x2, dx)

    q = x / ( 1 + x**(1/2.) / (C1  + C2 * x**(1/12.)) )

    y  = np.exp( - c0 * q )
    y  = y / sum(y)/dx
    
    return x,y


#---------------------------------------------------------------------

def modelQ(i):
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
