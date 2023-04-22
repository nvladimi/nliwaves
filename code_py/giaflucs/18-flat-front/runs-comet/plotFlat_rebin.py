#import importlib ; import plotFlat_rebin
#importlib.reload(plotFlat_rebin) ; plotFlat_rebin.Plot()


#--------------------------------------------------------------------


def Plot(c1=9, c2=11):
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatRebin.pdf'
    
    fbase = ( "DATA/flat_dx1n04_dz1nz100", 
              "DATA/flat_dx1n12_dz1nz100",
              "DATA_trim/flat_dx1n12_dz1nz100",
              "DATA_rebin01/flat_dx1n12_dz1nz100",            
              "DATA_rebin02/flat_dx1n12_dz1nz100"              
    )

    ML  =   ("base", '3Lx3L', 'trim', 'rebin01', 'rebin02')
    MT  =   ("or", "or", "oy", "ob", 'ob')
    MFC = 'none'                        # color - marker face color
    MS  = 0.5

    dz=1.062944882793539

    iz= (( 1,  5, 10, 20), (30, 40, 50, 60), (70, 80, 90, 100))

    xIhi = ((7, 60, 60, 40), (40, 35, 30, 30), (25, 25, 25, 25)) 
    yPlo = ((1e-10, 1e-10, 1e-10, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10)) 
     

    #-- set canvas --

    fig, ax = plt.subplots(ncols=4, nrows=3, constrained_layout=True, figsize=(12,9))


    for i in range (0,4):        
        for j in range (0,3):

            ax[j,i].set(title='z = {:.2f} km'.format(iz[j][i]*dz))               
 
            ax[j,i].set(xlabel = "$I$")
            ax[j,i].set(ylabel = "PDF")
            ax[j,i].set_yscale('log')
            ax[j,i].set_xlim(0, xIhi[j][i])
            ax[j,i].set_ylim(yPlo[j][i], 10)
            #ax[j,i].grid()
  
 
    #-- add data --

    for f in (2, 4):
        for i in (0, 1, 2, 3):
            for j in range (0,3):
           
                fname = fbase[f] + "_iz" + str(iz[j][i]).zfill(4) + "_probI.dat"   
                dat = np.loadtxt(fname)
                ax[j,i].plot(dat[:,0], dat[:,1], MT[f], label=ML[f],
                             markersize=MS,  mfc=MFC, linewidth=1)



            
    #-- add model curves --

    x, y = modelLinear(1.0, 0, 25)

    for i in range (0,4):
        for j in range (0,3):
            ax[j,i].plot(x, y, '--k', linewidth=0.5, label = 'exp(-I)')

    
    
    #for i in range (0,4):

    #    x,y,lbl = model(i)
    #    ax[0,i].plot(x, y, '-k', linewidth=1, label=lbl)
    
    #-- show and export --

    for i in range (0,4):

        ax[0,i].legend()
        ax[1,i].legend()
        ax[2,i].legend()
    
    plt.tight_layout()
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
    
    dx = (x2-x1)/100
    x  = np.arange(x1, x2, dx)
    y  = np.exp( - x/a )
    y  = y / sum(y)/dx
    
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
