#import importlib ; import plot_flatIring
#importlib.reload(plot_flatIring); plot_flatIring.Plot();



#---------------------------------------------------------------------


def Plot():
 
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatIring.pdf'
    
    run = ( (  "DATA_fold/flat_dx1n04_dz1nz100",     '',             ), #0
            (  "DATA_rebin11/flat_dx1n08_dz1nz200b", 'inf, 8k', "or" ), #1
            (  "DATA_rebin12/flat_dx1n12_dz1nz200b", '',              ), #2
            (  "DATA_rebin13/flat_dx1n16_dz1nz200b", 'inf, 16k', "-r" ), #3
            (  "DATA_rebin14/fr08_dx1n04_dz1nz200",  ''               ), #4
            (  "DATA_rebin15/fr08_dx1n08_dz1nz200",  ''               ), #5
            (  "DATA_rebin16/fr08_dx1n16_dz1nz200",  '',              ), #6
            (  "DATA_rebin17/fr13_dx1n04_dz1nz200",  '',              ), #7
            (  "DATA_rebin18/fr26_dx1n04_dz1nz200",  '',              ), #8
            (  "DATA_rebin19/fr51_dx1n04_dz1nz200",  '',              ), #9
            (  "DATA_rebin20/fr13_dx1n08_dz1nz200",  '', "-b" ),  #10
            (  "DATA_rebin21/fr51_dx1n08_dz1nz200",  '', "-g" ),  #11
            (  "DATA_rebin22/fr80_dx1n04_dz1nz200",  '', "ob" ),  #12
            (  "DATA_rebin23/f120_dx1n04_dz1nz200",  '', "og" ),  #13
            (  "DATA_rebin24/f160_dx1n04_dz1nz200",  '', "oy" ),  #14
            (  "DATA_rebin25/f200_dx1n04_dz1nz200",  '', "om" ),  #15
            (  "DATA_rebin26/fr80_dx1n08_dz1nz200",  '', "-b" ),  #16
            (  "DATA_rebin27/f120_dx1n08_dz1nz200",  '', "-g" ),  #17
            (  "DATA_rebin28/f160_dx1n08_dz1nz200",  '', "-y" ),  #18
            (  "DATA_rebin29/f200_dx1n08_dz1nz200",  '200 cm, 8k', "og" ),  #19
            (  "DATA_rebin30/f600_dx1n08_dz1nz200",  '600 cm, 8k', "ob" ),  #20
            (  "DATA_rebin31/f200_dx1n16_dz1nz200",  '200 cm, 16k', "-g" ),  #21
            (  "DATA_rebin32/f200_dx1n12_dz1nz200",  '200 cm, 12k', "og" )   #22
    )

    MFC = 'none'                        # color - marker face color
    MS  = 2

    dz=1.062944882793539

    iz= ((5, 10, 20, 40),  (50, 60, 80, 100),  (120, 140, 160, 180))

    xIhi = ((60, 60, 40, 35),  ( 30, 30, 25, 25),  (25, 25, 25, 25)) 
#   xIhi = ((5, 5, 5, 5),  ( 5, 5, 5, 5),  (8, 8, 8, 8)) 
    yPlo = ((1e-10, 1e-10, 1e-10, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10)) 
     

    #-- set canvas --

    fig, ax = plt.subplots(ncols=4, nrows=3, tight_layout=True, figsize=(12,9))


    for i in range (0,4):        
        for j in range (0,3):

            ax[j,i].set(title='z = {:.2f} km'.format(iz[j][i]*dz))               
 
            ax[j,i].set(xlabel = "$I$")
            ax[j,i].set(ylabel = "PDF")
            ax[j,i].set_yscale('log')
            ax[j,i].set_xlim(0, xIhi[j][i])
            ax[j,i].set_ylim(yPlo[j][i], 10)
            ax[j,i].grid()
  
 
    #-- add data --

    for f in (1, 3, 20, 19, 21): 
       for i in (0, 1, 2, 3):
            for j in range (0,3):

                try:
                    fname = run[f][0] + "_iz" + str(iz[j][i]).zfill(4) + "_probI.dat"
                    #print(fname)
                    dat = np.loadtxt(fname)
                    ax[j,i].plot(dat[:,0], dat[:,1], run[f][2], markersize=MS,  mfc=MFC, label=run[f][1])
                except:          
                    pass

            
    #-- add model curves --

    x, y = modelLinear(1.0, 0, 25)

    for i in range (0,4):
        for j in range (0,3):
            ax[j,i].plot(x, y, '--k', linewidth=1) # label = 'c = 1.0')

    
    
    #for i in range (0,4):

    #    x,y,lbl = model(i)
    #    ax[0,i].plot(x, y, '-k', linewidth=1, label=lbl)
    
    #-- show and export --

    for i in range (0,4):

        ax[0,i].legend()
        ax[1,i].legend()
        ax[2,i].legend()
    
    #plt.tight_layout()
    plt.savefig(fnameout, pad_inches=0)

    plt.close()



#---------------------------------------------------------------------


    
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
