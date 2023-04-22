#import importlib ; import gfProbFlatTheory
#importlib.reload(gfProbFlatTheory); gfProbFlatTheory.Plot();
#importlib.reload(gfProbFlatTheory); gfProbFlatTheory.PlotLargeZ()
#importlib.reload(gfProbFlatTheory); gfProbFlatTheory.PlotLargeL()


def Plot(c1=9, c2=11):
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatItheory.pdf'
    
    fbase = ( "dx1n04_dz1nz100", "dx1n08_dz1nz100",
              "dx2n08_dz1nz100", "dx1n04_dz2nz200",
              "dx1n04_dz8nz040", "dx1n04_dzAnz040",
              "dx1n12_dz1nz100")

    ML  =   ("base", "2Lx2L", "dx/2", "dz/2",  'dz/8', 'dz/16', '3Lx3L')
    MT  =   ("or", "sg", "vb", "^c", "om-", "oy-", "y-")
    MFC = 'none'                        # color - marker face color
    MS  = 2

    dz=1.062944882793539

    iz0=(1, 2, 3, 4)
    iz1=(8, 12, 24, 40)
    iz2=(10, 20, 60, 100)
#    iz2=(40, 60, 80, 100)
     
    xIhi = ((2, 2, 2, 2), (5, 8, 25, 60), (70, 50, 30, 30)) 
    yPlo = ((1e-8, 1e-8, 1e-8, 1e-8), (1e-8, 1e-8, 1e-8, 1e-10), (1e-10, 1e-10, 1e-10, 1e-10)) 
     

    #-- set canvas --

    fig, ax = plt.subplots(ncols=4, nrows=3, constrained_layout=True, figsize=(12,9))


    for i in range (0,4):

        ax[0,i].set(title='z = {:.2f} km'.format(iz0[i]*dz/16))               
        ax[1,i].set(title='z = {:.2f} km'.format(iz1[i]*dz/8))        
        ax[2,i].set(title='z = {:.2f} km'.format(iz2[i]*dz)) 
        
        for j in range (0,3):
    
            ax[j,i].set(xlabel = "$I$")
            ax[j,i].set(ylabel = "PDF")
            ax[j,i].set_yscale('log')
            ax[j,i].set_xlim(0, xIhi[j][i])
            ax[j,i].set_ylim(yPlo[j][i], 10)
            ax[j,i].grid()
  
 
    #-- add data --

    for i in range (0,4):

        f=5;  iiz = iz0[i]
            
        fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
        dat = np.loadtxt(fname)
        ax[0,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])


        f=4; iiz = iz1[i]
            
        fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
        dat = np.loadtxt(fname)
        ax[1,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])


    f=4;

    i=1; iiz=1; 
    fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
    dat = np.loadtxt(fname)
    ax[0,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])

    i=3; iiz=2; 
    fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
    dat = np.loadtxt(fname)
    ax[0,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])


        

    for f in (0, 1, 2, 3, 6):

        i= 0
        iiz = 1
        if f==3: iiz = 2*iiz
        
        fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
        dat = np.loadtxt(fname)
        ax[1,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])

        
        i= 3
        iiz = 5
        if f==3: iiz = 2*iiz
        
        fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
        dat = np.loadtxt(fname)
        ax[1,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])


        for i in (0, 1, 2, 3):

            iiz = iz2[i]
            if f==3: iiz = 2*iiz
            
            fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iiz).zfill(4) + "_probI.dat"   
            dat = np.loadtxt(fname)
            ax[2,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])



            
    #-- add model curves --

    x, y = modelCubic(0.06, 0.7, 1.3)
    ax[0,0].plot(x, y, '--k', linewidth=1, label = 'a = 0.06')

    x, y = modelCubic(0.09, 0.5, 1.5, b=1.0)
    ax[0,1].plot(x, y, '--k', linewidth=1, label = 'a = 0.09')

    x, y = modelCubic(0.12, 0.5, 1.7, b=0.3)
    ax[0,2].plot(x, y, '--k', linewidth=1, label = 'a = 0.12')

    x, y = modelCubic(0.15, 0.4, 2.0, b=0.35)
    ax[0,3].plot(x, y, '--k', linewidth=1, label = 'a = 0.15')
    
    x, y = modelLinear(0.38, 0.75, 8)
    ax[1,1].plot(x, y, '--k', linewidth=1, label = 'c = 0.38')


    x, y = modelLinear(1.0, 0, 25)

    ax[1,2].plot(x, y, '--k', linewidth=1, label = 'c = 1.0')
    ax[1,3].plot(x, y, '--k', linewidth=1, label = 'c = 1.0')
    for i in range (0,4):
        ax[2,i].plot(x, y, '--k', linewidth=1, label = 'c = 1.0')

        
    x, y = modelTvelve(6.64, 0, 60, c1, c2)
    ax[1,3].plot(x, y, '-k', linewidth=1)

    x, y = modelTvelve(13.29, 0, 60, c1, c2)
    ax[2,0].plot(x, y, '-k', linewidth=1)

    x, y = modelTvelve(26.57, 0, 50, c1, c2)
    ax[2,1].plot(x, y, '-k', linewidth=1)

    
    
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

#--------------------------------------------------------------------


def PlotLargeZ(c1=9, c2=11):
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatI_largeZ.pdf'
    
    fbase = ( "dx1n04_dz1nz100", "dx1n08_dz1nz100",
              "dx2n08_dz1nz100", "dx1n04_dz2nz200",
              "dx1n04_dz8nz040", "dx1n04_dzAnz040",
              "dx1n12_dz1nz100")

    ML  =   ("base", "2Lx2L", "dx/2", "dz/2",  'dz/8', 'dz/16', '3Lx3L')
    MT  =   ("or", "sg", "vb", "^c", "om-", "oy-", "y-")
    MFC = 'none'                        # color - marker face color
    MS  = 2

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
            ax[j,i].grid()
  
 
    #-- add data --

    for f in (0, 1, 6):
        for i in (0, 1, 2, 3):
            for j in range (0,3):
           
                fname = 'DATA/flat_' + fbase[f] + "_iz" + str(iz[j][i]).zfill(4) + "_probI.dat"   
                dat = np.loadtxt(fname)
                ax[j,i].plot(dat[:,0], dat[:,1], MT[f], markersize=MS,  mfc=MFC, label=ML[f])



            
    #-- add model curves --

    x, y = modelLinear(1.0, 0, 25)

    for i in range (0,4):
        for j in range (0,3):
            ax[j,i].plot(x, y, '--k', linewidth=1, label = 'c = 1.0')

    
    
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


def PlotLargeL(c1=9, c2=11):
    
    import numpy as np
 
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt

    #-- definitions --

    fnameout = 'flatI_largeL.pdf'
    
    run = ( (  "DATA_fold/flat_dx1n04_dz1nz100",     "base",  "or" ), #0
            (  "DATA_fold/flat_dx1n08_dz1nz100",     "2Lx2L", "sy" ), #1
            (  "DATA_fold/flat_dx2n08_dz1nz100",     "dx/2",  "vb" ), #2
            (  "DATA_fold/flat_dx1n04_dz2nz200",     "dz/2",  "^c" ), #3
            (  "DATA_fold/flat_dx1n04_dz8nz040",     'dz/8',  "om" ), #4
            (  "DATA_fold/flat_dx1n04_dzAnz040",     'dz/16', "oy" ), #5
            (  "DATA_fold/flat_dx1n12_dz1nz100",     '3Lx3L', "oy" ), #6
            (  "DATA_rebin11/flat_dx1n08_dz1nz200b", '2Lx2L', "-y" ), #7
            (  "DATA_rebin12/flat_dx1n12_dz1nz200b", '3Lx3L', "sm" ), #8
            (  "DATA_rebin13/flat_dx1n16_dz1nz200b", '4Lx4L', "vg" ), #9
            (  "DATA_trim11/flat_dx1n08_dz1nz200b",  '2Lx2L', ".y" ), #10
            (  "DATA_trim12/flat_dx1n12_dz1nz200b",  '3Lx3L', ".y" ), #11
            (  "DATA_trim13/flat_dx1n16_dz1nz200b", ' 4Lx4L', ".y" )  #12
    )

    MFC = 'none'                        # color - marker face color
    MS  = 2

    dz=1.062944882793539

    iz= ((5, 10, 20, 40),  (50, 60, 80, 100),  (120, 140, 160, 180))

    xIhi = ((60, 60, 40, 35),  ( 30, 30, 25, 25),  (25, 25, 25, 25)) 
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

    for f in (0, 7, 8, 9):
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
            ax[j,i].plot(x, y, '--k', linewidth=1, label = 'c = 1.0')

    
    
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
