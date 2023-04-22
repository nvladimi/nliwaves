
#import importlib ; import fiboSigmaCorr
#importlib.reload(fiboSigmaCorr) ; fiboSigmaCorr.all()


def all():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    import pickle
    

    outfile = "fiboSigmaCorr.pdf"


    
    #-- read data  --


    with open('../PX/q7_v0dt1_s4_i12_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Ca1 = C/C[round(n/2),0]

    
    with open('../PX/q7_v0dt1_s4_i24_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Ca2 = C/C[round(n/2),0]

    
    with open('../PX/q7_v0dt1_s4_i36_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Ca3 = C/C[round(n/2),0]


    #--------
    
    with open('../PX/q7b_g3500dt5_s4_i12_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Cb1 = C/C[round(n/2),0]


    with open('../PX/q7b_g3500dt5_s4_i24_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Cb2 = C/C[round(n/2),0]

    
    with open('../PX/q7b_g3500dt5_s4_i36_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Cb3 = C/C[round(n/2),0]

    
    #-------
    
    tA1 = (np.arange(0,n) - round(n/2))*0.1
    tB1 = (np.arange(0,n) - round(n/2))*1e-5


    Fi = Fibonacci(60)

    Fa = np.power(Fi, -1/3.)
    Fb = np.power(Fi,  1/3.)



    Ca1s = (Ca1[1:,:] +  np.flipud(Ca1)[:-1,:])/2
    Ca2s = (Ca2[1:,:] +  np.flipud(Ca2)[:-1,:])/2
    Ca3s = (Ca3[1:,:] +  np.flipud(Ca3)[:-1,:])/2
    
    Cb1s = (Cb1[1:,:] +  np.flipud(Cb1)[:-1,:])/2
    Cb2s = (Cb2[1:,:] +  np.flipud(Cb2)[:-1,:])/2
    Cb3s = (Cb3[1:,:] +  np.flipud(Cb3)[:-1,:])/2


    Ca1u = (Ca1[1:,:] -  np.flipud(Ca1)[:-1,:])/2
    Ca2u = (Ca2[1:,:] -  np.flipud(Ca2)[:-1,:])/2
    Ca3u = (Ca3[1:,:] -  np.flipud(Ca3)[:-1,:])/2
    
    Cb1u = (Cb1[1:,:] -  np.flipud(Cb1)[:-1,:])/2
    Cb2u = (Cb2[1:,:] -  np.flipud(Cb2)[:-1,:])/2
    Cb3u = (Cb3[1:,:] -  np.flipud(Cb3)[:-1,:])/2



    

    
    #-- plot --


    
    fig, ax = plt.subplots(ncols=2, nrows=3, tight_layout=True, figsize=(10,10))

    MS = 2
    LW = 1

    y0 = 0.06
    
    def f(j):
        return ((-1)**j * np.exp(j))




    #----------------------------------------

    
    
    iax = ax[0,0]
 
    iax.set_xlim(-y0,y0)
    iax.set_ylim(-2,2)
    iax.set(xlabel='$\\tau  F_i^{-1/3}$', ylabel='$(-1)^j \exp(j)$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
    iax.axhline(lw=0.5, c="gray")
    iax.axvline(lw=0.5, c="gray")

    
    t = tA1 * Fa[12-1]
    
    iax.plot(t, Ca1[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);
    
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 24$', transform=iax.transAxes)

    
    t = tA1 * Fa[24-1]

    iax.plot(t, Ca2[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW);

    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 36$', transform=iax.transAxes)

    t = tA1 * Fa[36-1]

    iax.plot(t, Ca3[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(t, Ca3[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(t, Ca3[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(t, Ca3[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(t, Ca3[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");

    iax.legend(frameon=True)


    #------
    
    
    iax = ax[0,1]
 
    iax.set_xlim(-y0,y0)
    iax.set_ylim(-1,3)
    iax.set(xlabel='$\\tau F_i^{1/3}$', ylabel='$(-1)^j \exp(j)$')
    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 12$', transform=iax.transAxes)
    iax.axhline(lw=0.5, c="gray")
    iax.axvline(lw=0.5, c="gray")

    t = tB1 * Fb[12-1]
    
    iax.plot(t, Cb1[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    #iax.plot(t, Cb1[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);

    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 24$', transform=iax.transAxes)

    t = tB1 * Fb[24-1]
    
    iax.plot(t, Cb2[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(t, Cb2[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(t, Cb2[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(t, Cb2[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(t, Cb2[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");

    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 36$', transform=iax.transAxes)

    t = tB1 * Fb[36-1]
    
    iax.plot(t, Cb3[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    #iax.plot(t, Cb3[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);
    
    iax.legend(frameon=True)


    #-------------------------------------------
  
    
    iax = ax[1,0]
 
    iax.set_xlim(-y0,y0)
    iax.set_ylim(-1,2)
    iax.set(xlabel='$\\tau  F_i^{-1/3}$', ylabel='$(-1)^j \exp(j)$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
    iax.axhline(lw=0.5, c="gray")
    iax.axvline(lw=0.5, c="gray")


    t = tA1[1:] * Fa[12-1]
    
    iax.plot(t, Ca1s[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1s[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1s[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1s[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1s[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);
    

    
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 24$', transform=iax.transAxes)

    t = tA1[1:] * Fa[24-1]

    iax.plot(t, Ca2s[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(t, Ca2s[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(t, Ca2s[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(t, Ca2s[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(t, Ca2s[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");

    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 36$', transform=iax.transAxes)

    t = tA1[1:] * Fa[36-1]

    iax.plot(t, Ca3s[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca3s[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca3s[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca3s[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca3s[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);

    
    #------
    
    
    iax = ax[1,1]
 
    iax.set_xlim(-y0,y0)
    iax.set_ylim(-1,2)
    iax.set(xlabel='$\\tau F_i^{1/3}$', ylabel='$(-1)^j \exp(j)$')
    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 12$', transform=iax.transAxes)
    iax.axhline(lw=0.5, c="gray")
    iax.axvline(lw=0.5, c="gray")

    t = tB1[1:] * Fb[12-1]
    
    iax.plot(t, Cb1s[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1s[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1s[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1s[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    #iax.plot(t, Cb1s[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);

    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 24$', transform=iax.transAxes)

    t = tB1[1:] * Fb[24-1]
    
    iax.plot(t, Cb2s[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(t, Cb2s[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(t, Cb2s[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(t, Cb2s[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(t, Cb2s[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");

    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 36$', transform=iax.transAxes)

    t = tB1[1:] * Fb[36-1]
    
    iax.plot(t, Cb3s[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3s[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3s[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3s[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    #iax.plot(t, Cb3s[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);
    


  #-------------------------------------------
  
    
    iax = ax[2,0]
 
    iax.set_xlim(-y0,y0)
    iax.set_ylim(-1.5,1.5)
    iax.set(xlabel='$\\tau  F_i^{-1/3}$', ylabel='$(-1)^j \exp(j)$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
    iax.axhline(lw=0.5, c="gray")
    iax.axvline(lw=0.5, c="gray")


    t = tA1[1:] * Fa[12-1]
    
    iax.plot(t, Ca1u[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1u[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1u[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1u[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca1u[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);
    

    
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 24$', transform=iax.transAxes)

    t = tA1[1:] * Fa[24-1]

    iax.plot(t, Ca2u[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2u[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2u[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2u[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Ca2u[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW);

    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 36$', transform=iax.transAxes)

    t = tA1[1:] * Fa[36-1]

    iax.plot(t, Ca3u[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(t, Ca3u[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(t, Ca3u[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(t, Ca3u[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(t, Ca3u[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");

 

    
    #------
    
    
    iax = ax[2,1]
 
    iax.set_xlim(-y0,y0)
    iax.set_ylim(-1.5,1.5)
    iax.set(xlabel='$\\tau F_i^{1/3}$', ylabel='$(-1)^j \exp(j)$')
    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 12$', transform=iax.transAxes)
    iax.axhline(lw=0.5, c="gray")
    iax.axvline(lw=0.5, c="gray")

    t = tB1[1:] * Fb[12-1]
    
    iax.plot(t, Cb1u[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1u[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1u[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb1u[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    #iax.plot(t, Cb1u[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);

    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 24$', transform=iax.transAxes)

    t = tB1[1:] * Fb[24-1]
    
    iax.plot(t, Cb2u[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(t, Cb2u[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(t, Cb2u[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(t, Cb2u[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(t, Cb2u[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");

    #iax.text(0.1, 0.9, '$\\alpha=1, \, i = 36$', transform=iax.transAxes)

    t = tB1[1:] * Fb[36-1]
    
    iax.plot(t, Cb3u[:,0]*f(0), 'r--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3u[:,1]*f(1), 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3u[:,2]*f(2), 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(t, Cb3u[:,3]*f(3), 'y--',  mfc='none', ms=MS,  lw=LW);
    #iax.plot(t, Cb3u[:,4]*f(4), 'm--',  mfc='none', ms=MS,  lw=LW);
    

    

    #------

        
    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
 






    

#=====================================

def Fibonacci(n):

    import numpy as np
     
    Fi = np.zeros(n, 'uint64')
     
    Fi[0]=1
    Fi[1]=1

    for i in range(2,n):
        
        Fi[i] = Fi[i-1] + Fi[i-2]
        
    return(Fi)
  

#=========================================================
