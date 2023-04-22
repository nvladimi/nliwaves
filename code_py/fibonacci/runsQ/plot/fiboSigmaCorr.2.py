
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

    
    with open('../PX/q7_v0dt1_s4_i40_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Ca3 = C/C[round(n/2),0]


    #--------
    
    with open('../PX/q7b_g3500dt5_s4_i20_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Cb1 = C/C[round(n/2),0]


    with open('../PX/q7b_g3500dt5_s4_i32_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Cb2 = C/C[round(n/2),0]

    
    with open('../PX/q7b_g3500dt5_s4_i48_sigmacorr.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    C   = dat["corrfun"]
    n   = C.shape[0]
    Cb3 = C/C[round(n/2),0]

    
    tA1 = (np.arange(0,n) - round(n/2))*0.1
    tB1 = (np.arange(0,n) - round(n/2))*1e-5



    #-- plot --

    fig, ax = plt.subplots(ncols=2, nrows=3, tight_layout=True, figsize=(8,8))

    MS = 2
    LW = 1

    def f(j):
        return(-2)**j
    
    
    iax = ax[0,0]
 
    iax.set_xlim(-1,1)
    #iax.set_ylim(1e-6,1)
    iax.set(xlabel='$\\tau$', ylabel='$(-2)^j C_j$')
    iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
 
    iax.plot(tA1, Ca1[:,0]*f(0), 'ro-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(tA1, Ca1[:,1]*f(1), 'bo-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(tA1, Ca1[:,2]*f(2), 'go-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(tA1, Ca1[:,3]*f(3), 'yo-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(tA1, Ca1[:,4]*f(4), 'mo-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");
    
    iax.grid(); iax.legend(frameon=True)



    iax = ax[1,0]
 
    iax.set_xlim(-4,4)
    #iax.set_ylim(1e-6,1)
    iax.set(xlabel='$\\tau$', ylabel='$(-2)^j C_j$')
    iax.text(0.1, 0.9, '$\\alpha=0, \, i = 24$', transform=iax.transAxes)

    
    iax.plot(tA1, Ca2[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(tA1, Ca2[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(tA1, Ca2[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(tA1, Ca2[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(tA1, Ca2[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");
    iax.grid(); iax.legend(frameon=True)


    iax = ax[2,0]
 
    iax.set_xlim(-30,30)
    #iax.set_ylim(1e-6,1)
    iax.set(xlabel='$\\tau$', ylabel='$(-2)^j C_j$')
    iax.text(0.1, 0.9, '$\\alpha=0, \, i = 40$', transform=iax.transAxes)
    
    iax.plot(tA1, Ca3[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(tA1, Ca3[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(tA1, Ca3[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(tA1, Ca3[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(tA1, Ca3[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");
    iax.grid(); iax.legend(frameon=True)


    #------
    
    
    iax = ax[0,1]
 
    iax.set_xlim(-3e-3,3e-3)
    #iax.set_ylim(1e-6,1)
    iax.set(xlabel='$\\tau$', ylabel='$(-2)^j C_j$')
    iax.text(0.1, 0.9, '$\\alpha=1, \, i = 20$', transform=iax.transAxes)
 
    iax.plot(tB1, Cb1[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(tB1, Cb1[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(tB1, Cb1[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(tB1, Cb1[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(tB1, Cb1[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");
    iax.grid(); iax.legend(frameon=True)



    iax = ax[1,1]
 
    iax.set_xlim(-5e-4,5e-4)
    #iax.set_ylim(1e-6,1)
    iax.set(xlabel='$\\tau$', ylabel='$(-2)^j C_j$')
    iax.text(0.1, 0.9, '$\\alpha=1, \, i = 32$', transform=iax.transAxes)

    
    iax.plot(tB1, Cb2[:,0]*f(0), 'r-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(tB1, Cb2[:,1]*f(1), 'b-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(tB1, Cb2[:,2]*f(2), 'g-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(tB1, Cb2[:,3]*f(3), 'y-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(tB1, Cb2[:,4]*f(4), 'm-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");
    iax.grid(); iax.legend(frameon=True)


    iax = ax[2,1]
 
    iax.set_xlim(-1e-4,1e-4)
    #iax.set_ylim(1e-6,1)
    iax.set(xlabel='$\\tau$', ylabel='$(-2)^j C_j$')
    iax.text(0.1, 0.9, '$\\alpha=1, \, i = 48$', transform=iax.transAxes)
    
    iax.plot(tB1, Cb3[:,0]*f(0), 'ro-',  mfc='none', ms=MS,  lw=LW,  label="$j = 0$");
    iax.plot(tB1, Cb3[:,1]*f(1), 'bo-',  mfc='none', ms=MS,  lw=LW,  label="$j = 1$");
    iax.plot(tB1, Cb3[:,2]*f(2), 'go-',  mfc='none', ms=MS,  lw=LW,  label="$j = 2$");
    iax.plot(tB1, Cb3[:,3]*f(3), 'yo-',  mfc='none', ms=MS,  lw=LW,  label="$j = 3$");
    iax.plot(tB1, Cb3[:,4]*f(4), 'mo-',  mfc='none', ms=MS,  lw=LW,  label="$j = 4$");
    iax.grid(); iax.legend(frameon=True)




    #------

        
    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
 






    

#=====================================



#=========================================================
