
#import importlib ; import fiboPXprob
#importlib.reload(fiboPXprob) ; fiboPXprob.all()


def all():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    import pickle
    

    outfile = "fiboPXprob_orig.pdf"
    
    #-- read data --
    
    with open('../PX/q7_v0dt1_probX.pkl', 'rb') as handle:
        dat = pickle.load(handle)
        
    P0 = dat["Xhist"]
    x0 = dat["centers"]
    ntot = dat["ntot"]
    dx=x0[2]-x0[1]
   
    ncount = np.sum(P0,0)

    #print (ntot, ncount, dx)
    P0 =  P0/ntot/dx

       
    with open('../PX/q7b_g3500dt5_probX.pkl', 'rb') as handle:
        dat = pickle.load(handle)

    P1 = dat["Xhist"]
    x1 = dat["centers"]
    ntot = dat["ntot"]
    dx=x1[2]-x1[1]

    ncount = np.sum(P1,0)

    #print (ntot, ncount, dx)
    P1 =  P1/ntot/dx


    #-- plot --

    fig, ax = plt.subplots(ncols=2, nrows=1, tight_layout=True, figsize=(7,3.5))

    MS = 2
    LW = 1

    k1=46; k2=40; k3=30; k4=20
    
       
    iax = ax[0]
 
    iax.set_xlim(-0.4,0.1)
    iax.set_ylim(-0.7,0)
 
    
    iax.set(xlabel='$X_k / k$', ylabel='$k^{-1} \ln P(X_k)$')

        
    iax.plot(x0/k1, np.log(P0[:, 50-k1-1])/k1, 'b-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k1));
    iax.plot(x0/k2, np.log(P0[:, 50-k2-1])/k2, 'r-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k2));
    iax.plot(x0/k3, np.log(P0[:, 50-k3-1])/k3, 'g-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k3));
    iax.plot(x0/k4, np.log(P0[:, 50-k4-1])/k4, 'y-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k4));
    iax.plot(x1, 2*x1, 'k--',  mfc='none', ms=MS,  lw=LW/2,  label="$2X/k$");
    
    iax.grid(True)
    iax.legend(frameon=True)


    iax = ax[1]
 
    iax.set_xlim(-0.4,0.1)
    iax.set_ylim(-0.7,0)
 
    iax.set(xlabel='$X_k / k$', ylabel='$k^{-1} \ln P(X_k)$')
         

    iax.plot(x1/k1, np.log(P1[:, 10+k1-1])/k1, 'b-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k1));
    iax.plot(x1/k2, np.log(P1[:, 10+k2-1])/k2, 'r-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k2));
    iax.plot(x1/k3, np.log(P1[:, 10+k3-1])/k3, 'g-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k3));
    iax.plot(x1/k4, np.log(P1[:, 10+k4-1])/k4, 'y-',  mfc='none', ms=MS,  lw=LW, label= "$k = $" + " {}".format(k4));
    iax.plot(x1, 2*x1, 'k--',  mfc='none', ms=MS,  lw=LW/2,  label="$2X/k$");

    
    iax.grid(True)
    iax.legend(frameon=True)


    #------------------
    
    #plt.show()

    
    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
 


#=========================================================
