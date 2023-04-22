
#import importlib ; import fiboPXprob_1434
#importlib.reload(fiboPXprob_1434) ; fiboPXprob_1434.all()


def all():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    import pickle
    

    outfile = "fiboPXprob_1434.pdf"
    
    #-- read data ---------------------

    with open('../POSTFrisch/q7_v0dt1_probX.pkl', 'rb')  as handle:
    #with open('../PX/q7_v0dt1_probX.pkl', 'rb') as handle:
        dat = pickle.load(handle)
        
    P0 = dat["Xhist"]
    x0 = dat["centers"]
    ntot = dat["ntot"]
    dx=x0[2]-x0[1]
   
    ncount = np.sum(P0,0)

    #print (ntot, ncount, dx)
    P0 =  P0/ntot/dx



    with open('../POSTFrisch/q7b_g3500dt5_probX.pkl', 'rb') as handle:
    #with open('../PX/q7b_g3500dt5_probX.pkl', 'rb') as handle:
        dat = pickle.load(handle)

    P1 = dat["Xhist"]
    x1 = dat["centers"]
    ntot = dat["ntot"]
    dx=x1[2]-x1[1]

    ncount = np.sum(P1,0)

    #print (ntot, ncount, dx)
    P1 =  P1/ntot/dx


    with open('../POSTFrisch/q7_v14L10_probX.pkl', 'rb') as handle:
        dat = pickle.load(handle)

    P14 = dat["Xhist"]
    x14 = dat["centers"]
    ntot = dat["ntot"]
    dx=x14[2]-x14[1]

    ncount = np.sum(P14,0)

    #print (ntot, ncount, dx)
    P14 =  P14/ntot/dx
    

    with open('../POSTFrisch/q7_v34R28_probX.pkl', 'rb') as handle:
        dat = pickle.load(handle)

    P34 = dat["Xhist"]
    x34 = dat["centers"]
    ntot = dat["ntot"]
    dx=x34[2]-x34[1]

    ncount = np.sum(P34,0)

    #print (ntot, ncount, dx)
    P34 =  P34/ntot/dx
    

    
    #-- plot ----------------------------

    fig, ax = plt.subplots(ncols=2, nrows=1, tight_layout=True, figsize=(7,3.5))

    MS = 2
    LW = 1

    k1=46; k2=40; k3=30; k4=20
    
       
    iax = ax[0]
 
    iax.set_xlim(-0.8,0.2)
    iax.set_ylim(-0.8,0)
 
    
    iax.set(xlabel='$2 X_k / k$', ylabel='$k^{-1} \ln(P(X_k) /2)$')
      
    #iax.plot(2*x0/k1, np.log(0.5*P0[:, 50-k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x0/k2, np.log(0.5*P0[:, 50-k2-1])/k2, 'r:',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x0/k3, np.log(0.5*P0[:, 50-k3-1])/k3, 'r--',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x0/k4, np.log(0.5*P0[:, 50-k4-1])/k4, 'r-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 0$");

    #iax.plot(2*x1/k1, np.log(0.5*P1[:, 10+k1-1])/k1, 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x1/k2, np.log(0.5*P1[:, 10+k2-1])/k2, 'b:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x1/k3, np.log(0.5*P1[:, 10+k3-1])/k3, 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x1/k4, np.log(0.5*P1[:, 10+k4-1])/k4, 'b-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 1$");

    #iax.plot(2*x14/k1, np.log(0.5*P14[:, 50-k1-1])/k1, 'y-',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x14/k2, np.log(0.5*P14[:, 50-k2-1])/k2, 'y:',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x14/k3, np.log(0.5*P14[:, 50-k3-1])/k3, 'y--',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x14/k4, np.log(0.5*P14[:, 50-k4-1])/k4, 'y-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 1/4$" );

    #iax.plot(2*x34/k1, np.log(0.5*P34[:, 10+k1-1])/k1, 'g-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x34/k2, np.log(0.5*P34[:, 10+k2-1])/k2, 'g:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x34/k3, np.log(0.5*P34[:, 10+k3-1])/k3, 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x34/k4, np.log(0.5*P34[:, 10+k4-1])/k4, 'g-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 3/4$");

    
    iax.plot(x1, x1, 'k--',  mfc='none', ms=MS,  lw=LW/2,  label="$2X/k$");

    
    iax.grid(True)
    iax.legend(frameon=True)


    iax = ax[1]

    
    iax.set_xlim(0,0.2)
    iax.set_ylim(-0.7,0)
    #iax.set_yscale('log')
    #iax.set_xscale('log')


       
    iax.set(xlabel='$2 X_k / k$', ylabel='$k^{-1} \ln(P(X_k) /2)$')
      
    #iax.plot(2*x0/k1, np.log(0.5*P0[:, 50-k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x0/k2, np.log(0.5*P0[:, 50-k2-1])/k2, 'r:',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x0/k3, np.log(0.5*P0[:, 50-k3-1])/k3, 'r--',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x0/k4, np.log(0.5*P0[:, 50-k4-1])/k4, 'r-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 0$");

    #iax.plot(2*x1/k1, np.log(0.5*P1[:, 10+k1-1])/k1, 'b-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x1/k2, np.log(0.5*P1[:, 10+k2-1])/k2, 'b:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x1/k3, np.log(0.5*P1[:, 10+k3-1])/k3, 'b--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x1/k4, np.log(0.5*P1[:, 10+k4-1])/k4, 'b-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 1$");

    #iax.plot(2*x14/k1, np.log(0.5*P14[:, 50-k1-1])/k1, 'y-',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x14/k2, np.log(0.5*P14[:, 50-k2-1])/k2, 'y:',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x14/k3, np.log(0.5*P14[:, 50-k3-1])/k3, 'y--',  mfc='none', ms=MS,  lw=LW );
    iax.plot(2*x14/k4, np.log(0.5*P14[:, 50-k4-1])/k4, 'y-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 1/4$" );

    #iax.plot(2*x34/k1, np.log(0.5*P34[:, 10+k1-1])/k1, 'g-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x34/k2, np.log(0.5*P34[:, 10+k2-1])/k2, 'g:',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x34/k3, np.log(0.5*P34[:, 10+k3-1])/k3, 'g--',  mfc='none', ms=MS,  lw=LW);
    iax.plot(2*x34/k4, np.log(0.5*P34[:, 10+k4-1])/k4, 'g-',  mfc='none', ms=MS,  lw=LW, label= "$\\alpha = 3/4$");


    
 

    
    iax.grid(True)
    iax.legend(frameon=True)


    #------------------
    
    #plt.show()

    
    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
 


#=========================================================
