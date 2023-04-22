
#import importlib ; import fiboPXprob_v2
#importlib.reload(fiboPXprob_v2) ; fiboPXprob_v2.all()


def all():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    import pickle
    

    outfile = "fiboPXprob_v2.pdf"
    
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

    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6,5.8))
    
    inax00 = ax[0,0].inset_axes([0.15, 0.09, 0.47, 0.42])
    inax01 = ax[0,1].inset_axes([0.15, 0.09, 0.47, 0.42])
    inax10 = ax[1,0].inset_axes([0.15, 0.09, 0.47, 0.42])
    inax11 = ax[1,1].inset_axes([0.15, 0.09, 0.47, 0.42])


    
    MS = 2
    LW = 1
    LW2= 1.3

    #---------
           
    iax = ax[0][0]
    k1=42; k2=34; k3=26
    iax.plot(x0/k1, np.log(P0[:, 50-k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW, label= 'k = {}'.format(k1));
    iax.plot(x0/k2, np.log(P0[:, 50-k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW, label= 'k = {}'.format(k2));
    iax.plot(x0/k3, np.log(P0[:, 50-k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2, label= 'k = {}'.format(k3));

    iax = inax00
    iax.plot(x0/k1, np.log(P0[:, 50-k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x0/k2, np.log(P0[:, 50-k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x0/k3, np.log(P0[:, 50-k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2);
    iax.plot(x0/k3, 2*x0/k3, 'k-',  mfc='none', ms=MS,  lw=LW/2);


    #---------
    
    iax = ax[0][1]
    k1=46; k2=38; k3=30
    iax.plot(x1/k1, np.log(P1[:, 10+k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW, label='k = {}'.format(k1));
    iax.plot(x1/k2, np.log(P1[:, 10+k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW, label='k = {}'.format(k2));
    iax.plot(x1/k3, np.log(P1[:, 10+k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2, label='k = {}'.format(k3));

    iax = inax01
    iax.plot(x1/k1, np.log(P1[:, 10+k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x1/k2, np.log(P1[:, 10+k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x1/k3, np.log(P1[:, 10+k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2);
    iax.plot(x1/k3, 2*x0/k3, 'k-',  mfc='none', ms=MS,  lw=LW/2);

    #---------
    
    iax = ax[1][0]
    k1=46; k2=38; k3=30    
    iax.plot(x14/k1, np.log(P14[:, 50-k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW, label= 'k = {}'.format(k1));
    iax.plot(x14/k2, np.log(P14[:, 50-k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW, label= 'k = {}'.format(k2));
    iax.plot(x14/k3, np.log(P14[:, 50-k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2, label= 'k = {}'.format(k3));

    iax = inax10
    iax.plot(x14/k1, np.log(P14[:, 50-k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x14/k2, np.log(P14[:, 50-k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x14/k3, np.log(P14[:, 50-k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2);
    iax.plot(x14/k3, 2*x0/k3, 'k-',  mfc='none', ms=MS,  lw=LW/2);

    
    #---------
    
    iax = ax[1][1]
    k1=46; k2=38; k3=30
    iax.plot(x34/k1, np.log(P34[:, 10+k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW, label='k = {}'.format(k1));
    iax.plot(x34/k2, np.log(P34[:, 10+k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW, label='k = {}'.format(k2));
    iax.plot(x34/k3, np.log(P34[:, 10+k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2, label='k = {}'.format(k3));

    iax = inax11
    iax.plot(x34/k1, np.log(P34[:, 10+k1-1])/k1, 'r-',  mfc='none', ms=MS,  lw=LW);
    iax.plot(x34/k2, np.log(P34[:, 10+k2-1])/k2, 'b--', mfc='none', ms=MS,  lw=LW);
    iax.plot(x34/k3, np.log(P34[:, 10+k3-1])/k3, 'g:',  mfc='none', ms=MS,  lw=LW2);
    iax.plot(x34/k3, 2*x0/k3, 'k-',  mfc='none', ms=MS,  lw=LW/2);



    
    #-- adjustments ----------------------------------------------------------


    HL=1.3
    BP=0.1
    
    iax = ax[0][0]
    iax.set_xlim(0,0.08)
    iax.set_ylim(-0.5,0)
    iax.set(xlabel='$X_k / k$', ylabel='$k^{-1} \ln(P(X_k))$')
    iax.text(0.037, -0.035, "$\\alpha = 0$")
    iax.yaxis.labelpad = 0  
    iax.legend(frameon=False, handlelength=HL, borderpad=BP)
 

    iax = ax[0][1]
    iax.set_xlim(0,0.08)
    iax.set_ylim(-0.5,0)
    ticks = [-0.4, -0.3, -0.2, -0.1, 0]
    iax.set_yticks(ticks)
    iax.set(xlabel='$X_k / k$', ylabel='$k^{-1} \ln(P(X_k))$')
    iax.text(0.037, -0.035, "$\\alpha = 1$")
    iax.yaxis.labelpad = 0  
    iax.legend(frameon=False, handlelength=HL, borderpad=BP)

           
    iax = ax[1][0]
    iax.set_xlim(0,0.05)
    iax.set_ylim(-0.4,0)
    ticks = [-0.4, -0.3, -0.2, -0.1, 0]
    iax.set_yticks(ticks)
    iax.set(xlabel='$X_k / k$', ylabel='$k^{-1} \ln(P(X_k))$')
    iax.text(0.02, -0.028, "$\\alpha = 1/4$")
    iax.yaxis.labelpad = 0  
    iax.legend(frameon=False, loc='upper right', handlelength=HL, borderpad=BP)

    iax = ax[1][1]
    iax.set_xlim(0,0.05)
    iax.set_ylim(-0.4,0)
    ticks = [-0.4, -0.3, -0.2, -0.1, 0]
    iax.set_yticks(ticks)
    iax.set(xlabel='$X_k / k$', ylabel='$k^{-1} \ln(P(X_k))$')
    iax.text(0.02, -0.028, "$\\alpha = 3/4$")
    iax.yaxis.labelpad = 0  
    iax.legend(frameon=False, loc='upper right', handlelength=HL, borderpad=BP)



    #---------------------

    for iax in (inax00, inax10, inax01, inax11):
    
        iax.set_xlim(-0.25, 0.04)
        ticks=[-0.2, -0.1,  0]
        iax.set_xticks(ticks)
        iax.set_xticklabels(["$%g$"%y for y in ticks], fontsize=8, position=(0,0.03))   
        iax.set_ylim(-0.5, 0)
        ticks=[-0.4, -0.2, 0]
        iax.set_yticks(ticks)
        iax.set_yticklabels(["$%g$"%y for y in ticks], fontsize=8, position=(0.03,0))
        iax.axvline(x=0,color='gray', lw = 0.5)
  


    
    #--------------------------------------
    
 
    plt.subplots_adjust(left=0.10, right=0.98, top=0.98, bottom=0.07, hspace=0.25, wspace=0.3)

    plt.savefig(outfile, pad_inches=0)

    plt.close()
 


#=========================================================
