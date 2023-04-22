
#import importlib ; import fiboSigmaZeroFj
#importlib.reload(fiboSigmaZeroFj) ; fiboSigmaZeroFj.all()


def all():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    import pickle
    

    outfile = "fiboSigmaZeroFj.pdf"

    jmax = 11

    #fdir = '../PX/'
    fdir = '../POSTFrisch/'

    
    with open(fdir + 'q7_v0dt1_i12_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Ca1 = dat["Corr"]

    with open(fdir + 'q7_v0dt1_i24_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Ca2 = dat["Corr"]
    
    with open(fdir + 'q7_v0dt1_i36_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Ca3 = dat["Corr"]


    with open(fdir + 'q7b_g3500dt5_i12_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cb1 = dat["Corr"]

    with open(fdir + 'q7b_g3500dt5_i24_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cb2 = dat["Corr"]
    
    with open(fdir + 'q7b_g3500dt5_i36_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cb3 = dat["Corr"]
    

    fdir = '../POSTFrisch/'

    
    with open(fdir + 'q7_v14L10_i12_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cc1 = dat["Corr"]

    with open(fdir + 'q7_v14L10_i24_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cc2 = dat["Corr"]

    with open(fdir + 'q7_v14L10_i36_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cc3 = dat["Corr"]



    with open(fdir + 'q7_v34R28_i12_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cd1 = dat["Corr"]

    with open(fdir + 'q7_v34R28_i24_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cd2 = dat["Corr"]

    with open(fdir + 'q7_v34R28_i36_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cd3 = dat["Corr"]


    
    #----------------------------

    fig, ax = plt.subplots(ncols=4, nrows=2, figsize=(8,4))

    MS = 4
    LW = 1

    def f(j):
        return((-1)**j)

    j = np.arange(jmax)

    Fj = Fibonacci(jmax)
    Fj2 = Fj*Fj


    #-- plot data

    #-- alpha = 0

    
    iax = ax[0][0]
  
    iax.plot(j, Ca1, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Ca2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Ca3, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);

    
    iax = ax[1][0]
  
    iax.plot(j, Ca1*Fj2, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Ca2*Fj2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Ca3*Fj2, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);

    
    #-- alpha = 1/4

    
    iax = ax[0][1]
  
    iax.plot(j, Cc1, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cc2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cc3, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);

    
    iax = ax[1][1]
  
    iax.plot(j, Cc1*Fj2, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cc2*Fj2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cc3*Fj2, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);


    #-- alpha = 3/4

    
    iax = ax[0][2]
  
    iax.plot(j, Cd1, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cd2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cd3, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);

    
    iax = ax[1][2]
  
    iax.plot(j, Cd1*Fj2, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cd2*Fj2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cd3*Fj2, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);



   
    
    #-- alpha = 1

    
    iax = ax[0][3]
 
    iax.plot(j, Cb1, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cb2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cb3, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);

    
    iax = ax[1][3]
 
    iax.plot(j, Cb1*Fj2, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cb2*Fj2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cb3*Fj2, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);


    

    #-- set canvas --



    for c in (0,1,2,3):

        iax = ax[0][c]

        iax.set_xlim(0,10)
        ticks=[0, 2, 4, 6, 8, 10]
        iax.set_xticks(ticks)
        iax.set_xticklabels(["$%g$"%y for y in ticks]) # position=(0,0.05))
        iax.set_ylim(-0.5,1)
        ticks=[-0.5, 0, 0.5, 1]
        iax.set_yticks(ticks)
        iax.set_yticklabels(["$%g$"%y for y in ticks])
        iax.set(xlabel='$j$', ylabel='$C_j$')
        iax.grid(axis="x", lw=LW/2)
        iax.xaxis.labelpad = -2
        iax.yaxis.labelpad = -10


        iax = ax[1][c]

        iax.set_xlim(0,10)
        ticks=[0, 2, 4, 6, 8, 10]
        iax.set_xticks(ticks)
        iax.set_xticklabels(["$%g$"%y for y in ticks]) # position=(0,0.05))
        iax.set_ylim(-2,2)
        iax.set(xlabel='$j$', ylabel='$C_j  F^2_{j+1}$')
        iax.grid(axis="x", lw=LW/2)
        iax.xaxis.labelpad = -2
        iax.yaxis.labelpad = -4

        
    ax[0][0].legend(frameon=True)

    ax[0][0].text(4.0, 1.05, '$\\alpha=0$')
    ax[0][1].text(3.8, 1.05, '$\\alpha=1/4$')
    ax[0][2].text(3.8, 1.05, '$\\alpha=3/4$')
    ax[0][3].text(4.0, 1.05, '$\\alpha=1$')



    
    #------

    plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.08, hspace=0.28, wspace=0.34)
       
    plt.savefig(outfile, pad_inches=0)

    plt.close()
 




#===========================================



#=========================================================

def Fibonacci(n):

    import numpy as np
     
    Fi = np.zeros(n, 'uint64')
     
    Fi[0]=1
    Fi[1]=1

    for i in range(2,n):
        
        Fi[i] = Fi[i-1] + Fi[i-2]
        
    return(Fi)
  
#=========================================================


