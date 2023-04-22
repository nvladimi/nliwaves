
#import importlib ; import fiboSigmaZero
#importlib.reload(fiboSigmaZero) ; fiboSigmaZero.all()


def all():

    import numpy as np
    import matplotlib
    matplotlib.rcParams['text.usetex'] = True
    import matplotlib.pyplot as plt
    import pickle
    

    outfile = "fiboSigmaZero.pdf"

    jmax = 11

    with open('../PX/q7_v0dt1_i12_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Ca1 = dat["Corr"]

    with open('../PX/q7_v0dt1_i24_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Ca2 = dat["Corr"]
    
    with open('../PX/q7_v0dt1_i36_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Ca3 = dat["Corr"]




    with open('../PX/q7b_g3500dt5_i12_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cb1 = dat["Corr"]

    with open('../PX/q7b_g3500dt5_i24_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cb2 = dat["Corr"]
    
    with open('../PX/q7b_g3500dt5_i36_sigmacorrzero.pkl', 'rb') as handle:
        dat = pickle.load(handle)
    Cb3 = dat["Corr"]
    

    #----------------------------

    fig, ax = plt.subplots(ncols=2, nrows=2, tight_layout=True, figsize=(6,5.5))

    MS = 6
    LW = 1

    def f(j):
        return((-1)**j)

    j = np.arange(jmax)
    
    iax = ax[0,0]
 
    iax.set_xlim(0,10)
    iax.set_ylim(-0.5,1)
    iax.set(xlabel='$j$', ylabel='$C_j$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
 
    iax.plot(j, Ca1, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Ca2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Ca3, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);
    iax.legend(frameon=False)


    iax = ax[0,1]
 
    iax.set_xlim(0,10)
    iax.set_ylim(-0.5,1)
    iax.set(xlabel='$j$', ylabel='$C_j$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
 
    iax.plot(j, Cb1, 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cb2, 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cb3, 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, j*0, '-k',  lw=LW/2);
    iax.legend(frameon=False)


    
    iax = ax[1,0]
 
    iax.set_xlim(0,10)
    iax.set_ylim(5e-5,1)
    iax.set(xlabel='$j$', ylabel='$(-1)^{j+1} \,C_j$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
 
    
    iax.plot(j, -Ca1*f(j), 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, -Ca2*f(j), 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, -Ca3*f(j), 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, np.exp(-j), '--k',  lw=LW/2, label="exp(-j)");

    iax.legend(frameon=False)
    iax.set_yscale('log')

    iax = ax[1,1]
 
    iax.set_xlim(0,10)
    iax.set_ylim(5e-5,1)
    iax.set(xlabel='$j$', ylabel='$(-1)^j \, C_j$')
    #iax.text(0.1, 0.9, '$\\alpha=0, \, i = 12$', transform=iax.transAxes)
 
    iax.plot(j, Cb1*f(j), 'gv-',  mfc='none', ms=MS,  lw=LW,  label="$i = 12$");
    iax.plot(j, Cb2*f(j), 'bs-',  mfc='none', ms=MS,  lw=LW,  label="$i = 24$");
    iax.plot(j, Cb3*f(j), 'or-',  mfc='none', ms=MS,  lw=LW,  label="$i = 36$");
    iax.plot(j, np.exp(-j), '--k',  lw=LW/2, label="exp(-j)");
    iax.legend(frameon=False)
    iax.set_yscale('log')
    
    

    #------

        
    plt.tight_layout()
    plt.savefig(outfile, pad_inches=0)

    plt.close()
 




#===========================================



#=========================================================
