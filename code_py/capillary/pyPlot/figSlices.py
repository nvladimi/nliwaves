#import importlib ; import figSlices
#importlib.reload(figSlices) ; figSlices.Plot()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np


#==============================================================================

def Plot(run="a5"):

    fnameout = "figSlices_" + run + ".pdf"
    
    prefix  = "../angular/avgdata/"

    fbase = {'a1': "a1n256",
             'a2': "a2n256",
             'a5': "a5n256",
             'n512': "n512spc",
             'hse1': "a1_256_spc",
             'hse2': "a2_256_spc"}
    
    ymin = {'a1': 0.0002,
            'a2': 0.0005,
            'a5': 0.0015,
            'n512': 0.0002,
            'hse1': 0.0002,
            'hse2': 0.0005}
    
    ymax = {'a1': 0.0008,
            'a2': 0.0016,
            'a5': 0.0040,
            'n512': 0.0008,
            'hse1': 0.0008,
            'hse2': 0.0020}
    
    label = {'a1': 'forcing = 0.001',
             'a2': 'forcing = 0.002',
             'a5': 'forcing = 0.005',
             'n512': 'forcing = 0.001, n=512',
             'hse1': 'forcing = 0.001, HSE',
             'hse2': 'forcing = 0.002, HSE'}
    
    #-- data --

    a = np.load(prefix + fbase[run] + "_kxky.npy")
    b = np.load(prefix + fbase[run] + "_kphi.npy")
    a = np.load(prefix + fbase[run] + "_kxky_new.npy")
    b = np.load(prefix + fbase[run] + "_kphi_new.npy")
    

    (Nk, Nphi) = b.shape
    b = np.fft.fftshift(b,1)
    b0 = b[:,0].reshape(Nk,1)
    b = np.hstack((b, b0))

    b = b/2           # <-- to take into account sqrt(2) definition of a_k 
    
    #-- plots --
    
    fig, ax = plt.subplots(ncols=2, nrows=1, tight_layout=True, figsize=(5.5,2.75))
    #constrained_layout=True
    iax = ax[0]
    iax.ticklabel_format(axis='y', style='sci', scilimits =(-2,3)) 
    phi = np.arange(-Nphi/2, Nphi/2+1)/Nphi*2 
    for k in (8,16,24,32):
        lbl = '$k = {}$'.format(k)
        iax.plot(phi, b[k-1,:],label=lbl)
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$n_k k^4$')
    iax.set_xlim((-1,1))
    iax.set_ylim((ymin[run], ymax[run]))
    #iax.set_title("angular slices, " + label[run])
    iax.legend(frameon=False, labelspacing=0.2)
    
    iax = ax[1]
    iax.ticklabel_format(axis='y', style='sci', scilimits =(-2,3)) 
    k = np.arange(Nk)+1 
    for i in (np.arange(4,9)/8*Nphi).astype(int) :
        q = 2*i/Nphi-1
        if q == 0:
            lbl = "0"
        elif q == 0.5:
            lbl = "$\pi/2$"
        elif q == 1:
            lbl = "$\pi$"
        else:
            lbl = None
        iax.plot(k, b[:,i], label = lbl)
    iax.plot(k, np.average(b,1), '--k') #label = 'avg')
    
    iax.set_xlabel('$k$')
    iax.set_ylabel('$n_k k^4$')
    iax.set_xlim((0,Nk+1))
    iax.set_ylim((ymin[run], ymax[run]))
    #iax.set_title("radial slices, " + label[run])
    iax.legend(frameon=False, labelspacing=0.2)


    
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()
    


#==============================================================================




#===============================================================================



#===============================================================================






