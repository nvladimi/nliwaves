#import importlib ; import plotCap512spectrum
#importlib.reload(plotCap512spectrum) ; plotCap512spectrum.Plot()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np


#==============================================================================

def Plot():

    fnameout = "a1n512_spectra.pdf"
    a = np.load('../wexac/a1n512_20210123/n512spc_kxky.npy')
    b = np.load('../wexac/a1n512_20210123/n512spc_kphi.npy')

    print("min, max: ", np.min(a), np.max(a))
    print("min, max: ", np.min(b), np.max(b))
    vmax=2e-3

    N = a.shape[0]
    (Nk, Nphi) = b.shape
    b = np.fft.fftshift(b,1)
    b0 = b[:,0].reshape(Nk,1)
    b = np.hstack((b, b0))
     
    #-- images --
    
    fig, ax = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(8,8))

    iax = ax[0,0]
    extent=(-N/2, N/2,-N/2, N/2)
    iax.imshow(a, cmap='jet', origin = 'lower', vmin=0, vmax=vmax, extent=extent);
    iax.set_title("averaged spectrum, forcing 0.001")
    iax.set_xlabel('$k_x$')
    iax.set_ylabel('$k_y$')
    
    iax = ax[0,1]
    extent=(-1, 1, 1, Nk)
    iax.imshow(b, cmap='jet', origin = 'lower', vmin=0, vmax=vmax, extent=extent, aspect='auto');
    iax.set_title("averaged shifted spectrum, forcing 0.001")
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$|k|$')

        
    iax = ax[1,0]
    phi = np.arange(-Nphi/2, Nphi/2+1)/Nphi*2 
    for k in (4,8,16,32,64):
        lbl = '$k = {}$'.format(k)
        iax.plot(phi, b[k-1,:],label=lbl)
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$k^4 |a_k|^2$')
    iax.set_xlim((-1,1))
    iax.set_ylim((0,0.0025))
    iax.set_title("angular slices, forcing 0.001")
    iax.legend(frameon=False)

    iax = ax[1,1]
    k = np.arange(Nk)+1 
    for i in (np.arange(4,9)/8*Nphi).astype(int) :
        lbl = "${}\pi$".format(i/Nphi-0.5)
        iax.plot(k, b[:,i], label = lbl)
    iax.plot(k, np.average(b,1), '--k', label = "avg")
    iax.set_xlabel('$k$')
    iax.set_ylabel('$k^4 |a_k|^2$')
    iax.set_xlim((0,N/2))
    iax.set_ylim((0,0.0025))
    iax.set_title("radial slices, forcing 0.001")
    iax.legend(frameon=False)
        
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()
    


#==============================================================================




#===============================================================================



#===============================================================================






