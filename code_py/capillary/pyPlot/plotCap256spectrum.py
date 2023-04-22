#import importlib ; import plotCap256spectrum
#importlib.reload(plotCap256spectrum) ; plotCap256spectrum.Plot()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import numpy as np


#==============================================================================

def Plot(run="a1"):
    
    prefix = "../2020-01-06_klo/postdata/"
    suffix = "n256"
    
    a = np.load(prefix + run + suffix + "_kxky.npy")
    b = np.load(prefix + run + suffix + "_kphi.npy")

    print("min, max: ", np.min(a), np.max(a))
    print("min, max: ", np.min(b), np.max(b))

    vmax = {'a1': 2e-3,
            'a2': 4e-3,
            'a5': 1.2e-2}

    ymin = {'a1': 0.0004,
            'a2': 0.0010,
            'a5': 0.0030}

    ymax = {'a1': 0.0016,
            'a2': 0.0030,
            'a5': 0.0090}

    
    label = {'a1': 'forcing = 0.001',
             'a2': 'forcing = 0.002',
             'a5': 'forcing = 0.005'}


    fnameout = {'a1': 'a1n256_spectra.pdf',
                'a2': 'a2n256_spectra.pdf',
                'a5': 'a5n256_spectra.pdf'}


    c = 1
    a=a*c
    b=b*c
    
    N = a.shape[0]
    (Nk, Nphi) = b.shape
    b = np.fft.fftshift(b,1)
    b0 = b[:,0].reshape(Nk,1)
    b = np.hstack((b, b0))
     
    #-- images --
    
    fig, ax = plt.subplots(ncols=2, nrows=2, constrained_layout=True, figsize=(8,8))

    iax = ax[0,0]
    extent=(-N/2, N/2,-N/2, N/2)
    iax.imshow(a, cmap='jet', origin = 'lower', vmin=0, vmax=vmax[run]*c, extent=extent);
    iax.set_title("averaged spectrum, "+ label[run])
    iax.set_xlabel('$k_x$')
    iax.set_ylabel('$k_y$')
    
    iax = ax[0,1]
    extent=(-1, 1, 1, Nk)
    iax.imshow(b, cmap='jet', origin = 'lower', vmin=0, vmax=vmax[run]*c, extent=extent, aspect='auto');
    iax.set_title("averaged shifted spectrum")
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$|k|$')

        
    iax = ax[1,0]
    phi = np.arange(-Nphi/2, Nphi/2+1)/Nphi*2 
    for k in (8,16,24,32):
        lbl = '$k = {}$'.format(k)
        iax.plot(phi, b[k-1,:],label=lbl)
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$k^4 |a_k|^2$')
    iax.set_xlim((-1,1))
    iax.set_ylim((ymin[run]*c, ymax[run]*c))
    iax.set_title("angular slices, " + label[run])
    iax.legend(frameon=False)
    iax.grid("on")
    
    iax = ax[1,1]
    k = np.arange(Nk)+1 
    for i in (np.arange(4,9)/8*Nphi).astype(int) :
        lbl = "${}\pi$".format(i/Nphi-0.5)
        iax.plot(k, b[:,i], label = lbl)
    iax.plot(k, np.average(b,1), '--k', label = 'avg')
    
    iax.set_xlabel('$k$')
    iax.set_ylabel('$k^4 |a_k|^2$')
    iax.set_xlim((0,N/2))
    iax.set_ylim((ymin[run]*c, ymax[run]*c))
    iax.set_title("radial slices, " + label[run])
    iax.legend(frameon=False)
    iax.grid("on")


    
    plt.savefig(fnameout[run], pad_inches=0)

    plt.show()
    plt.close()
    


#==============================================================================




#===============================================================================



#===============================================================================






