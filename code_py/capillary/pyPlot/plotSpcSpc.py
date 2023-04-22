#import importlib ; import plotSpcSpc
#importlib.reload(plotSpcSpc) ; plotSpcSpc.Plot()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import nlwTools as nlw

import numpy as np


#==============================================================================

def Plot(run="a5"):

    fnameout = "plotSpcSpc_" + run + ".pdf"

    prefix  = "../angular/avgdata/"
    
    fbase = {'a1': "a1n256",
             'a2': "a2n256",
             'a5': "a5n256",
             'n512': "n512spc",
             'hse1': "a1_256_spc",
             'hse2': "a2_256_spc"}
             
    vmax = {'a1': 2e-3,
            'a2': 4e-3,
            'a5': 1.2e-2,
            'n512': 2e-3,
            'hse1': 2e-3,
            'hse2': 4e-3
    }

    
    ymin = {'a1': 0.0004,
            'a2': 0.0010,
            'a5': 0.0030}

    ymax = {'a1': 0.0016,
            'a2': 0.0030,
            'a5': 0.0090}

    
    label = {'a1': 'forcing = 0.001',
             'a2': 'forcing = 0.002',
             'a5': 'forcing = 0.005',
             'n512': 'forcing = 0.001, n=512',
             'hse1': 'forcing = 0.001, HSE',
             'hse2': 'forcing = 0.002, HSE'}

    #-- data --

    a = np.load(prefix + fbase[run] + "_kxky.npy")
    b = np.load(prefix + fbase[run] + "_kphi.npy")

    a = nlw.nkInterpolate(a, dNk=1, Nphi=128, AngleAvg = False, iShift = True)
    k = np.arange(1, b.shape[0]+1)

    print (a.shape, b.shape, k[0], k[-1])
 
    (Nk, Nphi) = b.shape
 
    a = np.fft.fftshift(a,1)
    b = np.fft.fftshift(b,1)

    fa = np.fft.fft(a,axis=1)/a.shape[1]
    fb = np.fft.fft(b,axis=1)/b.shape[1]

    fa = np.sqrt(np.real(fa * np.conjugate(fa)))
    fb = np.sqrt(np.real(fb * np.conjugate(fb)))

    fa = fa[:, 0:int(0.5*fa.shape[1])]
    fb = fb[:, 0:int(0.5*fb.shape[1])]
    
    print (a.shape, b.shape, fa.shape, fb.shape)

    
    #-- images --
    
    fig, ax = plt.subplots(ncols=2, nrows=3, constrained_layout=True, figsize=(7,9))

    iax = ax[0,0]
    extent=(-1, 1, 1, Nk)
    iax.imshow(a, cmap='jet', origin = 'lower', vmin=0, vmax=vmax[run], extent=extent, aspect='auto');
    iax.set_title("averaged spectrum, "+ label[run])
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$|k|$')

    iax = ax[0,1]
    extent=(-1, 1, 1, Nk)
    iax.imshow(b, cmap='jet', origin = 'lower', vmin=0, vmax=vmax[run], extent=extent, aspect='auto');
    iax.set_title("averaged shifted spectrum")
    iax.set_xlabel('$\phi / \pi$')
    iax.set_ylabel('$|k|$')

    iax = ax[1,0]
    extent=(0, fa.shape[1], 1, fa.shape[0])
    iax.imshow(np.log(fa), cmap='jet', origin = 'lower', extent=extent, aspect='auto');
    iax.set_title("Fourier of averaged spectrum")
    iax.set_xlabel('$\phi_m$')
    iax.set_ylabel('$|k|$')

    iax = ax[1,1]
    extent=(0, fb.shape[1], 1, fb.shape[0])
    iax.imshow(np.log(fb), cmap='jet', origin = 'lower', extent=extent, aspect='auto');
    iax.set_title("Fourier of shifted spectrum")
    iax.set_xlabel('$\phi_m$')
    iax.set_ylabel('$|k|$')

    
    iax = ax[2,0]
    k = np.arange(fa.shape[0])+1
    for m in (0,1,2,3):
        lbl = '$m = {}$'.format(m)
        iax.loglog(k, fa[:,m], '-o', linewidth=1, markersize=2, label=lbl)
    iax.loglog(k, vmax[run]*np.power(1.*k,-1), '--k', linewidth = 1) # label='$k^{-1}$')
    iax.loglog(k, vmax[run]*np.power(1.*k,-2), '--k', linewidth = 1) # label='$k^{-2}$')
    iax.set_ylim((1e-6, 1e-2))
    iax.set_xlabel('$|k|$')
    iax.set_ylabel('$|F(|a_k|^2 k^4)_m|$')
    iax.set_title("Angular modes of averaged spectrum")
    iax.legend(frameon=False)
    #iax.grid("on")


    iax = ax[2,1]
    k = np.arange(fb.shape[0])+1
    for m in (0,1,2,3):
        lbl = '$m = {}$'.format(m)
        iax.loglog(k, fb[:,m], '-o', linewidth=1, markersize=2, label=lbl)
    iax.loglog(k, vmax[run]*np.power(1.*k,-1), '--k', linewidth = 1)
    iax.loglog(k, vmax[run]*np.power(1.*k,-2), '--k', linewidth = 1)
    iax.set_ylim((1e-6, 1e-2))
    iax.set_xlabel('$|k|$')
    iax.set_ylabel('$|F(|a_k|^2 k^4)_m|$')
    iax.set_title("Angular modes of shifted spectrum")
    iax.legend(frameon=False)
    #iax.grid("on")
 

    
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()
    


#==============================================================================




#===============================================================================



#===============================================================================






