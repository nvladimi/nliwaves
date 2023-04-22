#import importlib ; import plotSpcSpcAll
#importlib.reload(plotSpcSpcAll) ; plotSpcSpcAll.Plot()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import nlwTools as nlw

import numpy as np


#==============================================================================

def Plot():

    fnameout = "plotSpcSpcAll.pdf"

        
    #-- all plots --
    
    fig, ax = plt.subplots(ncols=2, nrows=3, constrained_layout=True, figsize=(7,9))
    
    one_plot(ax[0,0], "a1")
    one_plot(ax[1,0], "a2")
    one_plot(ax[2,0], "a5")

    one_plot(ax[0,1], "hse1")
    one_plot(ax[1,1], "hse2")
    one_plot(ax[2,1], "n512")

    
    
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()






  #==============================================================================
  

def one_plot(iax, run):
    
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

    #print (a.shape, b.shape, k[0], k[-1])

    (Nk, Nphi) = b.shape

    a = np.fft.fftshift(a,1)
    b = np.fft.fftshift(b,1)

    fa = np.fft.fft(a,axis=1)/a.shape[1]
    fb = np.fft.fft(b,axis=1)/b.shape[1]

    fa = np.sqrt(np.real(fa * np.conjugate(fa)))
    fb = np.sqrt(np.real(fb * np.conjugate(fb)))

    fa = fa[:, 0:int(0.5*fa.shape[1])]
    fb = fb[:, 0:int(0.5*fb.shape[1])]

    #-- plot --

    k = np.arange(fa.shape[0])+1
    for m in (0,1,2):
        lbl = '$m = {}$'.format(m)
        iax.loglog(k, fb[:,m], '-o', linewidth=1, markersize=2, color='C'+str(m), label=lbl)
        iax.loglog(k, fa[:,m], '--', linewidth=1, markersize=2, color='C'+str(m))

    iax.loglog(k, vmax[run]*np.power(1.*k,-1), '--k', linewidth = 1) # label='$k^{-1}$')
    iax.loglog(k, vmax[run]*np.power(1.*k,-2), '--k', linewidth = 1) # label='$k^{-2}$')
    iax.set_ylim((1e-6, 1e-2))
    iax.set_xlabel('$|k|$')
    iax.set_ylabel('$|F(|a_k|^2 k^4)_m|$')
    iax.set_title(run + ": solid - shifted,  dashed - not")
    iax.legend(frameon=False)

 
    


#==============================================================================




#===============================================================================



#===============================================================================






