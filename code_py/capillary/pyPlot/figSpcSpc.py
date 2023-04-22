#import importlib ; import figSpcSpc
#importlib.reload(figSpcSpc) ; figSpcSpc.Plot()

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import nlwTools as nlw

import numpy as np


#==============================================================================

def Plot():

    fnameout = "figSpcSpc.pdf"
    
    prefix  = "../angular/avgdata/"

    fbase = {'a1': "a1n256",
             'a2': "a2n256",
             'a5': "a5n256",
             'n512': "n512spc",
             'hse1': "a1_256_spc",
             'hse2': "a2_256_spc"}

    colors = ("green", "red", "blue")
              
    #-- data --

    a = np.load(prefix + fbase["a2"] + "_kxky.npy")
    b = np.load(prefix + fbase["a2"] + "_kphi.npy")
    c = np.load(prefix + fbase["hse2"] + "_kphi_new.npy")

    a = nlw.nkInterpolate(a, dNk=1, Nphi=128, AngleAvg = False, iShift = True)

    #print (a.shape, b.shape, k[0], k[-1])

    fa = np.fft.fft(a,axis=1)/a.shape[1]
    fb = np.fft.fft(b,axis=1)/b.shape[1]
    fc = np.fft.fft(c,axis=1)/c.shape[1]

    fa = np.sqrt(np.real(fa * np.conjugate(fa)))
    fb = np.sqrt(np.real(fb * np.conjugate(fb)))
    fc = np.sqrt(np.real(fc * np.conjugate(fc)))

    fa = fa[:, 0:int(0.5*fa.shape[1])]
    fb = fb[:, 0:int(0.5*fb.shape[1])]
    fc = fc[:, 0:int(0.5*fc.shape[1])]

    #-- plot --


    #-- all plots --
    
    fig, iax = plt.subplots(ncols=1, nrows=1, constrained_layout=True, figsize=(3,2.8))

    k = np.arange(fa.shape[0])+1
    
    for m in (1,2):
        lbl = '$m = {}$'.format(m)
        iax.loglog(k, fc[:,m]/fc[:,0], '-o',  linewidth=1, markersize=4, color=colors[m], label=lbl, mfc='none')
        iax.loglog(k, fa[:,m]/fa[:,0], ':',   linewidth=1.5, markersize=1, color=colors[m])
        iax.loglog(k, fb[:,m]/fb[:,0], '--d', linewidth=1, markersize=3, color=colors[m], mfc='none')

    iax.loglog(k, np.power(1.*k,-1), '--k', linewidth = 1, label='$k^{-1}$')
    iax.loglog(k, np.power(1.*k,-2), '-.k', linewidth = 1, label='$k^{-2}$')
    iax.set_ylim((1e-3, 1))
    iax.set_xlim((1,40))
    iax.set_xlabel('$|k|$')
    iax.set_ylabel('$F_m/F_0$')


    ticks = 2**(np.arange(0,6))
    iax.set_xticks(ticks)
    iax.set_xticklabels(["$%g$"%y for y in ticks])
   
    iax.legend(frameon=False)

 
    plt.savefig(fnameout, pad_inches=0)

    plt.show()
    plt.close()
   


#==============================================================================




#===============================================================================



#===============================================================================






